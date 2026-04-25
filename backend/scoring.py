"""Transparent evidence-based scoring and tier assignment.

Scores translocation candidates using an explicit linear model whose
components are fully visible in the output JSON.  The score does NOT
rely on p-values — background p is used only as a small penalty factor.

Tier hierarchy (highest to lowest):
  CONFIRMED  — multi-evidence + reciprocal + stringent QC
  VALIDATED  — strong support from 2+ evidence types
  LIKELY     — moderate support, some evidence gaps
  CANDIDATE  — passes basic artifact filters, in top-N
  FILTERED   — rejected for clear artifact reasons (set by filters.py)

Every cluster records ``reject_reasons`` explaining why it was not
promoted to a higher tier, enabling near-miss reporting.
"""
from __future__ import annotations

from typing import Callable, Optional

from models import EvidenceCluster, Tier

# ---------------------------------------------------------------------------
# Evidence weights
# ---------------------------------------------------------------------------
W_SPLIT: float = 5.0
W_DISCORDANT: float = 3.0
W_CLIPPED: float = 2.0
W_RECIPROCAL: float = 4.0

# ---------------------------------------------------------------------------
# Bonuses
# ---------------------------------------------------------------------------
UNIQUE_START_BONUS_PER: float = 2.0
UNIQUE_START_CAP: float = 20.0
MAPQ_BONUS_MAX: float = 10.0

# ---------------------------------------------------------------------------
# Penalties (applied from filter_flags)
# ---------------------------------------------------------------------------
PENALTY_CENTROMERE: float = -30.0
PENALTY_TELOMERE: float = -30.0
PENALTY_BLACKLIST: float = -50.0
PENALTY_ACROCENTRIC: float = -30.0
PENALTY_SEGDUP: float = -20.0
PENALTY_SINGLE_ORIENT: float = -15.0
PENALTY_HIGH_LOCAL_COV: float = -10.0
PENALTY_UNEVEN: float = -5.0

# Severe flags that block CONFIRMED / VALIDATED tiers entirely
_SEVERE_FLAGS = frozenset({"blacklist"})


class ScoringEngine:
    """Score and assign confidence tiers to translocation evidence clusters."""

    def score_and_tier(
        self,
        clusters: list[EvidenceCluster],
        callback: Optional[Callable[[dict], None]] = None,
    ) -> list[EvidenceCluster]:
        """Score every non-filtered cluster and assign tiers.

        Clusters already marked FILTERED by hard filters are left
        unchanged.  All others receive a transparent evidence score
        with full component breakdown and a tier assignment.
        """
        if callback:
            callback({"type": "scan.stage_changed", "stage": "scoring"})

        total = len(clusters)
        for idx, cluster in enumerate(clusters):
            if cluster.tier == Tier.FILTERED:
                continue

            score, components = _compute_score(cluster)
            cluster.score = score
            cluster.score_components = components

            tier, label, reasons = _assign_tier(cluster)
            cluster.tier = tier
            cluster.evidence_label = label
            cluster.reject_reasons = reasons

            if callback and total > 0 and (
                idx % max(1, total // 20) == 0 or idx == total - 1
            ):
                callback({
                    "type": "scan.progress",
                    "stage": "scoring",
                    "pct": round(((idx + 1) / total) * 100, 1),
                })

        clusters.sort(
            key=lambda c: (c.tier != Tier.FILTERED, c.score), reverse=True
        )
        return clusters


# ---------------------------------------------------------------------------
# Transparent scoring
# ---------------------------------------------------------------------------

def _compute_score(cluster: EvidenceCluster) -> tuple[float, dict]:
    """Compute the evidence score with full component breakdown."""
    components: dict = {}

    # --- Evidence counts ---
    components["split"] = round(cluster.split_count * W_SPLIT, 2)
    components["discordant"] = round(cluster.discordant_count * W_DISCORDANT, 2)
    components["clipped"] = round(cluster.clipped_count * W_CLIPPED, 2)
    components["reciprocal"] = round(cluster.reciprocal_support * W_RECIPROCAL, 2)

    # --- Unique start bonus ---
    bonus_a = max(0, cluster.unique_starts_a - 1) * UNIQUE_START_BONUS_PER
    bonus_b = max(0, cluster.unique_starts_b - 1) * UNIQUE_START_BONUS_PER
    components["unique_starts"] = round(
        min(bonus_a + bonus_b, UNIQUE_START_CAP), 2
    )

    # --- MAPQ bonus (linear ramp 20→60 maps to 0→10) ---
    mq = cluster.median_mapq
    if mq >= 60:
        components["mapq_bonus"] = MAPQ_BONUS_MAX
    elif mq >= 20:
        components["mapq_bonus"] = round(MAPQ_BONUS_MAX * (mq - 20) / 40, 2)
    else:
        components["mapq_bonus"] = 0.0

    # --- Penalties from filter flags ---
    penalties: dict[str, float] = {}
    flags = set(cluster.filter_flags)

    _FLAG_PENALTY_MAP = {
        "centromere_proximity": ("centromere", PENALTY_CENTROMERE),
        "telomere_proximity": ("telomere", PENALTY_TELOMERE),
        "blacklist": ("blacklist", PENALTY_BLACKLIST),
        "acrocentric_parm": ("acrocentric", PENALTY_ACROCENTRIC),
        "segdup": ("segdup", PENALTY_SEGDUP),
        "single_orientation": ("single_orientation", PENALTY_SINGLE_ORIENT),
        "high_local_coverage": ("high_local_coverage", PENALTY_HIGH_LOCAL_COV),
        "uneven_support": ("uneven_support", PENALTY_UNEVEN),
    }
    for flag, (key, value) in _FLAG_PENALTY_MAP.items():
        if flag in flags:
            penalties[key] = value

    # --- Background p-value penalty (mild, not primary ranker) ---
    if cluster.background_p > 0.1:
        penalties["high_background"] = -20.0
    elif cluster.background_p > 0.01:
        penalties["high_background"] = -10.0
    elif cluster.background_p > 0.001:
        penalties["high_background"] = -5.0

    components["penalties"] = penalties

    total = (
        components["split"]
        + components["discordant"]
        + components["clipped"]
        + components["reciprocal"]
        + components["unique_starts"]
        + components["mapq_bonus"]
        + sum(penalties.values())
    )
    components["total"] = round(total, 2)

    return total, components


# ---------------------------------------------------------------------------
# Tier assignment with reject-reason tracking
# ---------------------------------------------------------------------------

def _assign_tier(
    cluster: EvidenceCluster,
) -> tuple[Tier, str, list[str]]:
    """Assign a confidence tier and record reasons for non-promotion."""
    reject_reasons: list[str] = []

    has_split = cluster.split_count >= 1
    has_disc = cluster.discordant_count >= 1
    has_clip = cluster.clipped_count >= 1
    n_evidence_types = sum([has_split, has_disc, has_clip])

    # --- Evidence label ---
    if has_split and has_disc:
        label = "multi-evidence"
    elif has_split and has_clip:
        label = "split+clip"
    elif has_disc and has_clip:
        label = "discordant+clip"
    elif has_split:
        label = "split-supported"
    elif has_disc:
        label = "discordant-supported"
    elif has_clip:
        label = "softclip-supported"
    else:
        label = "unknown"

    has_severe = bool(_SEVERE_FLAGS & set(cluster.filter_flags))

    # --- CONFIRMED ---
    # Multi-evidence with stringent quality
    confirmed_ok = True
    if cluster.split_count < 3:
        reject_reasons.append("insufficient_split_support")
        confirmed_ok = False
    if cluster.discordant_count < 5:
        reject_reasons.append("insufficient_discordant_support")
        confirmed_ok = False
    if cluster.median_mapq < 30:
        reject_reasons.append("low_mapq")
        confirmed_ok = False
    if cluster.unique_starts_a < 2 or cluster.unique_starts_b < 2:
        reject_reasons.append("one_sided_support")
        confirmed_ok = False
    if has_severe:
        reject_reasons.append("repeat_overlap")
        confirmed_ok = False
    if "single_orientation" in cluster.filter_flags:
        reject_reasons.append("orientation_inconsistent")
        confirmed_ok = False

    if confirmed_ok:
        return Tier.CONFIRMED, label, []

    # --- VALIDATED ---
    # Strong support from 2+ evidence types
    validated_ok = True
    if n_evidence_types < 2:
        validated_ok = False
    if cluster.total_support < 8:
        validated_ok = False
    if cluster.median_mapq < 25:
        validated_ok = False
    if has_severe:
        validated_ok = False

    if validated_ok:
        return Tier.VALIDATED, label, reject_reasons

    # --- LIKELY ---
    # Moderate support from any evidence type
    if (
        cluster.total_support >= 5
        and cluster.median_mapq >= 20
        and cluster.score > 15
    ):
        return Tier.LIKELY, label, reject_reasons

    # --- CANDIDATE ---
    return Tier.CANDIDATE, label, reject_reasons
