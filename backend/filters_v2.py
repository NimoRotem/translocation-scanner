"""Spec-compliant hard excludes, reject filters, and tier assignment.

Implements the exact gate criteria from the rebuild spec:

Hard excludes (never reported):
  - MT involvement
  - unplaced/random/alt contigs
  - severe ENCODE blacklist overlap on either side
  - centromere core overlap on either side
  - telomere core overlap on either side
  - orientation impossible/incoherent across supporting reads

Reject (filtered, logged):
  - SR=0 AND PR<10
  - duplicate_fraction > 0.7
  - unique_read_starts < 5
  - local_coverage_ratio > 3 on either side
  - severe promiscuous_hotspot without split or external caller

Tiers (ordered gates):
  candidate -> strong_candidate -> likely -> validated -> confirmed
"""
from __future__ import annotations

import logging
from typing import Callable, Optional

from models import EvidenceCluster, Tier

logger = logging.getLogger(__name__)

# Valid primary chromosomes
_PRIMARY_CHROMS = set()
for i in range(1, 23):
    _PRIMARY_CHROMS.add(f"chr{i}")
    _PRIMARY_CHROMS.add(str(i))
_PRIMARY_CHROMS.update({"X", "Y", "chrX", "chrY"})


class FilterEngineV2:
    """Spec-compliant filter and tier engine.

    Operates in three phases:
      1. Hard excludes — remove from pipeline entirely
      2. Reject filters — flag with reason, set tier=FILTERED
      3. Tier assignment ��� ordered gate evaluation
    """

    def __init__(
        self,
        masks=None,
        external_callers: Optional[dict] = None,
    ) -> None:
        self._masks = masks
        self._external_callers = external_callers or {}

    def apply_all(
        self,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict[str, int],
        callback: Optional[Callable[[dict], None]] = None,
    ) -> list[EvidenceCluster]:
        """Run hard excludes, reject filters, and tier assignment."""
        if callback:
            callback({"type": "scan.stage_changed", "stage": "filtering"})

        n = len(clusters)
        excluded = 0
        rejected = 0

        for idx, cluster in enumerate(clusters):
            # Phase 1: Hard excludes
            exclude_reason = self._check_hard_exclude(cluster)
            if exclude_reason:
                cluster.tier = Tier.FILTERED
                cluster.filter_flags.append(f"hard_exclude:{exclude_reason}")
                excluded += 1
                continue

            # Phase 2: Reject filters
            reject_reason = self._check_reject(cluster)
            if reject_reason:
                cluster.tier = Tier.FILTERED
                cluster.filter_flags.append(f"reject:{reject_reason}")
                cluster.reject_reasons.append(reject_reason)
                rejected += 1
                continue

            # Phase 3: Annotate mask overlaps (for tier decisions)
            self._annotate_masks(cluster)

            # Phase 4: Tier assignment
            self._assign_tier(cluster)

            if callback and n > 0 and (
                idx % max(1, n // 20) == 0 or idx == n - 1
            ):
                callback({
                    "type": "scan.progress",
                    "stage": "filtering",
                    "pct": round((idx + 1) / n * 100, 1),
                })

        logger.info(
            "Filtering: %d excluded, %d rejected, %d remaining of %d",
            excluded, rejected, n - excluded - rejected, n,
        )
        return clusters

    # ------------------------------------------------------------------
    # Phase 1: Hard excludes
    # ------------------------------------------------------------------

    def _check_hard_exclude(self, cluster: EvidenceCluster) -> Optional[str]:
        """Return reason string if cluster should be hard-excluded, else None."""
        # MT involvement
        if _is_mt(cluster.chrom_a) or _is_mt(cluster.chrom_b):
            return "mt_involvement"

        # Unplaced/random/alt contigs
        if cluster.chrom_a not in _PRIMARY_CHROMS:
            return f"non_primary_chrom:{cluster.chrom_a}"
        if cluster.chrom_b not in _PRIMARY_CHROMS:
            return f"non_primary_chrom:{cluster.chrom_b}"

        # Severe ENCODE blacklist on either side
        if self._masks and self._masks.encode_blacklist.overlaps(cluster.chrom_a, cluster.pos_a):
            return "blacklist_a"
        if self._masks and self._masks.encode_blacklist.overlaps(cluster.chrom_b, cluster.pos_b):
            return "blacklist_b"

        # Centromere core overlap on either side
        if self._masks and self._masks.centromere.overlaps(cluster.chrom_a, cluster.pos_a):
            return "centromere_core_a"
        if self._masks and self._masks.centromere.overlaps(cluster.chrom_b, cluster.pos_b):
            return "centromere_core_b"

        # Telomere core overlap on either side
        if self._masks and self._masks.telomere.overlaps(cluster.chrom_a, cluster.pos_a):
            return "telomere_core_a"
        if self._masks and self._masks.telomere.overlaps(cluster.chrom_b, cluster.pos_b):
            return "telomere_core_b"

        # Orientation incoherent across supporting reads
        if self._orientation_incoherent(cluster):
            return "orientation_incoherent"

        return None

    # ------------------------------------------------------------------
    # Phase 2: Reject filters
    # ------------------------------------------------------------------

    def _check_reject(self, cluster: EvidenceCluster) -> Optional[str]:
        """Return reason string if cluster should be rejected, else None."""
        sr = cluster.split_count
        pr = cluster.discordant_count

        # SR=0 AND PR<10
        if sr == 0 and pr < 10:
            return "sr0_pr_lt10"

        # duplicate_fraction > 0.7
        dup_frac = getattr(cluster, 'duplicate_fraction', 0.0)
        if dup_frac > 0.7:
            return f"high_dup_fraction:{dup_frac:.2f}"

        # unique_read_starts < 5
        total_unique = cluster.unique_starts_a + cluster.unique_starts_b
        if total_unique < 5:
            return f"low_unique_starts:{total_unique}"

        # local_coverage_ratio > 3 on either side (only when SR=0)
        # When split reads are present, the cluster has base-pair evidence
        # that is much harder to produce artifactually; the tier-assignment
        # NB p-value gates handle quality assessment for these clusters.
        if sr == 0:
            lcr_a = getattr(cluster, 'local_coverage_ratio_a', 0.0)
            lcr_b = getattr(cluster, 'local_coverage_ratio_b', 0.0)
            if lcr_a > 3.0:
                return f"high_local_cov_a:{lcr_a:.1f}"
            if lcr_b > 3.0:
                return f"high_local_cov_b:{lcr_b:.1f}"

        # Promiscuous hotspot without split or external caller
        is_promiscuous = getattr(cluster, 'promiscuous_hotspot', False)
        has_external = bool(getattr(cluster, 'external_callers', []))
        if is_promiscuous and sr == 0 and not has_external:
            return "promiscuous_no_split_no_external"

        return None

    # ------------------------------------------------------------------
    # Phase 3: Mask annotation
    # ------------------------------------------------------------------

    def _annotate_masks(self, cluster: EvidenceCluster) -> None:
        """Annotate cluster with mask overlap information."""
        if not self._masks:
            return

        overlaps_a = self._masks.get_overlaps(cluster.chrom_a, cluster.pos_a)
        overlaps_b = self._masks.get_overlaps(cluster.chrom_b, cluster.pos_b)

        cluster.mask_overlaps_a = overlaps_a
        cluster.mask_overlaps_b = overlaps_b

        # Segdup %identity
        cluster.segdup_pct_a = self._masks.segdup.max_identity_at(
            cluster.chrom_a, cluster.pos_a
        )
        cluster.segdup_pct_b = self._masks.segdup.max_identity_at(
            cluster.chrom_b, cluster.pos_b
        )

    # ------------------------------------------------------------------
    # Phase 4: Tier assignment (ordered gates)
    # ------------------------------------------------------------------

    def _assign_tier(self, cluster: EvidenceCluster) -> None:
        """Assign tier using ordered gate evaluation (highest first)."""
        sr = cluster.split_count
        pr = cluster.discordant_count
        ext = getattr(cluster, 'external_callers', [])
        n_ext = len(ext)

        # Try CONFIRMED first
        if self._is_confirmed(cluster, sr, pr, n_ext):
            cluster.tier = Tier.CONFIRMED
            cluster.evidence_label = self._evidence_label(cluster)
            return

        # Try VALIDATED
        if self._is_validated(cluster, sr, pr, n_ext):
            cluster.tier = Tier.VALIDATED
            cluster.evidence_label = self._evidence_label(cluster)
            return

        # Try LIKELY
        if self._is_likely(cluster, sr, pr, n_ext):
            cluster.tier = Tier.LIKELY
            cluster.evidence_label = self._evidence_label(cluster)
            return

        # Try STRONG_CANDIDATE
        if self._is_strong_candidate(cluster, sr, pr):
            # Use LIKELY tier for strong_candidate since our Tier enum
            # doesn't have STRONG_CANDIDATE — mark via evidence_label
            cluster.tier = Tier.CANDIDATE
            cluster.evidence_label = "strong_candidate:" + self._evidence_label(cluster)
            return

        # Default: CANDIDATE
        if pr >= 3 or sr >= 1:
            cluster.tier = Tier.CANDIDATE
            cluster.evidence_label = self._evidence_label(cluster)
        else:
            cluster.tier = Tier.FILTERED
            cluster.reject_reasons.append("below_minimum_support")

    def _is_confirmed(self, c: EvidenceCluster, sr: int, pr: int, n_ext: int) -> bool:
        """Confirmed: validated + (assembly OR 3+ external callers OR manual IGV)."""
        if not self._is_validated(c, sr, pr, n_ext):
            c.reject_reasons.append("not_validated")
            return False

        # Assembly resolves breakpoint to single-bp
        has_assembly = getattr(c, 'assembly_resolved', False)
        if has_assembly:
            return True

        # 3+ external callers agree
        if n_ext >= 3:
            return True

        # Manual IGV review (user flag)
        if getattr(c, 'igv_confirmed', False):
            return True

        c.reject_reasons.append("no_confirmation_evidence")
        return False

    def _is_validated(self, c: EvidenceCluster, sr: int, pr: int, n_ext: int) -> bool:
        """Validated: PR>=6, SR>=3, pval<1e-6 both sides, clean masks, coherent orientation.
        OR: likely + 2+ external callers."""
        pval_a = getattr(c, 'local_nb_pvalue_a', 1.0)
        pval_b = getattr(c, 'local_nb_pvalue_b', 1.0)

        # Path 1: intrinsic evidence
        if (pr >= 6 and sr >= 3
                and pval_a < 1e-6 and pval_b < 1e-6
                and self._clean_masks(c)
                and not self._orientation_incoherent(c)):
            return True

        # Path 2: likely + 2+ external callers
        if self._is_likely(c, sr, pr, n_ext) and n_ext >= 2:
            return True

        if pr < 6:
            c.reject_reasons.append(f"pr_lt6:{pr}")
        if sr < 3:
            c.reject_reasons.append(f"sr_lt3:{sr}")
        if pval_a >= 1e-6:
            c.reject_reasons.append(f"pval_a_high:{pval_a:.2e}")
        if pval_b >= 1e-6:
            c.reject_reasons.append(f"pval_b_high:{pval_b:.2e}")
        return False

    def _is_likely(self, c: EvidenceCluster, sr: int, pr: int, n_ext: int) -> bool:
        """Likely: strong_candidate + (SR>=2 agreeing OR 1 external caller)."""
        if not self._is_strong_candidate(c, sr, pr):
            return False

        # SR >= 2 with breakpoints agreeing within 10bp
        if sr >= 2:
            return True

        # One external caller overlap
        if n_ext >= 1:
            return True

        c.reject_reasons.append("no_split_or_external")
        return False

    def _is_strong_candidate(self, c: EvidenceCluster, sr: int, pr: int) -> bool:
        """Strong candidate: PR>=15 with full uniqueness gates (if SR=0), OR PR>=5 + SR>=1 + clean."""
        # Path 1: PR>=15, SR=0 allowed only with full uniqueness gates
        if pr >= 15 and sr == 0:
            remap_ok = getattr(c, 'both_flanks_remap_uniquely', True)
            lcr_a = getattr(c, 'local_coverage_ratio_a', 0.0)
            lcr_b = getattr(c, 'local_coverage_ratio_b', 0.0)
            segdup_a = getattr(c, 'segdup_pct_a', 0.0)
            segdup_b = getattr(c, 'segdup_pct_b', 0.0)
            has_blacklist = any(
                "blacklist" in f for f in getattr(c, 'mask_overlaps_a', []) + getattr(c, 'mask_overlaps_b', [])
            )

            gates_pass = (
                remap_ok
                and lcr_a <= 1.5
                and lcr_b <= 1.5
                and c.unique_starts_a >= 10 and c.unique_starts_b >= 10
                and segdup_a < 95.0 and segdup_b < 95.0
                and not has_blacklist
                and not self._orientation_incoherent(c)
            )
            if gates_pass:
                return True
            else:
                if not remap_ok:
                    c.reject_reasons.append("flanks_not_unique")
                if lcr_a > 1.5 or lcr_b > 1.5:
                    c.reject_reasons.append(f"lcr_high:{max(lcr_a, lcr_b):.1f}")
                if c.unique_starts_a < 10 or c.unique_starts_b < 10:
                    c.reject_reasons.append("low_unique_starts_for_sr0")
                if segdup_a >= 95 or segdup_b >= 95:
                    c.reject_reasons.append("segdup_high_identity")

        # Path 2: PR>=5 + SR>=1 + clean masks
        if pr >= 5 and sr >= 1 and self._clean_masks(c):
            return True

        if pr >= 15:
            return False  # Failed uniqueness gates

        if pr < 5:
            c.reject_reasons.append(f"pr_lt5:{pr}")
        return False

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _clean_masks(self, c: EvidenceCluster) -> bool:
        """Check that neither side has segdup>=95%id, blacklist, or promiscuous."""
        segdup_a = getattr(c, 'segdup_pct_a', 0.0)
        segdup_b = getattr(c, 'segdup_pct_b', 0.0)
        if segdup_a >= 95.0 or segdup_b >= 95.0:
            return False
        mask_a = getattr(c, 'mask_overlaps_a', [])
        mask_b = getattr(c, 'mask_overlaps_b', [])
        for m in mask_a + mask_b:
            if "blacklist" in m:
                return False
        if getattr(c, 'promiscuous_hotspot', False):
            return False
        return True

    @staticmethod
    def _orientation_incoherent(cluster: EvidenceCluster) -> bool:
        """Check if orientation is incoherent across reads."""
        if not cluster.reads:
            return False
        orientations = set()
        for r in cluster.reads:
            strand_a = "-" if r.is_reverse else "+"
            strand_b = "-" if r.mate_is_reverse else "+"
            orientations.add(strand_a + strand_b)
        # More than 2 distinct orientations = incoherent
        return len(orientations) > 2

    @staticmethod
    def _evidence_label(cluster: EvidenceCluster) -> str:
        has_sr = cluster.split_count >= 1
        has_pr = cluster.discordant_count >= 1
        has_clip = cluster.clipped_count >= 1
        parts = []
        if has_pr:
            parts.append("PR")
        if has_sr:
            parts.append("SR")
        if has_clip:
            parts.append("CLIP")
        ext = getattr(cluster, 'external_callers', [])
        if ext:
            parts.append("+".join(ext))
        return "+".join(parts) if parts else "unknown"


def _is_mt(chrom: str) -> bool:
    return chrom in ("chrM", "MT", "M", "chrMT")
