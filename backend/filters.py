"""Hard and soft filters for translocation evidence clusters.

Discovery-mode philosophy: most genomic-region filters are SOFT — they
add flags that feed into score penalties but do NOT reject clusters.
Only truly unrecoverable artifacts are hard-filtered:

  Hard filters (tier → FILTERED):
    - total_support < 2
    - median_mapq < 10

  Soft filters (flag added, scoring applies penalty):
    - centromere proximity
    - telomere proximity
    - acrocentric p-arm
    - ENCODE blacklist
    - segmental duplication
    - single orientation
    - high local coverage
    - uneven support
    - high background p-value

GRCh38 coordinates are hardcoded for region-based filtering.
"""
from __future__ import annotations

from typing import Callable, Optional

from models import EvidenceCluster, Tier

# ---------------------------------------------------------------------------
# GRCh38 approximate centromere midpoint positions (bp).
# ---------------------------------------------------------------------------
CENTROMERE_POSITIONS: dict[str, int] = {
    "chr1":  124_000_000,
    "chr2":   93_300_000,
    "chr3":   91_000_000,
    "chr4":   50_400_000,
    "chr5":   48_400_000,
    "chr6":   61_000_000,
    "chr7":   59_900_000,
    "chr8":   45_600_000,
    "chr9":   49_000_000,
    "chr10":  40_200_000,
    "chr11":  53_700_000,
    "chr12":  35_800_000,
    "chr13":  17_900_000,
    "chr14":  17_600_000,
    "chr15":  19_000_000,
    "chr16":  36_600_000,
    "chr17":  24_000_000,
    "chr18":  17_200_000,
    "chr19":  26_500_000,
    "chr20":  27_500_000,
    "chr21":  13_200_000,
    "chr22":  14_700_000,
    "chrX":   60_600_000,
    "chrY":   12_500_000,
}

# ---------------------------------------------------------------------------
# GRCh38 acrocentric chromosome p-arm boundaries.
# ---------------------------------------------------------------------------
ACROCENTRIC_P_ARMS: dict[str, tuple[int, int]] = {
    "chr13": (0, 17_900_000),
    "chr14": (0, 17_600_000),
    "chr15": (0, 19_000_000),
    "chr21": (0, 13_200_000),
    "chr22": (0, 14_700_000),
}

# ---------------------------------------------------------------------------
# Telomere margin
# ---------------------------------------------------------------------------
TELOMERE_MARGIN: int = 500_000  # 500 kb


class FilterEngine:
    """Apply hard and soft filters to translocation evidence clusters.

    Hard filters set the cluster tier to FILTERED.
    Soft filters append descriptive strings to ``filter_flags`` — the
    scoring engine then applies corresponding penalties.
    """

    def __init__(
        self,
        centromere_margin: int = 1_000_000,
        bg_pvalue_threshold: float = 0.001,
        blacklist_index: Optional[object] = None,
        segdup_index: Optional[object] = None,
    ) -> None:
        self.centromere_margin = centromere_margin
        self.bg_pvalue_threshold = bg_pvalue_threshold
        self._blacklist = blacklist_index
        self._segdup = segdup_index

    def apply_filters(
        self,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict[str, int],
        callback: Optional[Callable[[dict], None]] = None,
    ) -> list[EvidenceCluster]:
        """Run all filters over *clusters* in-place."""
        if callback:
            callback({"type": "scan.stage_changed", "stage": "filtering"})

        self._chrom_lengths = chrom_lengths
        total = len(clusters)
        for idx, cluster in enumerate(clusters):
            self._apply_hard_filters(cluster)
            # Always run soft filters — even on filtered clusters, for
            # near-miss reporting.
            self._apply_soft_filters(cluster)

            if callback and total > 0 and (
                idx % max(1, total // 20) == 0 or idx == total - 1
            ):
                callback({
                    "type": "scan.progress",
                    "stage": "filtering",
                    "pct": round(((idx + 1) / total) * 100, 1),
                })

        return clusters

    # ------------------------------------------------------------------
    # Hard filters — only truly unrecoverable artifacts
    # ------------------------------------------------------------------

    def _apply_hard_filters(self, cluster: EvidenceCluster) -> None:
        """Hard-reject only the most obvious artifacts."""
        if cluster.total_support < 2:
            cluster.tier = Tier.FILTERED
            cluster.filter_flags.append("insufficient_support")
            return

        if cluster.median_mapq < 10:
            cluster.tier = Tier.FILTERED
            cluster.filter_flags.append("very_low_mapq")
            return

    # ------------------------------------------------------------------
    # Soft filters — flag for scoring penalties, do NOT reject
    # ------------------------------------------------------------------

    def _apply_soft_filters(self, cluster: EvidenceCluster) -> None:
        """Flag clusters for downstream penalty scoring."""
        # Region-based flags
        if self._either_near_centromere(cluster):
            cluster.filter_flags.append("centromere_proximity")

        if self._either_in_acrocentric_parm(cluster):
            cluster.filter_flags.append("acrocentric_parm")

        if self._either_near_telomere(cluster):
            cluster.filter_flags.append("telomere_proximity")

        if self._in_blacklist(cluster):
            cluster.filter_flags.append("blacklist")

        if self._in_segdup(cluster):
            cluster.filter_flags.append("segdup")

        # Evidence-quality flags
        if self._single_orientation_only(cluster):
            cluster.filter_flags.append("single_orientation")

        if self._high_local_coverage(cluster):
            cluster.filter_flags.append("high_local_coverage")

        if self._uneven_support(cluster):
            cluster.filter_flags.append("uneven_support")

        # Background p-value flag (for scoring penalty)
        if cluster.background_p > self.bg_pvalue_threshold:
            cluster.filter_flags.append("high_background_p")

    # ------------------------------------------------------------------
    # Region predicates
    # ------------------------------------------------------------------

    def _either_near_centromere(self, cluster: EvidenceCluster) -> bool:
        cen_a = CENTROMERE_POSITIONS.get(cluster.chrom_a)
        cen_b = CENTROMERE_POSITIONS.get(cluster.chrom_b)
        if cen_a is not None and abs(cluster.pos_a - cen_a) <= self.centromere_margin:
            return True
        if cen_b is not None and abs(cluster.pos_b - cen_b) <= self.centromere_margin:
            return True
        return False

    @staticmethod
    def _either_in_acrocentric_parm(cluster: EvidenceCluster) -> bool:
        for chrom_attr, pos_attr in [
            (cluster.chrom_a, cluster.pos_a),
            (cluster.chrom_b, cluster.pos_b),
        ]:
            arm = ACROCENTRIC_P_ARMS.get(chrom_attr)
            if arm is not None and arm[0] <= pos_attr <= arm[1]:
                return True
        return False

    def _either_near_telomere(self, cluster: EvidenceCluster) -> bool:
        for chrom_attr, pos_attr in [
            (cluster.chrom_a, cluster.pos_a),
            (cluster.chrom_b, cluster.pos_b),
        ]:
            chrom_len = self._chrom_lengths.get(chrom_attr)
            if chrom_len is None:
                continue
            if pos_attr <= TELOMERE_MARGIN:
                return True
            if pos_attr >= chrom_len - TELOMERE_MARGIN:
                return True
        return False

    def _in_blacklist(self, cluster: EvidenceCluster) -> bool:
        if self._blacklist is None:
            return False
        return (
            self._blacklist.overlaps(cluster.chrom_a, cluster.pos_a)
            or self._blacklist.overlaps(cluster.chrom_b, cluster.pos_b)
        )

    def _in_segdup(self, cluster: EvidenceCluster) -> bool:
        if self._segdup is None:
            return False
        return (
            self._segdup.overlaps(cluster.chrom_a, cluster.pos_a)
            or self._segdup.overlaps(cluster.chrom_b, cluster.pos_b)
        )

    @staticmethod
    def _single_orientation_only(cluster: EvidenceCluster) -> bool:
        if not cluster.reads:
            return False
        strands = {r.is_reverse for r in cluster.reads}
        return len(strands) == 1

    @staticmethod
    def _high_local_coverage(cluster: EvidenceCluster) -> bool:
        return (cluster.discordant_count + cluster.split_count) > 100

    @staticmethod
    def _uneven_support(cluster: EvidenceCluster) -> bool:
        recip = cluster.reciprocal_support
        non_recip = cluster.total_support - recip
        if recip == 0 or non_recip == 0:
            return False
        ratio = max(recip, non_recip) / max(min(recip, non_recip), 1)
        return ratio > 3.0
