"""Chimerism background model for translocation scoring.

Estimates the genome-wide background rate of chimeric (interchromosomal)
read pairs, then scores each evidence cluster against a Poisson null model
to produce calibrated p-values with Benjamini-Hochberg correction.

Includes bin-pair noise estimation to detect and suppress globally noisy
loci that produce many interchromosomal pairs (e.g. segdups, repeats).
"""
from __future__ import annotations

import logging
from typing import Callable, Optional

import numpy as np
from scipy.stats import poisson

from models import EvidenceCluster

logger = logging.getLogger(__name__)


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction to a list of p-values.

    Args:
        pvalues: Raw p-values, one per hypothesis test.

    Returns:
        Adjusted p-values in the same order as the input. Values are
        capped at 1.0 and monotonicity is enforced.
    """
    n = len(pvalues)
    if n == 0:
        return []

    pv = np.asarray(pvalues, dtype=np.float64)
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=np.float64)
    ranks[order] = np.arange(1, n + 1, dtype=np.float64)

    adjusted = pv * n / ranks

    sorted_adjusted = adjusted[order]
    cummin = sorted_adjusted[-1]
    for i in range(n - 1, -1, -1):
        cummin = min(cummin, sorted_adjusted[i])
        sorted_adjusted[i] = cummin

    result = np.empty(n, dtype=np.float64)
    result[order] = sorted_adjusted

    np.minimum(result, 1.0, out=result)

    return result.tolist()


class BackgroundModel:
    """Estimate chimeric read-pair background and score clusters against it.

    The model divides the genome into fixed-size bins (default 100 kb),
    counts coverage per bin, and uses the genome-wide chimeric rate to
    compute the expected number of chimeric pairs linking any two bins
    under the null hypothesis (random ligation / mapping artifacts).

    Additionally, it estimates per-bin-pair noise levels to detect
    globally noisy loci and suppress them.

    Args:
        bin_size: Genomic bin width in base pairs. Default 100,000 (100 kb).
    """

    def __init__(self, bin_size: int = 100_000) -> None:
        self.bin_size = bin_size
        self.chimeric_rate: float = 0.0
        self._coverage_bins: dict[tuple[str, int], int] = {}
        self._bin_pair_counts: dict[tuple[str, str, int, int], int] = {}
        self._noise_floor: float = 10.0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def estimate_chimeric_rate(
        self,
        discordant_reads: list,
        chrom_lengths: dict[str, int],
        total_reads: int = 0,
    ) -> float:
        """Estimate the genome-wide interchromosomal chimeric rate.

        Uses the number of non-empty coverage bins for normalization
        (instead of total possible bins) to avoid underestimating the
        per-bin-pair expected count.

        Args:
            discordant_reads: List of SVRead objects whose mates map to
                a different chromosome.
            chrom_lengths: Mapping of chromosome name to length in bp.
            total_reads: Total reads processed (for improved rate estimation).

        Returns:
            Estimated chimeric rate as a float.
        """
        total_genome_bp = sum(chrom_lengths.values())
        if total_genome_bp == 0 or len(discordant_reads) == 0:
            self.chimeric_rate = 0.0
            return self.chimeric_rate

        interchrom_count = 0
        for read in discordant_reads:
            mate_chrom = getattr(read, "mate_chrom", None)
            if mate_chrom is not None and mate_chrom != read.chrom:
                interchrom_count += 1

        # Use non-empty bins for normalization — this avoids the problem
        # where total_bins^2 is so large that expected values are tiny,
        # making everything "significant".
        if self._coverage_bins:
            non_empty_bins = max(1, len(self._coverage_bins))
        else:
            non_empty_bins = max(1, total_genome_bp // self.bin_size)

        self.chimeric_rate = interchrom_count / (non_empty_bins * non_empty_bins)

        logger.info(
            "Background chimeric rate estimated: %.3e "
            "(%d interchromosomal pairs, %d non-empty bins, %d total reads)",
            self.chimeric_rate,
            interchrom_count,
            non_empty_bins,
            total_reads,
        )
        return self.chimeric_rate

    def compute_coverage_bins(
        self,
        all_reads: list,
        chrom_lengths: dict[str, int],
    ) -> dict[tuple[str, int], int]:
        """Bin all reads into fixed-size genomic windows and count coverage."""
        bins: dict[tuple[str, int], int] = {}

        for read in all_reads:
            chrom = read.chrom
            pos = read.pos
            chrom_len = chrom_lengths.get(chrom)
            if chrom_len is None:
                continue
            pos = max(0, min(pos, chrom_len - 1))
            bin_idx = pos // self.bin_size
            key = (chrom, bin_idx)
            bins[key] = bins.get(key, 0) + 1

        self._coverage_bins = bins
        logger.info(
            "Coverage bins computed: %d non-empty bins across %d chromosomes",
            len(bins),
            len({k[0] for k in bins}),
        )
        return bins

    def compute_bin_pair_noise(self, discordant_reads: list) -> None:
        """Count reads per (chrA, chrB, binA, binB) tuple for noise estimation.

        Identifies globally noisy loci (e.g. segdup-to-segdup) by computing
        per-bin-pair counts and deriving a noise floor (median + 3*MAD).
        Clusters at noisy loci that don't exceed 2x the local noise level
        are not enriched and will be assigned p=1.0.
        """
        bin_pair_counts: dict[tuple[str, str, int, int], int] = {}
        for read in discordant_reads:
            mate_chrom = getattr(read, "mate_chrom", None)
            mate_pos = getattr(read, "mate_pos", None)
            if mate_chrom is None or mate_pos is None:
                continue
            bin_a = read.pos // self.bin_size
            bin_b = mate_pos // self.bin_size
            # Canonicalize
            if read.chrom <= mate_chrom:
                key = (read.chrom, mate_chrom, bin_a, bin_b)
            else:
                key = (mate_chrom, read.chrom, bin_b, bin_a)
            bin_pair_counts[key] = bin_pair_counts.get(key, 0) + 1

        self._bin_pair_counts = bin_pair_counts

        # Compute noise floor: median + 3*MAD of bin-pair counts
        if bin_pair_counts:
            counts = np.array(list(bin_pair_counts.values()), dtype=np.float64)
            median = float(np.median(counts))
            mad = float(np.median(np.abs(counts - median)))
            self._noise_floor = median + 3 * max(mad, 1.0)
        else:
            self._noise_floor = 10.0

        logger.info(
            "Bin-pair noise: %d non-empty bin pairs, noise floor=%.1f",
            len(bin_pair_counts),
            self._noise_floor,
        )

    def score_clusters(
        self,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict[str, int],
        callback: Optional[Callable] = None,
    ) -> list[EvidenceCluster]:
        """Score each cluster against the Poisson background model.

        For each cluster, the expected chimeric count under the null is:

            expected = coverage(bin_A) * coverage(bin_B) * chimeric_rate
                       / total_non_empty_bins

        Additionally, clusters at globally noisy loci (where the bin-pair
        count exceeds the noise floor and the cluster doesn't show 2x
        enrichment) are not considered significant.

        All p-values are corrected for multiple testing using
        Benjamini-Hochberg FDR control.
        """
        if callback is not None:
            callback({
                "type": "scan.stage_changed",
                "stage": "background_model",
            })

        if not clusters:
            logger.info("No clusters to score against background model.")
            return clusters

        non_empty_bins = max(1, len(self._coverage_bins)) if self._coverage_bins else 1
        total_genome_bins = max(
            1,
            sum(
                max(1, length // self.bin_size)
                for length in chrom_lengths.values()
            ),
        )

        n_clusters = len(clusters)
        raw_pvalues: list[float] = []
        noise_suppressed = 0

        for idx, cluster in enumerate(clusters):
            bin_a_idx = max(0, cluster.pos_a) // self.bin_size
            bin_b_idx = max(0, cluster.pos_b) // self.bin_size
            bin_a = (cluster.chrom_a, bin_a_idx)
            bin_b = (cluster.chrom_b, bin_b_idx)

            cov_a = self._coverage_bins.get(bin_a, 0)
            cov_b = self._coverage_bins.get(bin_b, 0)

            if self.chimeric_rate <= 0 or cov_a == 0 or cov_b == 0:
                expected = 0.0
            else:
                expected = (
                    cov_a * cov_b * self.chimeric_rate / non_empty_bins
                )

            observed = cluster.total_support

            # Local noise check: if this bin-pair is globally noisy and
            # the cluster doesn't show enrichment over local noise, it's
            # not a real signal.
            if self._bin_pair_counts:
                if cluster.chrom_a <= cluster.chrom_b:
                    bp_key = (cluster.chrom_a, cluster.chrom_b, bin_a_idx, bin_b_idx)
                else:
                    bp_key = (cluster.chrom_b, cluster.chrom_a, bin_b_idx, bin_a_idx)
                local_noise = self._bin_pair_counts.get(bp_key, 0)
                if local_noise >= self._noise_floor and observed <= local_noise * 2:
                    # Not enriched over local noise — suppress
                    raw_pvalues.append(1.0)
                    noise_suppressed += 1
                    continue

            # Poisson p-value: P(X >= observed)
            if expected <= 0:
                pval = 0.0 if observed > 0 else 1.0
            else:
                pval = float(poisson.sf(observed - 1, expected))

            raw_pvalues.append(pval)

            if callback is not None and (
                idx % max(1, n_clusters // 20) == 0
                or idx == n_clusters - 1
            ):
                pct = round((idx + 1) / n_clusters * 100, 1)
                callback({
                    "type": "scan.progress",
                    "stage": "background_model",
                    "pct": pct,
                    "detail": (
                        f"Scored {idx + 1}/{n_clusters} clusters "
                        f"(latest expected={expected:.2e}, "
                        f"observed={observed})"
                    ),
                })

        # Apply Benjamini-Hochberg correction
        adjusted_pvalues = benjamini_hochberg(raw_pvalues)

        for cluster, adj_p in zip(clusters, adjusted_pvalues):
            cluster.background_p = adj_p

        sig_count = sum(1 for p in adjusted_pvalues if p < 0.05)
        logger.info(
            "Background scoring complete: %d/%d clusters significant "
            "(BH-adjusted p < 0.05), %d noise-suppressed",
            sig_count,
            n_clusters,
            noise_suppressed,
        )

        return clusters
