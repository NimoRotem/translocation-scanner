"""Local negative binomial background model for translocation scoring.

Replaces the global Poisson model with per-chromosome sliding-window
negative binomial estimation.  Reports local_nb_pvalue per side for
each evidence cluster.

Algorithm:
  1. Compute discordant-pair rate in sliding 10kb windows per chromosome.
  2. Fit negative binomial per chromosome (or per chromosome pair if
     enough data).
  3. For each candidate cluster, report local NB p-value against
     local rate at BOTH breakpoints.
  4. Report as "local_nb_pvalue_a" and "local_nb_pvalue_b" — no combined
     p-value unless combination method is documented.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Callable, Optional

import numpy as np
from scipy.stats import nbinom
from scipy.special import digamma

from models import EvidenceCluster

logger = logging.getLogger(__name__)


def _fit_nb_moments(counts: np.ndarray) -> tuple[float, float]:
    """Fit negative binomial parameters using method of moments.

    Returns (r, p) where r is the dispersion and p is the success probability.
    NB parameterization: mean = r*(1-p)/p, var = r*(1-p)/p^2

    Falls back to Poisson-like (r=very large) if variance <= mean.
    """
    if len(counts) < 3:
        return 1e6, 0.5  # Effectively Poisson

    mean = float(np.mean(counts))
    var = float(np.var(counts))

    if mean <= 0:
        return 1e6, 0.999

    # If variance <= mean, data is under-dispersed; use Poisson approximation
    if var <= mean:
        r = 1e6
        p = r / (r + mean)
        return r, p

    # Method of moments: r = mean^2 / (var - mean), p = mean / var
    r = mean * mean / (var - mean)
    p = mean / var

    # Clamp to reasonable range
    r = max(0.01, min(r, 1e8))
    p = max(1e-10, min(p, 1.0 - 1e-10))

    return r, p


class BackgroundModelV2:
    """Local negative binomial background model.

    Estimates per-chromosome discordant-pair rates using sliding
    windows, fits NB distributions, and scores clusters against
    local background.

    Args:
        window_size: Size of sliding windows for rate estimation (bp).
            Default 10,000 (10kb).
        bin_size: Bin size for counting reads. Default 1,000 (1kb).
    """

    def __init__(
        self,
        window_size: int = 10_000,
        bin_size: int = 1_000,
    ) -> None:
        self.window_size = window_size
        self.bin_size = bin_size
        # Per-chromosome bin counts: chrom -> array of counts per bin
        self._chrom_bins: dict[str, np.ndarray] = {}
        # Per-chromosome NB parameters: chrom -> (r, p)
        self._chrom_nb_params: dict[str, tuple[float, float]] = {}
        # Global chimeric rate (for compatibility)
        self.chimeric_rate: float = 0.0

    def compute_local_rates(
        self,
        discordant_reads: list,
        chrom_lengths: dict[str, int],
        total_reads: int = 0,
    ) -> None:
        """Compute per-chromosome discordant-pair rates in sliding windows.

        Bins all interchromosomal discordant reads by position, then
        fits a NB distribution per chromosome.
        """
        # Count discordant reads per bin
        bin_counts: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        total_interchrom = 0

        for read in discordant_reads:
            mate_chrom = getattr(read, "mate_chrom", None)
            if mate_chrom is not None and mate_chrom != read.chrom:
                total_interchrom += 1
                bin_idx = read.pos // self.bin_size
                bin_counts[read.chrom][bin_idx] += 1
                # Also count the mate side
                mate_pos = getattr(read, "mate_pos", 0) or 0
                mate_bin = mate_pos // self.bin_size
                bin_counts[mate_chrom][mate_bin] += 1

        # Build per-chromosome count arrays
        for chrom, length in chrom_lengths.items():
            n_bins = max(1, length // self.bin_size)
            counts = np.zeros(n_bins, dtype=np.float64)
            for bin_idx, count in bin_counts.get(chrom, {}).items():
                if 0 <= bin_idx < n_bins:
                    counts[bin_idx] = count
            self._chrom_bins[chrom] = counts

            # Fit NB on non-zero bins to get local background parameters
            # Use windowed counts for better NB fit
            window_bins = self.window_size // self.bin_size
            if window_bins > 1 and len(counts) > window_bins:
                # Compute windowed sums
                cumsum = np.cumsum(counts)
                windowed = np.zeros(len(counts), dtype=np.float64)
                for i in range(len(counts)):
                    start = max(0, i - window_bins // 2)
                    end = min(len(counts), i + window_bins // 2)
                    windowed[i] = cumsum[end - 1] - (cumsum[start - 1] if start > 0 else 0)
                r, p = _fit_nb_moments(windowed)
            else:
                r, p = _fit_nb_moments(counts)

            self._chrom_nb_params[chrom] = (r, p)

        # Global chimeric rate for compatibility
        total_bins = sum(max(1, l // self.bin_size) for l in chrom_lengths.values())
        non_empty = sum(1 for c in self._chrom_bins.values() for x in c if x > 0)
        if total_reads > 0:
            self.chimeric_rate = total_interchrom / total_reads
        elif non_empty > 0:
            self.chimeric_rate = total_interchrom / (non_empty * non_empty)

        logger.info(
            "Local background model: %d chromosomes, %d interchrom pairs, "
            "chimeric rate %.4e",
            len(self._chrom_bins), total_interchrom, self.chimeric_rate,
        )

    def local_rate_at(self, chrom: str, pos: int) -> float:
        """Return the local discordant-pair rate in a window around pos."""
        counts = self._chrom_bins.get(chrom)
        if counts is None or len(counts) == 0:
            return 0.0

        center_bin = pos // self.bin_size
        half_window = (self.window_size // self.bin_size) // 2

        start = max(0, center_bin - half_window)
        end = min(len(counts), center_bin + half_window + 1)

        if start >= end:
            return 0.0

        return float(np.sum(counts[start:end]))

    def nb_pvalue_at(self, chrom: str, pos: int, observed: int) -> float:
        """Compute NB p-value for observing >= `observed` at position.

        Uses the per-chromosome NB parameters and local windowed rate.

        Args:
            chrom: Chromosome name.
            pos: Genomic position.
            observed: Number of supporting reads in the cluster at this position.

        Returns:
            P-value (probability of seeing >= observed under NB null).
        """
        r, p = self._chrom_nb_params.get(chrom, (1e6, 0.999))

        local_rate = self.local_rate_at(chrom, pos)

        if local_rate <= 0 and observed > 0:
            return 0.0  # No background, any signal is significant

        if observed <= 0:
            return 1.0

        # Adjust NB parameters to local rate
        # mean = r*(1-p)/p, so for a given local_rate:
        # We scale r to match the local rate while keeping the dispersion ratio
        if r < 1e5:  # Actually dispersed
            # Keep the same dispersion ratio (var/mean), adjust mean
            original_mean = r * (1 - p) / p if p > 0 else local_rate
            if original_mean > 0:
                scale = local_rate / original_mean
                r_local = r * scale
                r_local = max(0.01, r_local)
                p_local = r_local / (r_local + local_rate) if (r_local + local_rate) > 0 else 0.999
            else:
                r_local = max(0.01, local_rate)
                p_local = r_local / (r_local + local_rate) if (r_local + local_rate) > 0 else 0.999
        else:
            # Poisson-like: r → ∞
            r_local = 1e6
            p_local = r_local / (r_local + local_rate) if local_rate > 0 else 0.999

        try:
            # P(X >= observed) = 1 - P(X < observed) = 1 - CDF(observed - 1)
            pval = float(nbinom.sf(observed - 1, r_local, p_local))
        except (ValueError, OverflowError):
            pval = 0.0 if observed > local_rate * 3 else 1.0

        return pval

    def score_clusters(
        self,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict[str, int],
        callback: Optional[Callable] = None,
    ) -> list[EvidenceCluster]:
        """Score each cluster with local NB p-values per side.

        Sets:
          - cluster.local_nb_pvalue_a: p-value at breakpoint A
          - cluster.local_nb_pvalue_b: p-value at breakpoint B
          - cluster.background_p: min(pvalue_a, pvalue_b) for compatibility
          - cluster.local_rate_a: local discordant rate at side A
          - cluster.local_rate_b: local discordant rate at side B

        Does NOT compute a combined p-value — these are per-side.
        """
        if callback:
            callback({
                "type": "scan.stage_changed",
                "stage": "background_model",
            })

        n = len(clusters)
        for idx, cluster in enumerate(clusters):
            # Count evidence touching each side
            observed = cluster.total_support

            pval_a = self.nb_pvalue_at(cluster.chrom_a, cluster.pos_a, observed)
            pval_b = self.nb_pvalue_at(cluster.chrom_b, cluster.pos_b, observed)

            rate_a = self.local_rate_at(cluster.chrom_a, cluster.pos_a)
            rate_b = self.local_rate_at(cluster.chrom_b, cluster.pos_b)

            # Store per-side values
            cluster.local_nb_pvalue_a = pval_a
            cluster.local_nb_pvalue_b = pval_b
            cluster.local_rate_a = rate_a
            cluster.local_rate_b = rate_b

            # Compatibility: use min for downstream scoring
            cluster.background_p = min(pval_a, pval_b)

            if callback and n > 0 and (
                idx % max(1, n // 20) == 0 or idx == n - 1
            ):
                callback({
                    "type": "scan.progress",
                    "stage": "background_model",
                    "pct": round((idx + 1) / n * 100, 1),
                    "detail": (
                        f"Scored {idx + 1}/{n} clusters "
                        f"(local_rate_a={rate_a:.1f}, local_rate_b={rate_b:.1f})"
                    ),
                })

        sig_a = sum(1 for c in clusters if getattr(c, 'local_nb_pvalue_a', 1.0) < 0.05)
        sig_b = sum(1 for c in clusters if getattr(c, 'local_nb_pvalue_b', 1.0) < 0.05)
        logger.info(
            "Local NB scoring: %d/%d significant side A, %d/%d significant side B",
            sig_a, n, sig_b, n,
        )

        return clusters
