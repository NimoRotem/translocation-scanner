"""Breakpoint clustering engine for the translocation scanner.

Groups SV-relevant reads (discordant pairs, split reads, soft-clip pileups)
into evidence clusters that represent candidate translocation breakpoints.

Algorithm overview:
  1. Pre-aggregate discordant reads into 1MB bins; discard singleton bins.
  2. Canonicalize chromosome pairs (chrom_a <= chrom_b lexicographically).
  3. Bucket reads by (canonical pair, orientation).
  4. Within each bucket, sort by genomic position and merge reads whose
     positions fall within a merge distance (500 bp for discordant, 2 bp
     for split reads).
  5. Compute per-cluster statistics: median position, confidence interval,
     support counts, median MAPQ.
  6. Search for reciprocal evidence (A->B paired with B->A).
  7. Assign sequential cluster IDs and return populated EvidenceCluster objects.
"""

from __future__ import annotations

import logging
import resource
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from models import ClipPileup, CompactEvidence, EvidenceCluster, EvidenceType, SVRead

logger = logging.getLogger(__name__)

# Map EvidenceType enum to compact int for CompactEvidence
_EVTYPE_TO_INT = {
    EvidenceType.DISCORDANT: 0,
    EvidenceType.SPLIT: 1,
    EvidenceType.CLIPPED: 2,
}


@dataclass
class _RawCluster:
    """Internal intermediate cluster before promotion to EvidenceCluster.

    Stores per-read data as parallel arrays (same index = same read).
    This avoids holding full SVRead objects in memory while preserving
    all fields needed for downstream filters and audit.
    """

    chrom_a: str
    chrom_b: str
    orientation: str
    # Core parallel arrays (one entry per read)
    positions_a: list[int] = field(default_factory=list)
    positions_b: list[int] = field(default_factory=list)
    mapqs: list[int] = field(default_factory=list)
    # Extended parallel arrays for compact evidence
    strands_a: list[bool] = field(default_factory=list)
    strands_b: list[bool] = field(default_factory=list)
    evidence_types: list[int] = field(default_factory=list)  # 0/1/2
    flags: list[int] = field(default_factory=list)
    read_hashes: list[int] = field(default_factory=list)
    # Cluster-level metadata
    evidence_type: EvidenceType = EvidenceType.DISCORDANT
    clip_pileups: list[ClipPileup] = field(default_factory=list)
    # Evidence counters (exact, derived from evidence_types array)
    discordant_n: int = 0
    split_n: int = 0
    clipped_n: int = 0


class ClusterEngine:
    """Groups SV reads into evidence clusters supporting translocation calls.

    Parameters
    ----------
    merge_distance : int
        Maximum distance (bp) between positions for two discordant reads
        to be merged into the same cluster. Default 500.
    refined_distance : int
        Maximum distance (bp) for merging split-read evidence, which has
        base-pair resolution. Default 2.
    pre_agg_bin_size : int
        Bin size for pre-clustering aggregation (bp). Default 1,000,000.
    pre_agg_min_count : int
        Minimum reads per bin to keep during pre-aggregation. Default 2.
    """

    def __init__(
        self,
        merge_distance: int = 500,
        refined_distance: int = 2,
        pre_agg_bin_size: int = 1_000_000,
        pre_agg_min_count: int = 2,
        min_cluster_support: int = 3,
    ) -> None:
        self.merge_distance = merge_distance
        self.refined_distance = refined_distance
        self._pre_agg_bin_size = pre_agg_bin_size
        self._pre_agg_min_count = pre_agg_min_count
        self._min_cluster_support = min_cluster_support

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def cluster(
        self,
        discordant_reads: list[SVRead],
        split_reads: list[SVRead],
        clip_pileups: list[ClipPileup],
        callback: Optional[Callable[[dict], None]] = None,
        cancel_event: Optional[object] = None,
    ) -> list[EvidenceCluster]:
        """Run the full clustering pipeline and return evidence clusters."""
        t_total = time.monotonic()

        self._emit(callback, {
            "type": "scan.stage_changed",
            "stage": "clustering",
        })
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 0,
            "detail": "Grouping reads by chromosome pair and orientation",
        })

        logger.info(
            "Clustering input: %d discordant, %d split, %d clip pileups",
            len(discordant_reads), len(split_reads), len(clip_pileups),
        )

        def _check_cancel():
            if cancel_event is not None and hasattr(cancel_event, 'is_set') and cancel_event.is_set():
                raise RuntimeError("Clustering cancelled")

        # Step 1: Build raw clusters per evidence type
        t0 = time.monotonic()
        discordant_clusters = self._group_and_merge(
            discordant_reads,
            EvidenceType.DISCORDANT,
            self.merge_distance,
            cancel_event=cancel_event,
        )
        logger.info(
            "Discordant clustering: %.2fs -> %d clusters",
            time.monotonic() - t0, len(discordant_clusters),
        )
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 25,
            "detail": f"Merged {len(discordant_reads)} discordant reads "
                      f"into {len(discordant_clusters)} clusters",
        })
        _check_cancel()

        t0 = time.monotonic()
        split_clusters = self._group_and_merge(
            split_reads,
            EvidenceType.SPLIT,
            self.refined_distance,
            cancel_event=cancel_event,
        )
        logger.info(
            "Split clustering: %.2fs -> %d clusters",
            time.monotonic() - t0, len(split_clusters),
        )
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 50,
            "detail": f"Merged {len(split_reads)} split reads "
                      f"into {len(split_clusters)} clusters",
        })
        _check_cancel()

        t0 = time.monotonic()
        clip_clusters = self._build_clip_clusters(clip_pileups)
        logger.info(
            "Clip clustering: %.2fs -> %d clusters",
            time.monotonic() - t0, len(clip_clusters),
        )
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 65,
            "detail": f"Built {len(clip_clusters)} clusters from "
                      f"{len(clip_pileups)} clip pileups",
        })
        _check_cancel()

        # Step 1b: Prune tiny clusters before expensive cross-merge
        t0 = time.monotonic()
        all_raw = discordant_clusters + split_clusters + clip_clusters
        pre_prune = len(all_raw)
        if self._min_cluster_support > 1:
            all_raw = [
                rc for rc in all_raw
                if len(rc.positions_a) >= self._min_cluster_support
            ]
        logger.info(
            "Cluster pruning (min_support=%d): %d -> %d clusters, %.2fs",
            self._min_cluster_support, pre_prune, len(all_raw),
            time.monotonic() - t0,
        )

        # Step 2: Cross-merge clusters from different evidence types
        t0 = time.monotonic()
        merged = self._cross_merge(all_raw)
        logger.info(
            "Cross-merge: %.2fs, %d -> %d clusters",
            time.monotonic() - t0, len(all_raw), len(merged),
        )
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 80,
            "detail": f"Cross-merged into {len(merged)} unified clusters",
        })

        # Step 3: Convert to EvidenceCluster objects
        t0 = time.monotonic()
        evidence_clusters = self._to_evidence_clusters(merged)
        logger.info(
            "Conversion to EvidenceCluster: %.2fs -> %d clusters",
            time.monotonic() - t0, len(evidence_clusters),
        )

        # Step 4: Detect reciprocal support
        self._annotate_reciprocal(evidence_clusters)
        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 95,
            "detail": "Reciprocal support annotated",
        })

        # Step 5: Assign sequential IDs and sort by total support
        evidence_clusters.sort(key=lambda c: c.total_support, reverse=True)
        for idx, cluster in enumerate(evidence_clusters, start=1):
            cluster.cluster_id = f"CLU_{idx:03d}"

        self._emit(callback, {
            "type": "scan.progress",
            "stage": "clustering",
            "pct": 100,
            "detail": f"Clustering complete: {len(evidence_clusters)} clusters",
        })

        # Summary timing
        try:
            peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        except Exception:
            peak_mb = 0
        logger.info(
            "Clustering total: %.2fs, %d final clusters, peak RSS %.0f MB",
            time.monotonic() - t_total, len(evidence_clusters), peak_mb,
        )
        return evidence_clusters

    # ------------------------------------------------------------------
    # Pre-clustering aggregation
    # ------------------------------------------------------------------

    def _pre_aggregate(
        self,
        reads: list[SVRead],
        bin_size: int,
        min_count: int,
    ) -> list[SVRead]:
        """Filter out low-signal bins before fine-clustering.

        Groups reads by (chrA, chrB, orientation, binA, binB) where bin
        is position // bin_size. Only reads belonging to bins with
        count >= min_count are kept.
        """
        if not reads:
            return reads

        t0 = time.monotonic()

        # Phase 1: Count reads per bin
        bin_counts: dict[tuple, int] = {}
        read_bins: list[tuple] = []
        for read in reads:
            if read.mate_chrom is None or read.mate_pos is None:
                read_bins.append(())
                continue
            canon_a, canon_b = _canonical_pair(read.chrom, read.mate_chrom)
            orient = _orientation_from_read(read)
            if canon_a != read.chrom:
                orient = orient[1] + orient[0]
            bin_a = read.pos // bin_size
            bin_b = read.mate_pos // bin_size
            key = (canon_a, canon_b, orient, bin_a, bin_b)
            read_bins.append(key)
            bin_counts[key] = bin_counts.get(key, 0) + 1

        # Phase 2: Build set of passing bins
        passing_bins = {k for k, v in bin_counts.items() if v >= min_count}

        # Phase 3: Filter reads
        filtered = [
            read for read, key in zip(reads, read_bins)
            if key and key in passing_bins
        ]

        logger.info(
            "Pre-aggregate: %d -> %d reads (%.1f%% kept), %d/%d bins passed "
            "threshold %d, %.2fs",
            len(reads), len(filtered),
            100 * len(filtered) / len(reads) if reads else 0,
            len(passing_bins), len(bin_counts), min_count,
            time.monotonic() - t0,
        )
        return filtered

    # ------------------------------------------------------------------
    # Grouping and merging
    # ------------------------------------------------------------------

    def _group_and_merge(
        self,
        reads: list[SVRead],
        evidence_type: EvidenceType,
        distance: int,
        cancel_event: Optional[object] = None,
    ) -> list[_RawCluster]:
        """Group reads by canonical pair + orientation, then merge nearby positions."""

        # Pre-aggregate discordant reads to eliminate noise bins
        if evidence_type == EvidenceType.DISCORDANT and len(reads) > 100_000:
            reads = self._pre_aggregate(
                reads,
                bin_size=self._pre_agg_bin_size,
                min_count=self._pre_agg_min_count,
            )

        t0 = time.monotonic()

        # Bucket: (chrom_a, chrom_b, orientation) -> list[SVRead]
        buckets: dict[tuple[str, str, str], list[SVRead]] = defaultdict(list)

        for i, read in enumerate(reads):
            # Periodic cancel check inside the bucketing loop
            if i % 500_000 == 0 and i > 0 and cancel_event is not None:
                if hasattr(cancel_event, 'is_set') and cancel_event.is_set():
                    raise RuntimeError("Clustering cancelled during bucketing")

            if read.mate_chrom is None or read.mate_pos is None:
                continue

            canon_a, canon_b = _canonical_pair(read.chrom, read.mate_chrom)
            orient = _orientation_from_read(read)

            # If the canonical ordering swapped the chromosomes, we must
            # mirror the orientation string so it stays consistent.
            if canon_a != read.chrom:
                orient = orient[1] + orient[0]

            buckets[(canon_a, canon_b, orient)].append(read)

        logger.info(
            "Bucketing %d %s reads: %.2fs, %d buckets",
            len(reads), evidence_type.value,
            time.monotonic() - t0, len(buckets),
        )

        raw_clusters: list[_RawCluster] = []

        for (chrom_a, chrom_b, orient), bucket_reads in buckets.items():
            merged = self._merge_reads_in_bucket(
                bucket_reads,
                chrom_a,
                chrom_b,
                orient,
                evidence_type,
                distance,
            )
            raw_clusters.extend(merged)

        return raw_clusters

    def _merge_reads_in_bucket(
        self,
        reads: list[SVRead],
        chrom_a: str,
        chrom_b: str,
        orientation: str,
        evidence_type: EvidenceType,
        distance: int,
    ) -> list[_RawCluster]:
        """Merge reads within a single (pair, orientation) bucket by position.

        Uses numpy arrays for position extraction + sorting to avoid
        Python-level comparison overhead on millions of reads. Extracts
        all per-read fields needed for CompactEvidence into parallel arrays.
        """
        if not reads:
            return []

        n = len(reads)
        evtype_int = _EVTYPE_TO_INT.get(evidence_type, 0)

        # Extract all per-read fields into numpy arrays
        pos_a_arr = np.empty(n, dtype=np.int64)
        pos_b_arr = np.empty(n, dtype=np.int64)
        mapq_arr = np.empty(n, dtype=np.int32)
        strand_a_arr = np.empty(n, dtype=np.bool_)
        strand_b_arr = np.empty(n, dtype=np.bool_)
        evtype_arr = np.full(n, evtype_int, dtype=np.uint8)
        flag_arr = np.empty(n, dtype=np.uint16)
        hash_arr = np.empty(n, dtype=np.uint64)

        for i, read in enumerate(reads):
            if read.chrom == chrom_a:
                pos_a_arr[i] = read.pos
                pos_b_arr[i] = read.mate_pos  # type: ignore
                strand_a_arr[i] = read.is_reverse
                strand_b_arr[i] = read.mate_is_reverse
            else:
                pos_a_arr[i] = read.mate_pos  # type: ignore
                pos_b_arr[i] = read.pos
                strand_a_arr[i] = read.mate_is_reverse
                strand_b_arr[i] = read.is_reverse
            mapq_arr[i] = read.mapq
            flag_arr[i] = read.flag & 0xFFFF
            try:
                hash_arr[i] = int(read.read_name_hash, 16) if read.read_name_hash else 0
            except (ValueError, TypeError):
                hash_arr[i] = 0

        # Sort by pos_a using numpy argsort (much faster than Python sort)
        order = np.argsort(pos_a_arr, kind='quicksort')
        pos_a_sorted = pos_a_arr[order]
        pos_b_sorted = pos_b_arr[order]
        mapq_sorted = mapq_arr[order]
        strand_a_sorted = strand_a_arr[order]
        strand_b_sorted = strand_b_arr[order]
        evtype_sorted = evtype_arr[order]
        flag_sorted = flag_arr[order]
        hash_sorted = hash_arr[order]

        # Find cluster boundaries using running max
        running_max = np.maximum.accumulate(pos_a_sorted)
        breaks = np.empty(n, dtype=bool)
        breaks[0] = True
        breaks[1:] = (pos_a_sorted[1:] - running_max[:-1]) > distance

        # Assign cluster labels
        cluster_ids = np.cumsum(breaks) - 1
        num_clusters = int(cluster_ids[-1]) + 1

        # Build clusters using numpy boolean indexing
        clusters: list[_RawCluster] = []
        for cid in range(num_clusters):
            mask = cluster_ids == cid
            count = int(mask.sum())

            # Exact counts from evidence type
            sub_evtypes = evtype_sorted[mask]
            disc_n = int((sub_evtypes == 0).sum())
            split_n = int((sub_evtypes == 1).sum())
            clip_n = int((sub_evtypes == 2).sum())

            rc = _RawCluster(
                chrom_a=chrom_a,
                chrom_b=chrom_b,
                orientation=orientation,
                evidence_type=evidence_type,
                positions_a=pos_a_sorted[mask].tolist(),
                positions_b=pos_b_sorted[mask].tolist(),
                mapqs=mapq_sorted[mask].tolist(),
                strands_a=strand_a_sorted[mask].tolist(),
                strands_b=strand_b_sorted[mask].tolist(),
                evidence_types=evtype_sorted[mask].tolist(),
                flags=flag_sorted[mask].tolist(),
                read_hashes=hash_sorted[mask].tolist(),
                discordant_n=disc_n,
                split_n=split_n,
                clipped_n=clip_n,
            )
            clusters.append(rc)

        # Second pass: split clusters that are too spread on the B side
        refined: list[_RawCluster] = []
        for cluster in clusters:
            refined.extend(
                self._split_on_pos_b(cluster, distance)
            )

        return refined

    def _split_on_pos_b(
        self,
        cluster: _RawCluster,
        distance: int,
    ) -> list[_RawCluster]:
        """Further split a cluster if pos_b values span more than ``distance``."""
        if len(cluster.positions_b) <= 1:
            return [cluster]

        # Sort indices by pos_b
        indices = sorted(range(len(cluster.positions_b)),
                         key=lambda i: cluster.positions_b[i])

        sub_clusters: list[_RawCluster] = []
        current_indices: list[int] = [indices[0]]

        for i in range(1, len(indices)):
            prev_pos = cluster.positions_b[indices[i - 1]]
            curr_pos = cluster.positions_b[indices[i]]

            if (curr_pos - prev_pos) > distance:
                sub_clusters.append(
                    self._subset_raw_cluster(cluster, current_indices)
                )
                current_indices = [indices[i]]
            else:
                current_indices.append(indices[i])

        sub_clusters.append(
            self._subset_raw_cluster(cluster, current_indices)
        )
        return sub_clusters

    @staticmethod
    def _subset_raw_cluster(
        cluster: _RawCluster,
        indices: list[int],
    ) -> _RawCluster:
        """Create a new _RawCluster from a subset of indices with exact counts."""
        # Slice all parallel arrays at the same indices
        sub_evtypes = [cluster.evidence_types[i] for i in indices] if cluster.evidence_types else []

        # Exact counts from evidence_types array
        disc_n = sum(1 for e in sub_evtypes if e == 0)
        split_n = sum(1 for e in sub_evtypes if e == 1)
        clip_n = sum(1 for e in sub_evtypes if e == 2)

        return _RawCluster(
            chrom_a=cluster.chrom_a,
            chrom_b=cluster.chrom_b,
            orientation=cluster.orientation,
            positions_a=[cluster.positions_a[i] for i in indices],
            positions_b=[cluster.positions_b[i] for i in indices],
            mapqs=[cluster.mapqs[i] for i in indices],
            strands_a=[cluster.strands_a[i] for i in indices] if cluster.strands_a else [],
            strands_b=[cluster.strands_b[i] for i in indices] if cluster.strands_b else [],
            evidence_types=sub_evtypes,
            flags=[cluster.flags[i] for i in indices] if cluster.flags else [],
            read_hashes=[cluster.read_hashes[i] for i in indices] if cluster.read_hashes else [],
            evidence_type=cluster.evidence_type,
            clip_pileups=[
                cluster.clip_pileups[i]
                for i in indices
                if i < len(cluster.clip_pileups)
            ],
            discordant_n=disc_n,
            split_n=split_n,
            clipped_n=clip_n,
        )

    # ------------------------------------------------------------------
    # Clip pileup handling
    # ------------------------------------------------------------------

    def _build_clip_clusters(
        self,
        pileups: list[ClipPileup],
    ) -> list[_RawCluster]:
        """Convert clip pileups with resolved partner loci into raw clusters."""
        buckets: dict[tuple[str, str], list[ClipPileup]] = defaultdict(list)

        for pileup in pileups:
            if pileup.partner_chrom is None or pileup.partner_pos is None:
                continue
            canon_a, canon_b = _canonical_pair(pileup.chrom, pileup.partner_chrom)
            buckets[(canon_a, canon_b)].append(pileup)

        clusters: list[_RawCluster] = []

        for (chrom_a, chrom_b), group in buckets.items():
            group_sorted = sorted(group, key=lambda p: p.pos)

            current: Optional[_RawCluster] = None
            current_max_pos: int = 0

            for pileup in group_sorted:
                # Determine positions relative to canonical order
                if pileup.chrom == chrom_a:
                    pos_a, pos_b = pileup.pos, pileup.partner_pos  # type: ignore[assignment]
                else:
                    pos_a, pos_b = pileup.partner_pos, pileup.pos  # type: ignore[assignment]

                if current is None or (pos_a - current_max_pos) > self.refined_distance:
                    current = _RawCluster(
                        chrom_a=chrom_a,
                        chrom_b=chrom_b,
                        orientation="++",  # Clip pileups lack strand; default
                        evidence_type=EvidenceType.CLIPPED,
                    )
                    clusters.append(current)
                    current_max_pos = pos_a

                current.positions_a.append(pos_a)
                current.positions_b.append(pos_b)
                current.mapqs.append(0)  # MAPQ not tracked for clips
                current.strands_a.append(False)
                current.strands_b.append(False)
                current.evidence_types.append(2)  # CLIPPED
                current.flags.append(0)
                current.read_hashes.append(0)
                current.clip_pileups.append(pileup)
                current.clipped_n += 1
                current_max_pos = max(current_max_pos, pos_a)

        return clusters

    # ------------------------------------------------------------------
    # Cross-merge across evidence types
    # ------------------------------------------------------------------

    def _cross_merge(
        self,
        raw_clusters: list[_RawCluster],
    ) -> list[_RawCluster]:
        """Merge raw clusters from different evidence types that overlap spatially."""
        if not raw_clusters:
            return []

        # Bucket by (chrom_a, chrom_b, orientation)
        buckets: dict[tuple[str, str, str], list[_RawCluster]] = defaultdict(list)
        for rc in raw_clusters:
            buckets[(rc.chrom_a, rc.chrom_b, rc.orientation)].append(rc)

        merged_all: list[_RawCluster] = []

        for key, group in buckets.items():
            merged_all.extend(self._merge_clusters(group, self.merge_distance))

        return merged_all

    def _merge_clusters(
        self,
        clusters: list[_RawCluster],
        distance: int,
    ) -> list[_RawCluster]:
        """Merge overlapping position windows within a group of raw clusters."""
        if len(clusters) <= 1:
            return list(clusters)

        def _median_pos(positions: list[int]) -> int:
            return int(np.median(positions))

        clusters_sorted = sorted(
            clusters, key=lambda c: _median_pos(c.positions_a)
        )

        merged: list[_RawCluster] = [clusters_sorted[0]]

        for candidate in clusters_sorted[1:]:
            last = merged[-1]

            med_a_last = _median_pos(last.positions_a)
            med_b_last = _median_pos(last.positions_b)
            med_a_cand = _median_pos(candidate.positions_a)
            med_b_cand = _median_pos(candidate.positions_b)

            if (abs(med_a_cand - med_a_last) <= distance
                    and abs(med_b_cand - med_b_last) <= distance):
                # Merge into last — extend all parallel arrays
                last.positions_a.extend(candidate.positions_a)
                last.positions_b.extend(candidate.positions_b)
                last.mapqs.extend(candidate.mapqs)
                last.strands_a.extend(candidate.strands_a)
                last.strands_b.extend(candidate.strands_b)
                last.evidence_types.extend(candidate.evidence_types)
                last.flags.extend(candidate.flags)
                last.read_hashes.extend(candidate.read_hashes)
                last.clip_pileups.extend(candidate.clip_pileups)
                last.discordant_n += candidate.discordant_n
                last.split_n += candidate.split_n
                last.clipped_n += candidate.clipped_n
            else:
                merged.append(candidate)

        return merged

    # ------------------------------------------------------------------
    # Conversion to EvidenceCluster
    # ------------------------------------------------------------------

    def _to_evidence_clusters(
        self,
        raw_clusters: list[_RawCluster],
    ) -> list[EvidenceCluster]:
        """Convert internal raw clusters to public EvidenceCluster objects."""
        results: list[EvidenceCluster] = []

        for rc in raw_clusters:
            if not rc.positions_a or not rc.positions_b:
                continue

            pos_a_arr = np.array(rc.positions_a)
            pos_b_arr = np.array(rc.positions_b)
            mapq_arr = np.array([m for m in rc.mapqs if m > 0]) if rc.mapqs else np.array([0])

            median_a = int(np.median(pos_a_arr))
            median_b = int(np.median(pos_b_arr))

            # Confidence interval: range of observed positions relative to median
            ci_a = (int(pos_a_arr.min() - median_a), int(pos_a_arr.max() - median_a))
            ci_b = (int(pos_b_arr.min() - median_b), int(pos_b_arr.max() - median_b))

            # If CI is zero (single read), use evidence-type defaults
            if ci_a == (0, 0):
                if rc.evidence_type == EvidenceType.DISCORDANT:
                    ci_a = (-500, 500)
                else:
                    ci_a = (-2, 2)
            if ci_b == (0, 0):
                if rc.evidence_type == EvidenceType.DISCORDANT:
                    ci_b = (-500, 500)
                else:
                    ci_b = (-2, 2)

            # Exact counts from counters
            discordant_count = rc.discordant_n
            split_count = rc.split_n
            clipped_count = len(rc.clip_pileups) + rc.clipped_n

            median_mapq = float(np.median(mapq_arr)) if len(mapq_arr) > 0 else 0.0

            # Build compact evidence records from parallel arrays
            n_reads = len(rc.positions_a)
            reads: list = []
            if n_reads > 0 and rc.strands_a:
                reads = [
                    CompactEvidence(
                        pos_a=rc.positions_a[i],
                        pos_b=rc.positions_b[i],
                        mapq=rc.mapqs[i],
                        is_reverse=rc.strands_a[i],
                        mate_is_reverse=rc.strands_b[i],
                        evidence_type=rc.evidence_types[i] if rc.evidence_types else 0,
                        flag=rc.flags[i] if rc.flags else 0,
                        read_name_hash=rc.read_hashes[i] if rc.read_hashes else 0,
                    )
                    for i in range(n_reads)
                ]

            # Unique start positions on each side
            unique_a = len(set(rc.positions_a))
            unique_b = len(set(rc.positions_b))

            ec = EvidenceCluster(
                chrom_a=rc.chrom_a,
                chrom_b=rc.chrom_b,
                pos_a=median_a,
                pos_b=median_b,
                ci_a=ci_a,
                ci_b=ci_b,
                orientation=rc.orientation,
                discordant_count=discordant_count,
                split_count=split_count,
                clipped_count=clipped_count,
                median_mapq=median_mapq,
                reads=reads,
                unique_starts_a=unique_a,
                unique_starts_b=unique_b,
            )
            results.append(ec)

        return results

    # ------------------------------------------------------------------
    # Reciprocal evidence detection
    # ------------------------------------------------------------------

    @staticmethod
    def _annotate_reciprocal(clusters: list[EvidenceCluster]) -> None:
        """Mark clusters that have reciprocal support (A->B and B->A)."""
        _reciprocal_orient = {
            "++": "--",
            "--": "++",
            "+-": "-+",
            "-+": "+-",
        }

        index: dict[tuple[str, str], list[EvidenceCluster]] = defaultdict(list)
        for cluster in clusters:
            index[(cluster.chrom_a, cluster.chrom_b)].append(cluster)

        for cluster in clusters:
            if cluster.reciprocal_support > 0:
                continue

            expected_orient = _reciprocal_orient.get(cluster.orientation)
            if expected_orient is None:
                continue

            candidates = index.get((cluster.chrom_a, cluster.chrom_b), [])

            for other in candidates:
                if other is cluster:
                    continue
                if other.orientation != expected_orient:
                    continue

                compat_distance = 1000
                if (abs(other.pos_a - cluster.pos_a) <= compat_distance
                        and abs(other.pos_b - cluster.pos_b) <= compat_distance):
                    cluster.reciprocal_support = other.total_support
                    other.reciprocal_support = cluster.total_support
                    break

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _emit(callback: Optional[Callable[[dict], None]], event: dict) -> None:
        """Safely invoke the progress callback."""
        if callback is not None:
            try:
                callback(event)
            except Exception:
                logger.warning("Clustering callback raised an exception", exc_info=True)


# ======================================================================
# Module-level helper functions
# ======================================================================


def _canonical_pair(chrom_a: str, chrom_b: str) -> tuple[str, str]:
    """Return a canonicalized chromosome pair (lexicographically smaller first)."""
    if chrom_a <= chrom_b:
        return (chrom_a, chrom_b)
    return (chrom_b, chrom_a)


def _orientation_from_read(read: SVRead) -> str:
    """Derive the orientation string for a discordant/split read."""
    strand_self = "-" if read.is_reverse else "+"
    strand_mate = "-" if read.mate_is_reverse else "+"
    return strand_self + strand_mate
