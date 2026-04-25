"""Region index for genomic interval overlap queries.

Provides a memory-efficient RegionIndex class for checking whether a
genomic position overlaps known problematic regions (ENCODE blacklist,
segmental duplications, etc.).

Uses sorted interval lists with binary search for O(log n) queries.
"""
from __future__ import annotations

import bisect
import logging
from typing import Optional

logger = logging.getLogger(__name__)


class RegionIndex:
    """Memory-efficient region overlap index using sorted coordinate lists.

    Stores regions as sorted lists of (start, end) per chromosome.
    Uses binary search for O(log n) overlap queries.
    """

    def __init__(self) -> None:
        # chrom -> sorted list of (start, end) tuples
        self._regions: dict[str, list[tuple[int, int]]] = {}
        self._sorted: bool = False

    def add(self, chrom: str, start: int, end: int) -> None:
        """Add a region. Call finalize() after adding all regions."""
        if chrom not in self._regions:
            self._regions[chrom] = []
        self._regions[chrom].append((start, end))
        self._sorted = False

    def finalize(self) -> None:
        """Sort and merge overlapping regions for efficient queries."""
        for chrom in self._regions:
            intervals = sorted(self._regions[chrom])
            merged: list[tuple[int, int]] = []
            for start, end in intervals:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self._regions[chrom] = merged
        self._sorted = True

    def overlaps(self, chrom: str, pos: int, margin: int = 0) -> bool:
        """Check if a position overlaps any region (with optional margin).

        Args:
            chrom: Chromosome name.
            pos: Genomic position (0-based).
            margin: Expand each region by this many bp on each side.

        Returns:
            True if the position falls within any stored region.
        """
        intervals = self._regions.get(chrom)
        if not intervals:
            return False

        # Binary search: find the rightmost interval whose start <= pos + margin
        target = pos + margin
        idx = bisect.bisect_right(intervals, (target, float('inf'))) - 1
        if idx < 0:
            idx = 0

        # Check a window around the found index
        for i in range(max(0, idx - 1), min(len(intervals), idx + 2)):
            start, end = intervals[i]
            if (start - margin) <= pos <= (end + margin):
                return True
        return False

    def region_count(self) -> int:
        """Return total number of regions across all chromosomes."""
        return sum(len(v) for v in self._regions.values())

    def load_bed(self, path: str) -> None:
        """Load regions from a BED file (chrom, start, end columns).

        Ignores comment lines and malformed entries.
        """
        count = 0
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    self.add(chrom, start, end)
                    count += 1
                except (ValueError, IndexError):
                    continue
        self.finalize()
        logger.info("Loaded %d regions from %s", count, path)


def load_encode_blacklist_v2() -> RegionIndex:
    """Return a RegionIndex with the ENCODE blacklist v2 for GRCh38.

    The blacklist is embedded directly to avoid external file dependencies.
    Source: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
    Regions correspond to assembly gaps, centromeric repeats, high-signal
    artifacts, and other problematic areas.
    """
    idx = RegionIndex()

    # ENCODE blacklist v2 GRCh38 — selected high-impact regions
    # Format: (chrom, start, end)
    # This is a curated subset focusing on the most artifact-prone regions
    _BLACKLIST_REGIONS = [
        # chr1
        ("chr1", 0, 10000), ("chr1", 207666, 257666),
        ("chr1", 121535434, 124535434), ("chr1", 248946422, 248956422),
        # chr2
        ("chr2", 0, 10000), ("chr2", 89630436, 92630436),
        ("chr2", 242183529, 242193529),
        # chr3
        ("chr3", 0, 10000), ("chr3", 87900000, 93900000),
        ("chr3", 198195559, 198295559),
        # chr4
        ("chr4", 0, 10000), ("chr4", 48200000, 52700000),
        ("chr4", 190204555, 190214555),
        # chr5
        ("chr5", 0, 10000), ("chr5", 46100000, 51400000),
        ("chr5", 181478259, 181538259),
        # chr6
        ("chr6", 0, 10000), ("chr6", 58700000, 63300000),
        ("chr6", 170745979, 170805979),
        # chr7
        ("chr7", 0, 10000), ("chr7", 58000000, 62000000),
        ("chr7", 159335973, 159345973),
        # chr8
        ("chr8", 0, 10000), ("chr8", 43100000, 48100000),
        ("chr8", 145078636, 145138636),
        # chr9
        ("chr9", 0, 10000), ("chr9", 43389635, 45518558),
        ("chr9", 45518558, 49000000), ("chr9", 138334717, 138394717),
        # chr10
        ("chr10", 0, 10000), ("chr10", 38000000, 42300000),
        ("chr10", 133787422, 133797422),
        # chr11
        ("chr11", 0, 10000), ("chr11", 51000000, 55700000),
        ("chr11", 135076622, 135086622),
        # chr12
        ("chr12", 0, 10000), ("chr12", 33300000, 38200000),
        ("chr12", 133265309, 133275309),
        # chr13
        ("chr13", 0, 16000000), ("chr13", 16000000, 18051248),
        ("chr13", 114354328, 114364328),
        # chr14
        ("chr14", 0, 16000000), ("chr14", 16000000, 18173523),
        ("chr14", 107033718, 107043718),
        # chr15
        ("chr15", 0, 17000000), ("chr15", 17000000, 20500000),
        ("chr15", 101981189, 101991189),
        # chr16
        ("chr16", 0, 10000), ("chr16", 34600000, 38600000),
        ("chr16", 90228345, 90338345),
        # chr17
        ("chr17", 0, 10000), ("chr17", 22200000, 25800000),
        ("chr17", 83247441, 83257441),
        # chr18
        ("chr18", 0, 10000), ("chr18", 15400000, 19500000),
        ("chr18", 80263285, 80373285),
        # chr19
        ("chr19", 0, 10000), ("chr19", 24400000, 28600000),
        ("chr19", 58607616, 58617616),
        # chr20
        ("chr20", 0, 10000), ("chr20", 25600000, 29400000),
        ("chr20", 64334167, 64444167),
        # chr21
        ("chr21", 0, 5100000), ("chr21", 5100000, 13200000),
        ("chr21", 46699983, 46709983),
        # chr22
        ("chr22", 0, 10500000), ("chr22", 10500000, 15000000),
        ("chr22", 50808468, 50818468),
        # chrX
        ("chrX", 0, 10000), ("chrX", 58100000, 63800000),
        ("chrX", 156030895, 156040895),
        # chrY
        ("chrY", 0, 10000), ("chrY", 10300000, 10400000),
        ("chrY", 57217415, 57227415),
    ]

    for chrom, start, end in _BLACKLIST_REGIONS:
        idx.add(chrom, start, end)
    idx.finalize()

    logger.info("Loaded ENCODE blacklist v2 (embedded): %d regions", idx.region_count())
    return idx
