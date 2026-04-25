"""Mask management for the translocation scanner.

Provides version-locked genomic mask tracks (centromere, telomere,
ENCODE blacklist, segmental duplications, etc.) that are downloaded
once and read from cache at runtime.

Usage::

    from masks import MaskSet
    masks = MaskSet.load("/path/to/masks/data")
    if masks.centromere.overlaps("chr9", 49000000):
        ...
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from typing import Optional

from regions import RegionIndex

logger = logging.getLogger(__name__)

MANIFEST_FILE = "mask_manifest.json"


@dataclass
class SegDupRegion:
    """A segmental duplication region with percent identity."""
    chrom: str
    start: int
    end: int
    other_locus: str
    pct_identity: float


class SegDupIndex:
    """Segmental duplication index with percent identity tracking."""

    def __init__(self) -> None:
        self._base = RegionIndex()
        self._regions: list[SegDupRegion] = []
        # Per-position max identity lookup (chrom, bin_1kb) -> max_pct
        self._identity_map: dict[tuple[str, int], float] = {}

    def add(self, chrom: str, start: int, end: int,
            other_locus: str = "", pct_identity: float = 0.0) -> None:
        self._base.add(chrom, start, end)
        self._regions.append(SegDupRegion(chrom, start, end, other_locus, pct_identity))
        # Track max identity per 1kb bin
        for pos in range(start // 1000, end // 1000 + 1):
            key = (chrom, pos)
            cur = self._identity_map.get(key, 0.0)
            if pct_identity > cur:
                self._identity_map[key] = pct_identity

    def finalize(self) -> None:
        self._base.finalize()

    def overlaps(self, chrom: str, pos: int, margin: int = 0) -> bool:
        return self._base.overlaps(chrom, pos, margin)

    def max_identity_at(self, chrom: str, pos: int) -> float:
        """Return the maximum segdup %identity at a position."""
        key = (chrom, pos // 1000)
        return self._identity_map.get(key, 0.0)

    def region_count(self) -> int:
        return self._base.region_count()


@dataclass
class MaskSet:
    """Collection of all genomic mask tracks."""
    centromere: RegionIndex = field(default_factory=RegionIndex)
    telomere: RegionIndex = field(default_factory=RegionIndex)
    encode_blacklist: RegionIndex = field(default_factory=RegionIndex)
    segdup: SegDupIndex = field(default_factory=SegDupIndex)
    simple_repeat: RegionIndex = field(default_factory=RegionIndex)
    gap: RegionIndex = field(default_factory=RegionIndex)
    acrocentric_p_arm: RegionIndex = field(default_factory=RegionIndex)
    manifest: dict = field(default_factory=dict)
    manifest_version: str = ""

    @classmethod
    def load(cls, mask_dir: str) -> MaskSet:
        """Load all mask tracks from a directory.

        Args:
            mask_dir: Path to the mask data directory containing BED
                files and mask_manifest.json.

        Returns:
            A populated MaskSet instance.
        """
        ms = cls()

        manifest_path = os.path.join(mask_dir, MANIFEST_FILE)
        if os.path.isfile(manifest_path):
            with open(manifest_path) as f:
                ms.manifest = json.load(f)
            ms.manifest_version = ms.manifest.get("created_at", "unknown")
        else:
            logger.warning("No mask manifest found at %s", manifest_path)

        # Load each track
        _load_bed(ms.centromere, mask_dir, "centromeres.bed", "centromere")
        _load_bed(ms.telomere, mask_dir, "telomeres.bed", "telomere")
        _load_bed(ms.encode_blacklist, mask_dir, "encode_blacklist_v2.bed", "encode_blacklist")
        _load_segdup(ms.segdup, mask_dir, "segdups.bed")
        _load_bed(ms.simple_repeat, mask_dir, "simple_repeats.bed", "simple_repeat")
        _load_bed(ms.gap, mask_dir, "gaps.bed", "gap")
        _load_bed(ms.acrocentric_p_arm, mask_dir, "acrocentric_p_arms.bed", "acrocentric_p_arm")

        total = sum([
            ms.centromere.region_count(),
            ms.telomere.region_count(),
            ms.encode_blacklist.region_count(),
            ms.segdup.region_count(),
            ms.simple_repeat.region_count(),
            ms.gap.region_count(),
            ms.acrocentric_p_arm.region_count(),
        ])
        logger.info(
            "MaskSet loaded: %d total regions from %s (version: %s)",
            total, mask_dir, ms.manifest_version,
        )
        return ms

    def get_overlaps(self, chrom: str, pos: int) -> list[str]:
        """Return list of mask names that overlap a position."""
        overlaps = []
        if self.centromere.overlaps(chrom, pos):
            overlaps.append("centromere")
        if self.telomere.overlaps(chrom, pos):
            overlaps.append("telomere")
        if self.encode_blacklist.overlaps(chrom, pos):
            overlaps.append("encode_blacklist")
        if self.segdup.overlaps(chrom, pos):
            pct = self.segdup.max_identity_at(chrom, pos)
            overlaps.append(f"segdup_{pct:.0f}pct")
        if self.simple_repeat.overlaps(chrom, pos):
            overlaps.append("simple_repeat")
        if self.gap.overlaps(chrom, pos):
            overlaps.append("gap")
        if self.acrocentric_p_arm.overlaps(chrom, pos):
            overlaps.append("acrocentric_p_arm")
        return overlaps


def _load_bed(index: RegionIndex, mask_dir: str, filename: str, name: str) -> None:
    """Load a BED file into a RegionIndex."""
    path = os.path.join(mask_dir, filename)
    if not os.path.isfile(path):
        logger.warning("Mask file not found: %s", path)
        return
    index.load_bed(path)
    logger.info("Loaded %s: %d regions", name, index.region_count())


def _load_segdup(index: SegDupIndex, mask_dir: str, filename: str) -> None:
    """Load segdup BED with %identity into SegDupIndex."""
    path = os.path.join(mask_dir, filename)
    if not os.path.isfile(path):
        logger.warning("Segdup mask file not found: %s", path)
        return
    count = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            other = parts[3] if len(parts) > 3 else ""
            pct = float(parts[4]) if len(parts) > 4 else 0.0
            index.add(chrom, start, end, other, pct)
            count += 1
    index.finalize()
    logger.info("Loaded segdup: %d regions", count)
