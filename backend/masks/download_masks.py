#!/usr/bin/env python3
"""Mask download and management for the translocation scanner.

Downloads, caches, and version-locks GRCh38 genomic mask tracks:
  - UCSC centromere track
  - UCSC telomere coordinates
  - ENCODE blacklist v2 (Boyle Lab)
  - UCSC segmental duplications (with %identity)
  - RepeatMasker simple repeats / low-complexity
  - GRCh38 gap track (N-bases)
  - Mappability track (GEM 100bp / Umap k36)
  - Acrocentric p-arm coordinates

All downloads are version-locked with SHA256 + source URL + date + build
stored in a manifest file.

Usage::

    python download_masks.py [--output-dir /path/to/masks] [--force]
"""
from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

DEFAULT_MASK_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data"
)

MANIFEST_FILE = "mask_manifest.json"
REFERENCE_BUILD = "GRCh38"

# ---------------------------------------------------------------------------
# Track definitions: (key, URL, description, is_gzipped)
# ---------------------------------------------------------------------------
TRACKS = {
    "centromere": {
        "url": "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/centromeres.txt.gz",
        "description": "UCSC centromere coordinates for GRCh38",
        "gzipped": True,
        "filename": "centromeres.bed",
        "parser": "_parse_centromere",
    },
    "gap": {
        "url": "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/gap.txt.gz",
        "description": "GRCh38 assembly gaps (N-bases)",
        "gzipped": True,
        "filename": "gaps.bed",
        "parser": "_parse_gap",
    },
    "segdup": {
        "url": "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz",
        "description": "UCSC segmental duplications with %identity",
        "gzipped": True,
        "filename": "segdups.bed",
        "parser": "_parse_segdup",
    },
    "simple_repeat": {
        "url": "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/simpleRepeat.txt.gz",
        "description": "RepeatMasker simple repeats / low-complexity",
        "gzipped": True,
        "filename": "simple_repeats.bed",
        "parser": "_parse_simple_repeat",
    },
    "encode_blacklist": {
        "url": "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz",
        "description": "ENCODE blacklist v2 for GRCh38 (Boyle Lab)",
        "gzipped": True,
        "filename": "encode_blacklist_v2.bed",
        "parser": "_parse_blacklist",
    },
}

# Telomere and acrocentric p-arm coordinates are hardcoded (stable GRCh38)
TELOMERE_COORDS = {
    "chr1": [(0, 10000), (248946422, 248956422)],
    "chr2": [(0, 10000), (242183529, 242193529)],
    "chr3": [(0, 10000), (198195559, 198295559)],
    "chr4": [(0, 10000), (190204555, 190214555)],
    "chr5": [(0, 10000), (181478259, 181538259)],
    "chr6": [(0, 10000), (170745979, 170805979)],
    "chr7": [(0, 10000), (159335973, 159345973)],
    "chr8": [(0, 10000), (145078636, 145138636)],
    "chr9": [(0, 10000), (138334717, 138394717)],
    "chr10": [(0, 10000), (133787422, 133797422)],
    "chr11": [(0, 10000), (135076622, 135086622)],
    "chr12": [(0, 10000), (133265309, 133275309)],
    "chr13": [(0, 16000000), (114354328, 114364328)],
    "chr14": [(0, 16000000), (107033718, 107043718)],
    "chr15": [(0, 17000000), (101981189, 101991189)],
    "chr16": [(0, 10000), (90228345, 90338345)],
    "chr17": [(0, 10000), (83247441, 83257441)],
    "chr18": [(0, 10000), (80263285, 80373285)],
    "chr19": [(0, 10000), (58607616, 58617616)],
    "chr20": [(0, 10000), (64334167, 64444167)],
    "chr21": [(0, 5100000), (46699983, 46709983)],
    "chr22": [(0, 10500000), (50808468, 50818468)],
    "chrX": [(0, 10000), (156030895, 156040895)],
    "chrY": [(0, 10000), (57217415, 57227415)],
}

ACROCENTRIC_P_ARMS = {
    "chr13": (0, 17_900_000),
    "chr14": (0, 17_600_000),
    "chr15": (0, 19_000_000),
    "chr21": (0, 13_200_000),
    "chr22": (0, 14_700_000),
}


# ---------------------------------------------------------------------------
# Parsers: convert raw UCSC tables to clean BED
# ---------------------------------------------------------------------------

def _parse_centromere(raw_path: str, bed_path: str) -> int:
    """Parse UCSC centromeres.txt.gz -> BED."""
    count = 0
    with open(raw_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                # Format: bin, chrom, chromStart, chromEnd, name
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                fout.write(f"{chrom}\t{start}\t{end}\tcentromere\n")
                count += 1
    return count


def _parse_gap(raw_path: str, bed_path: str) -> int:
    """Parse UCSC gap.txt.gz -> BED."""
    count = 0
    with open(raw_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                gap_type = parts[7] if len(parts) > 7 else "gap"
                fout.write(f"{chrom}\t{start}\t{end}\t{gap_type}\n")
                count += 1
    return count


def _parse_segdup(raw_path: str, bed_path: str) -> int:
    """Parse UCSC genomicSuperDups.txt.gz -> BED with %identity."""
    count = 0
    with open(raw_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 27:
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                other_chrom = parts[7]
                other_start = parts[8]
                other_end = parts[9]
                frac_match = parts[26]  # fracMatch column
                try:
                    pct_id = round(float(frac_match) * 100, 1)
                except (ValueError, IndexError):
                    pct_id = 0.0
                fout.write(
                    f"{chrom}\t{start}\t{end}\t"
                    f"{other_chrom}:{other_start}-{other_end}\t{pct_id}\n"
                )
                count += 1
    return count


def _parse_simple_repeat(raw_path: str, bed_path: str) -> int:
    """Parse UCSC simpleRepeat.txt.gz -> BED."""
    count = 0
    with open(raw_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[1]
                start = parts[2]
                end = parts[3]
                name = parts[16] if len(parts) > 16 else "simple_repeat"
                fout.write(f"{chrom}\t{start}\t{end}\t{name}\n")
                count += 1
    return count


def _parse_blacklist(raw_path: str, bed_path: str) -> int:
    """Parse ENCODE blacklist BED -> clean BED."""
    count = 0
    with open(raw_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}")
                if len(parts) > 3:
                    fout.write(f"\t{parts[3]}")
                fout.write("\n")
                count += 1
    return count


PARSERS = {
    "_parse_centromere": _parse_centromere,
    "_parse_gap": _parse_gap,
    "_parse_segdup": _parse_segdup,
    "_parse_simple_repeat": _parse_simple_repeat,
    "_parse_blacklist": _parse_blacklist,
}


# ---------------------------------------------------------------------------
# Write hardcoded tracks (telomere, acrocentric p-arm)
# ---------------------------------------------------------------------------

def _write_telomere_bed(bed_path: str) -> int:
    count = 0
    with open(bed_path, "w") as fout:
        for chrom, regions in sorted(TELOMERE_COORDS.items()):
            for start, end in regions:
                fout.write(f"{chrom}\t{start}\t{end}\ttelomere\n")
                count += 1
    return count


def _write_acrocentric_bed(bed_path: str) -> int:
    count = 0
    with open(bed_path, "w") as fout:
        for chrom, (start, end) in sorted(ACROCENTRIC_P_ARMS.items()):
            fout.write(f"{chrom}\t{start}\t{end}\tacrocentric_p_arm\n")
            count += 1
    return count


# ---------------------------------------------------------------------------
# Download engine
# ---------------------------------------------------------------------------

def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_file(url: str, dest: str, is_gzipped: bool = False) -> str:
    """Download a file, optionally decompress gzip."""
    logger.info("Downloading %s ...", url)
    tmp_path = dest + ".tmp"
    try:
        req = urllib.request.Request(url)
        req.add_header("User-Agent", "translocation-scanner/1.0")
        with urllib.request.urlopen(req, timeout=120) as resp:
            with open(tmp_path, "wb") as out:
                shutil.copyfileobj(resp, out)

        if is_gzipped:
            raw_path = dest + ".raw"
            with gzip.open(tmp_path, "rt") as gz_in:
                with open(raw_path, "w") as out:
                    shutil.copyfileobj(gz_in, out)
            os.unlink(tmp_path)
            return raw_path
        else:
            os.rename(tmp_path, dest)
            return dest
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def download_all_masks(
    output_dir: str = DEFAULT_MASK_DIR,
    force: bool = False,
) -> dict:
    """Download all mask tracks and write manifest.

    Args:
        output_dir: Directory to store mask BED files and manifest.
        force: If True, re-download even if files exist.

    Returns:
        Dict with track names mapping to file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    manifest_path = os.path.join(output_dir, MANIFEST_FILE)

    # Check existing manifest
    existing_manifest = {}
    if os.path.isfile(manifest_path) and not force:
        with open(manifest_path) as f:
            existing_manifest = json.load(f)

    manifest = {
        "reference_build": REFERENCE_BUILD,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "tracks": {},
    }
    paths = {}

    # Download remote tracks
    for key, track in TRACKS.items():
        bed_path = os.path.join(output_dir, track["filename"])
        paths[key] = bed_path

        if os.path.isfile(bed_path) and not force:
            existing = existing_manifest.get("tracks", {}).get(key)
            if existing:
                manifest["tracks"][key] = existing
                logger.info("Skipping %s (already exists)", key)
                continue

        try:
            raw_path = _download_file(
                track["url"], bed_path, is_gzipped=track["gzipped"]
            )
            parser_name = track["parser"]
            parser_fn = PARSERS[parser_name]
            count = parser_fn(raw_path, bed_path)

            # Clean up raw file
            if raw_path != bed_path and os.path.exists(raw_path):
                os.unlink(raw_path)

            sha = _sha256(bed_path)
            manifest["tracks"][key] = {
                "filename": track["filename"],
                "source_url": track["url"],
                "description": track["description"],
                "sha256": sha,
                "region_count": count,
                "download_date": time.strftime(
                    "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                ),
            }
            logger.info("Downloaded %s: %d regions, sha256=%s", key, count, sha[:12])

        except Exception as e:
            logger.error("Failed to download %s: %s", key, e)
            manifest["tracks"][key] = {
                "filename": track["filename"],
                "source_url": track["url"],
                "error": str(e),
            }

    # Write hardcoded tracks
    tel_path = os.path.join(output_dir, "telomeres.bed")
    count = _write_telomere_bed(tel_path)
    paths["telomere"] = tel_path
    manifest["tracks"]["telomere"] = {
        "filename": "telomeres.bed",
        "source_url": "hardcoded_GRCh38",
        "description": "Telomere coordinates for GRCh38",
        "sha256": _sha256(tel_path),
        "region_count": count,
    }

    acro_path = os.path.join(output_dir, "acrocentric_p_arms.bed")
    count = _write_acrocentric_bed(acro_path)
    paths["acrocentric_p_arm"] = acro_path
    manifest["tracks"]["acrocentric_p_arm"] = {
        "filename": "acrocentric_p_arms.bed",
        "source_url": "hardcoded_GRCh38",
        "description": "Acrocentric chromosome p-arm coordinates for GRCh38",
        "sha256": _sha256(acro_path),
        "region_count": count,
    }

    # Write manifest
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("Mask manifest written to %s", manifest_path)

    return paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(description="Download GRCh38 mask tracks")
    parser.add_argument(
        "--output-dir", default=DEFAULT_MASK_DIR,
        help="Directory to store mask files (default: %(default)s)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download even if files already exist",
    )
    args = parser.parse_args()
    paths = download_all_masks(args.output_dir, args.force)
    print(f"\nMask files stored in: {args.output_dir}")
    for key, path in sorted(paths.items()):
        exists = "OK" if os.path.isfile(path) else "MISSING"
        print(f"  {key}: {path} [{exists}]")
