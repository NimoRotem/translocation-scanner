"""SV-read pre-extractor module.

Performs a single streaming pass over a BAM/CRAM file, extracting
structurally-variant-relevant reads: discordant pairs, split reads,
and soft-clipped reads.  Emits real-time telemetry via a callback
for the SSE live-evidence stream.

Usage::

    extractor = SVExtractor("sample.bam", "reference.fasta", callback=emit_sse)
    result = extractor.extract()
    # result keys: discordant_reads, split_reads, clip_pileups,
    #              chrom_progress, library_stats
"""

from __future__ import annotations

import hashlib
import logging
import multiprocessing as mp
import queue as _queue_mod
import threading
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]

try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None  # type: ignore[assignment]

from models import SVRead, EvidenceType, ChromProgress, ClipPileup

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MIN_MAPQ = 20
MIN_CLIP_LEN = 20
MIN_SPLIT_ALIGNED = 30
CLIP_PILEUP_RADIUS = 5
CLIP_PILEUP_MIN_DEPTH = 4
TELEMETRY_INTERVAL_SEC = 0.250
INSERT_SIZE_SAMPLE_CAP = 10_000_000
BIN_SIZE = 100_000  # 100 kb bins for bin_update tracking
BAD_FLAGS = 0xF00  # secondary (0x100) | supplementary (0x800) | duplicate (0x400)
READ_NAME_HASH_LEN = 12

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _normalize_chrom(name: str, *, use_chr_prefix: bool = True) -> str:
    """Normalize a chromosome name to a consistent format.

    Strips or adds the ``chr`` prefix depending on *use_chr_prefix*.
    Also normalizes common aliases (``MT`` <-> ``chrM``).

    Args:
        name: Raw chromosome/contig name from the BAM header or read.
        use_chr_prefix: When ``True`` (default), output uses ``chr1`` style.
            When ``False``, output uses bare ``1`` style.

    Returns:
        The normalized chromosome name string.
    """
    stripped = name
    had_chr = False
    if stripped.startswith("chr"):
        stripped = stripped[3:]
        had_chr = True

    # MT <-> M alias
    if stripped == "MT":
        stripped = "M"
    elif stripped == "M" and not had_chr:
        # bare "M" stays "M"; will become "chrM" if prefix requested
        pass

    if use_chr_prefix:
        return f"chr{stripped}"
    return stripped


def _parse_sa_tag(sa_string: str) -> list[dict[str, Any]]:
    """Parse a ``SA:Z:`` supplementary alignment tag.

    The SA tag is a semicolon-delimited list of supplementary alignments
    in the format ``rname,pos,strand,CIGAR,mapQ,NM;``.

    Args:
        sa_string: The raw SA tag value (without the ``SA:Z:`` prefix).

    Returns:
        A list of dicts with keys ``chrom``, ``pos``, ``strand``,
        ``cigar``, ``mapq``, ``nm``.
    """
    alignments: list[dict[str, Any]] = []
    for entry in sa_string.rstrip(";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(",")
        if len(parts) < 6:
            continue
        try:
            alignments.append(
                {
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "strand": parts[2],
                    "cigar": parts[3],
                    "mapq": int(parts[4]),
                    "nm": int(parts[5]),
                }
            )
        except (ValueError, IndexError):
            continue
    return alignments


def _cigar_aligned_length(cigar_str: str) -> int:
    """Return the total aligned (consumed-reference) length from a CIGAR string.

    Counts M, D, N, =, X operations (opcodes that consume reference bases).
    """
    import re

    total = 0
    for length_str, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in ("M", "D", "N", "=", "X"):
            total += int(length_str)
    return total


def _get_clip_info(
    read: Any,
) -> tuple[Optional[str], int, Optional[str]]:
    """Extract soft-clip information from a pysam read's CIGAR.

    Inspects the first and last CIGAR operations for soft clips (op 4).
    Returns information about the *longer* clip if both ends are clipped.

    Args:
        read: A ``pysam.AlignedSegment`` object.

    Returns:
        A tuple of ``(clip_side, clip_len, clip_seq)`` where *clip_side*
        is ``"left"`` or ``"right"``, *clip_len* is the number of
        soft-clipped bases, and *clip_seq* is the clipped sequence
        (or ``None`` if unavailable).  Returns ``(None, 0, None)`` if
        no significant soft clip is present.
    """
    cigar = read.cigartuples
    if not cigar:
        return None, 0, None

    left_clip = 0
    right_clip = 0

    # CIGAR op 4 = soft clip
    if cigar[0][0] == 4:
        left_clip = cigar[0][1]
    if cigar[-1][0] == 4:
        right_clip = cigar[-1][1]

    if left_clip == 0 and right_clip == 0:
        return None, 0, None

    seq = read.query_sequence or ""

    if left_clip >= right_clip:
        clip_seq = seq[:left_clip] if seq else None
        return "left", left_clip, clip_seq
    else:
        clip_seq = seq[-right_clip:] if seq else None
        return "right", right_clip, clip_seq


def _hash_read_name(name: str) -> str:
    """Return a privacy-preserving hash of a read name.

    Uses the first 12 hex characters of an MD5 digest.

    Args:
        name: The original QNAME string.

    Returns:
        A 12-character hexadecimal string.
    """
    return hashlib.md5(name.encode(), usedforsecurity=False).hexdigest()[
        :READ_NAME_HASH_LEN
    ]


def _canonical_chrom_pair(
    chrom_a: str, chrom_b: str
) -> tuple[str, str]:
    """Return the chromosome pair in canonical (sorted) order."""
    if chrom_a <= chrom_b:
        return chrom_a, chrom_b
    return chrom_b, chrom_a


# ---------------------------------------------------------------------------
# Parallel extraction worker
# ---------------------------------------------------------------------------

_worker_queue: Any = None
_worker_bam_path: Optional[str] = None
_worker_ref_path: Optional[str] = None
_worker_cancel: Any = None


def _init_chrom_worker(
    queue: Any, bam_path: str, ref_path: str, cancel_event: Any = None
) -> None:
    """Pool initializer: store shared state in each worker process."""
    global _worker_queue, _worker_bam_path, _worker_ref_path, _worker_cancel
    _worker_queue = queue
    _worker_bam_path = bam_path
    _worker_ref_path = ref_path
    _worker_cancel = cancel_event


def _extract_single_chrom(args: dict[str, Any]) -> dict[str, Any]:
    """Worker: extract SV reads from one chromosome.

    Opens its own pysam file handle, fetches a single chromosome via
    ``fetch(contig=raw_chrom)``, classifies reads, and returns results.
    Sends periodic ``chrom.progress`` events to the shared queue.
    """
    chrom: str = args["chrom"]
    raw_chrom: str = args["raw_chrom"]
    length: int = args["length"]
    use_chr_prefix: bool = args["use_chr_prefix"]
    valid_chroms: set[str] = args["valid_chroms"]
    min_mapq: int = args["min_mapq"]
    min_clip_length: int = args["min_clip_length"]
    min_split_aligned: int = args["min_split_aligned"]
    insert_size_cap: int = args.get("insert_size_cap", INSERT_SIZE_SAMPLE_CAP)

    debug_regions: list[tuple[str, int, int]] = args.get("debug_regions", [])
    debug_evidence: list[dict[str, Any]] = []

    discordant_reads: list[SVRead] = []
    split_reads: list[SVRead] = []
    clip_reads: list[SVRead] = []
    insert_sizes: list[int] = []
    bin_counts: dict[tuple[str, int], dict[str, int]] = {}
    pair_density: dict[tuple[str, str], int] = {}
    pair_density_last: dict[tuple[str, str], int] = {}
    dirty_bins: set[tuple[str, int]] = set()
    reads_by_pair: dict[tuple[str, str], list[SVRead]] = {}
    _MILESTONES = {1, 5, 10, 25, 50, 100}
    # Read-name deduplication: prevent same read contributing multiple SVReads
    seen_disc: set[str] = set()
    seen_split: set[str] = set()

    def _in_debug_region(c: str, p: int) -> bool:
        for dc, ds, de in debug_regions:
            if c == dc and ds <= p <= de:
                return True
        return False

    def _log_debug(
        read: Any, c: str, mc: str, mp: int, etype: str,
        accepted: bool, reason: str = "", sa_mq: int = 0,
    ) -> None:
        if not debug_regions:
            return
        if not (_in_debug_region(c, read.reference_start) or _in_debug_region(mc, mp)):
            return
        cigar_str = read.cigarstring or ""
        orient = ("-" if read.is_reverse else "+") + (
            "-" if getattr(read, 'mate_is_reverse', False) else "+"
        )
        debug_evidence.append({
            "read_hash": _hash_read_name(read.query_name),
            "chrom_a": c, "pos_a": read.reference_start,
            "chrom_b": mc, "pos_b": mp,
            "mapq": read.mapping_quality, "sa_mapq": sa_mq,
            "evidence_type": etype, "cigar": cigar_str[:80],
            "orientation": orient,
            "accepted": accepted, "rejection_reason": reason,
        })

    reads_processed = 0
    bytes_processed = 0
    last_progress_time = time.monotonic()

    aln_file = pysam.AlignmentFile(
        _worker_bam_path, reference_filename=_worker_ref_path
    )

    try:
        for read in aln_file.fetch(contig=raw_chrom):
            if read.is_unmapped:
                continue

            reads_processed += 1
            bytes_processed += read.query_length or 100

            # ---- Insert sizes ----
            if (
                len(insert_sizes) < insert_size_cap
                and read.is_proper_pair
                and not read.mate_is_unmapped
                and read.reference_name == read.next_reference_name
                and (read.flag & BAD_FLAGS) == 0
                and read.template_length > 0
            ):
                insert_sizes.append(abs(read.template_length))

            flag = read.flag
            is_bad = (flag & BAD_FLAGS) != 0

            # ---- Discordant ----
            if (
                not is_bad
                and read.mapping_quality >= min_mapq
                and not read.mate_is_unmapped
            ):
                mate_raw = read.next_reference_name
                if mate_raw is not None:
                    mate_chrom = _normalize_chrom(
                        mate_raw, use_chr_prefix=use_chr_prefix
                    )
                    if mate_chrom != chrom and mate_chrom in valid_chroms:
                        rn_hash = _hash_read_name(read.query_name)
                        dedup_key = f"{rn_hash}_{chrom}_{mate_chrom}"
                        if dedup_key in seen_disc:
                            _log_debug(
                                read, chrom, mate_chrom, read.next_reference_start,
                                "discordant", accepted=False, reason="dedup",
                            )
                            continue
                        seen_disc.add(dedup_key)
                        sv_read = SVRead(
                            read_name_hash=rn_hash,
                            chrom=chrom,
                            pos=read.reference_start,
                            mapq=read.mapping_quality,
                            evidence_type=EvidenceType.DISCORDANT,
                            mate_chrom=mate_chrom,
                            mate_pos=read.next_reference_start,
                            is_reverse=read.is_reverse,
                            mate_is_reverse=read.mate_is_reverse,
                            insert_size=abs(read.template_length),
                            flag=read.flag,
                        )
                        discordant_reads.append(sv_read)
                        _log_debug(
                            read, chrom, mate_chrom, read.next_reference_start,
                            "discordant", accepted=True,
                        )

                        bin_key = (chrom, read.reference_start // BIN_SIZE)
                        if bin_key not in bin_counts:
                            bin_counts[bin_key] = {
                                "discordant": 0, "split": 0, "clip": 0
                            }
                        bin_counts[bin_key]["discordant"] += 1
                        dirty_bins.add(bin_key)

                        pair_key = _canonical_chrom_pair(chrom, mate_chrom)
                        pair_density[pair_key] = (
                            pair_density.get(pair_key, 0) + 1
                        )
                        reads_by_pair.setdefault(pair_key, []).append(sv_read)

                        if (
                            pair_density[pair_key] in _MILESTONES
                            and _worker_queue is not None
                        ):
                            try:
                                _worker_queue.put_nowait({
                                    "type": "evidence.highlight",
                                    "evidence_type": EvidenceType.DISCORDANT.value,
                                    "chrom_a": pair_key[0],
                                    "pos_a": read.reference_start,
                                    "chrom_b": pair_key[1],
                                    "pos_b": read.next_reference_start,
                                })
                            except Exception:
                                pass

            # ---- Split ----
            if (
                read.mapping_quality >= min_mapq
                and (read.flag & 0xD00) == 0  # primary only: skip secondary|duplicate|supplementary
                and read.has_tag("SA")
            ):
                sa_raw = read.get_tag("SA")
                sa_entries = _parse_sa_tag(sa_raw)
                read_len = read.query_length or read.infer_read_length() or 0

                for sa in sa_entries:
                    sa_chrom = _normalize_chrom(
                        sa["chrom"], use_chr_prefix=use_chr_prefix
                    )
                    if sa_chrom == chrom:
                        continue
                    if sa_chrom not in valid_chroms:
                        _log_debug(
                            read, chrom, sa_chrom, sa["pos"],
                            "split", accepted=False, sa_mq=sa["mapq"],
                            reason="non_primary_chrom",
                        )
                        continue
                    if sa["mapq"] < 20:
                        _log_debug(
                            read, chrom, sa_chrom, sa["pos"],
                            "split", accepted=False, sa_mq=sa["mapq"],
                            reason=f"sa_mapq_low({sa['mapq']})",
                        )
                        continue

                    sa_aligned = _cigar_aligned_length(sa["cigar"])
                    primary_aligned = read.query_alignment_length or 0
                    min_required = (
                        max(min_split_aligned, int(read_len * 0.30))
                        if read_len > 0
                        else min_split_aligned
                    )
                    if (
                        primary_aligned < min_required
                        or sa_aligned < min_required
                    ):
                        _log_debug(
                            read, chrom, sa_chrom, sa["pos"],
                            "split", accepted=False, sa_mq=sa["mapq"],
                            reason=f"aligned_too_short(pri={primary_aligned},sa={sa_aligned},min={min_required})",
                        )
                        continue

                    rn_hash = _hash_read_name(read.query_name)
                    dedup_key = f"{rn_hash}_{chrom}_{sa_chrom}"
                    if dedup_key in seen_split:
                        _log_debug(
                            read, chrom, sa_chrom, sa["pos"],
                            "split", accepted=False, sa_mq=sa["mapq"],
                            reason="dedup",
                        )
                        continue
                    seen_split.add(dedup_key)
                    sv_read = SVRead(
                        read_name_hash=rn_hash,
                        chrom=chrom,
                        pos=read.reference_start,
                        mapq=read.mapping_quality,
                        evidence_type=EvidenceType.SPLIT,
                        mate_chrom=sa_chrom,
                        mate_pos=sa["pos"],
                        mate_mapq=sa["mapq"],
                        sa_tag=sa_raw,
                        is_reverse=read.is_reverse,
                        mate_is_reverse=(sa["strand"] == "-"),
                        flag=read.flag,
                    )
                    split_reads.append(sv_read)
                    _log_debug(
                        read, chrom, sa_chrom, sa["pos"],
                        "split", accepted=True, sa_mq=sa["mapq"],
                    )

                    bin_key = (chrom, read.reference_start // BIN_SIZE)
                    if bin_key not in bin_counts:
                        bin_counts[bin_key] = {
                            "discordant": 0, "split": 0, "clip": 0
                        }
                    bin_counts[bin_key]["split"] += 1
                    dirty_bins.add(bin_key)

                    pair_key = _canonical_chrom_pair(chrom, sa_chrom)
                    pair_density[pair_key] = (
                        pair_density.get(pair_key, 0) + 1
                    )
                    reads_by_pair.setdefault(pair_key, []).append(sv_read)

                    if _worker_queue is not None:
                        try:
                            _worker_queue.put_nowait({
                                "type": "evidence.highlight",
                                "evidence_type": EvidenceType.SPLIT.value,
                                "chrom_a": chrom,
                                "pos_a": read.reference_start,
                                "chrom_b": sa_chrom,
                                "pos_b": sa["pos"],
                            })
                        except Exception:
                            pass

            # ---- Clipped ----
            clip_side, clip_len, clip_seq = _get_clip_info(read)
            if clip_side is not None and clip_len >= min_clip_length:
                sv_read = SVRead(
                    read_name_hash=_hash_read_name(read.query_name),
                    chrom=chrom,
                    pos=read.reference_start,
                    mapq=read.mapping_quality,
                    evidence_type=EvidenceType.CLIPPED,
                    clip_side=clip_side,
                    clip_len=clip_len,
                    clip_seq=clip_seq,
                    is_reverse=read.is_reverse,
                    flag=read.flag,
                )
                clip_reads.append(sv_read)

                bin_key = (chrom, read.reference_start // BIN_SIZE)
                if bin_key not in bin_counts:
                    bin_counts[bin_key] = {
                        "discordant": 0, "split": 0, "clip": 0
                    }
                bin_counts[bin_key]["clip"] += 1
                dirty_bins.add(bin_key)

            # ---- Periodic telemetry ----
            now = time.monotonic()
            if (
                now - last_progress_time >= 1.0
                and _worker_cancel is not None
                and _worker_cancel.is_set()
            ):
                break  # cancelled — exit read loop early
            if now - last_progress_time >= 1.0 and _worker_queue is not None:
                pct = (
                    min(100.0, (read.reference_start / length) * 100.0)
                    if length > 0
                    else 0.0
                )
                try:
                    _worker_queue.put_nowait({
                        "type": "chrom.progress",
                        "chrom": chrom,
                        "pct": round(pct, 1),
                        "reads": reads_processed,
                        "discordant": len(discordant_reads),
                        "split": len(split_reads),
                    })
                    # Pair density deltas
                    for pk, cnt in pair_density.items():
                        prev = pair_density_last.get(pk, 0)
                        if cnt != prev:
                            _worker_queue.put_nowait({
                                "type": "pair.density",
                                "chrom_a": pk[0],
                                "chrom_b": pk[1],
                                "count": cnt,
                                "delta": cnt - prev,
                            })
                            pair_density_last[pk] = cnt
                    # Dirty bin updates
                    for bk in dirty_bins:
                        bc = bin_counts.get(bk)
                        if bc:
                            _worker_queue.put_nowait({
                                "type": "chrom.bin_update",
                                "chrom": bk[0],
                                "bin_start": bk[1] * BIN_SIZE,
                                "bin_end": bk[1] * BIN_SIZE + BIN_SIZE,
                                "discordant": bc["discordant"],
                                "split": bc["split"],
                                "clip": bc["clip"],
                            })
                    dirty_bins.clear()
                except Exception:
                    pass
                last_progress_time = now
    finally:
        aln_file.close()

    return {
        "chrom": chrom,
        "discordant_reads": discordant_reads,
        "split_reads": split_reads,
        "clip_reads": clip_reads,
        "insert_sizes": insert_sizes,
        "bin_counts": bin_counts,
        "pair_density": pair_density,
        "reads_by_pair": reads_by_pair,
        "reads_processed": reads_processed,
        "bytes_processed": bytes_processed,
        "discordant_count": len(discordant_reads),
        "split_count": len(split_reads),
        "clip_count": len(clip_reads),
        "debug_evidence": debug_evidence,
    }


# ---------------------------------------------------------------------------
# Library statistics helper
# ---------------------------------------------------------------------------


@dataclass
class LibraryStats:
    """Insert-size distribution statistics estimated from the BAM."""

    median: float = 0.0
    mean: float = 0.0
    std: float = 0.0
    mad: float = 0.0
    q25: float = 0.0
    q75: float = 0.0
    n_sampled: int = 0


# ---------------------------------------------------------------------------
# SVExtractor
# ---------------------------------------------------------------------------


class SVExtractor:
    """Single-pass SV-relevant read extractor.

    Streams through a coordinate-sorted BAM/CRAM file, classifying each
    read as discordant, split, or clipped and collecting per-chromosome
    progress.  Emits telemetry events via *callback* every ~250 ms so
    a downstream SSE layer can push live updates to the frontend.

    Args:
        bam_path: Path to the input BAM or CRAM file (must be indexed).
        reference_path: Path to the reference FASTA (required for CRAM,
            recommended for BAM).
        callback: Optional callable that receives telemetry event dicts.
            Each dict has a ``"type"`` key identifying the event kind.

    Example::

        def on_event(evt):
            print(evt["type"], evt)

        ext = SVExtractor("tumor.bam", "ref.fa", callback=on_event)
        result = ext.extract()
    """

    def __init__(
        self,
        bam_path: str,
        reference_path: str,
        callback: Optional[Callable[[dict[str, Any]], None]] = None,
        min_mapq: int = MIN_MAPQ,
        min_clip_length: int = MIN_CLIP_LEN,
        min_split_aligned: int = MIN_SPLIT_ALIGNED,
        min_pileup_depth: int = CLIP_PILEUP_MIN_DEPTH,
        pileup_window: int = CLIP_PILEUP_RADIUS,
        exclude_chrM: bool = True,
    ) -> None:
        self.bam_path = bam_path
        self.reference_path = reference_path
        self.callback = callback
        self.min_mapq = min_mapq
        self.min_clip_length = min_clip_length
        self.min_split_aligned = min_split_aligned
        self.min_pileup_depth = min_pileup_depth
        self.pileup_window = pileup_window
        self.exclude_chrM = exclude_chrM
        self._debug_regions: list[tuple[str, int, int]] = []
        self._debug_evidence: list[dict[str, Any]] = []

        # Outputs
        self._discordant_reads: list[SVRead] = []
        self._split_reads: list[SVRead] = []
        self._clip_reads: list[SVRead] = []  # raw, before pileup grouping
        self._clip_pileups: list[dict[str, Any]] = []
        self._chrom_progress: dict[str, ChromProgress] = {}
        self._library_stats = LibraryStats()

        # Internal tracking
        self._insert_sizes: list[int] = []
        self._use_chr_prefix: bool = True
        self._chrom_lengths: dict[str, int] = {}
        self._chrom_order: list[str] = []

        # Bin-level tracking (100 kb bins)
        # Key: (chrom, bin_index)  Value: {discordant, split, clip}
        self._bin_counts: dict[tuple[str, int], dict[str, int]] = defaultdict(
            lambda: {"discordant": 0, "split": 0, "clip": 0}
        )
        # Dirty bins needing emission
        self._dirty_bins: set[tuple[str, int]] = set()

        # Pair density tracking: canonical (chrom_a, chrom_b) -> count
        self._pair_density: dict[tuple[str, str], int] = defaultdict(int)
        self._pair_density_last: dict[tuple[str, str], int] = defaultdict(int)
        self._last_topbins_time: float = 0.0

        # Reads grouped by canonical chrom pair for downstream clustering
        self._reads_by_pair: dict[
            tuple[str, str], list[SVRead]
        ] = defaultdict(list)

        # Read-name deduplication
        self._seen_disc: set[str] = set()
        self._seen_split: set[str] = set()

        # Telemetry state
        self._last_telemetry_time: float = 0.0
        self._telemetry_reads_last: int = 0
        self._telemetry_bytes_last: int = 0
        self._telemetry_disc_last: int = 0
        self._telemetry_split_last: int = 0
        self._total_reads: int = 0
        self._total_bytes: int = 0

    # ------------------------------------------------------------------
    # Debug region support
    # ------------------------------------------------------------------

    def set_debug_regions(
        self,
        region_a: Optional[str],
        region_b: Optional[str],
        margin: int = 2_000_000,
    ) -> None:
        """Configure debug regions for targeted evidence logging.

        Args:
            region_a: Region string like "chr9:100000000".
            region_b: Region string like "chr22:23000000".
            margin: Window size around each position (default 2MB).
        """
        self._debug_regions = []
        for region in [region_a, region_b]:
            if region:
                try:
                    chrom, pos_str = region.split(":")
                    pos = int(pos_str)
                    chrom = _normalize_chrom(chrom, use_chr_prefix=self._use_chr_prefix)
                    self._debug_regions.append((chrom, pos - margin, pos + margin))
                except (ValueError, IndexError):
                    logger.warning("Invalid debug region: %s", region)

    def _is_in_debug_region(self, chrom: str, pos: int) -> bool:
        """Check if a position falls within any debug region."""
        for d_chrom, d_start, d_end in self._debug_regions:
            if chrom == d_chrom and d_start <= pos <= d_end:
                return True
        return False

    def _log_debug_evidence(
        self,
        read: Any,
        chrom: str,
        mate_chrom: str,
        mate_pos: int,
        evidence_type: str,
        accepted: bool,
        rejection_reason: str = "",
        sa_mapq: int = 0,
    ) -> None:
        """Log evidence from a debug region for later output."""
        if not self._debug_regions:
            return
        in_a = self._is_in_debug_region(chrom, read.reference_start)
        in_b = self._is_in_debug_region(mate_chrom, mate_pos)
        if not (in_a or in_b):
            return
        cigar_str = read.cigarstring or ""
        orient = ("-" if read.is_reverse else "+") + ("-" if getattr(read, 'mate_is_reverse', False) else "+")
        self._debug_evidence.append({
            "read_hash": _hash_read_name(read.query_name),
            "chrom_a": chrom,
            "pos_a": read.reference_start,
            "chrom_b": mate_chrom,
            "pos_b": mate_pos,
            "mapq": read.mapping_quality,
            "sa_mapq": sa_mapq,
            "evidence_type": evidence_type,
            "cigar": cigar_str[:80],
            "orientation": orient,
            "accepted": accepted,
            "rejection_reason": rejection_reason,
        })

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def extract(self) -> dict[str, Any]:
        """Run the single-pass extraction over the BAM/CRAM.

        Returns:
            A dict with the following keys:

            - ``discordant_reads``: list of :class:`SVRead` (discordant pairs)
            - ``split_reads``: list of :class:`SVRead` (split/supplementary)
            - ``clip_pileups``: list of clip-pileup dicts
            - ``chrom_progress``: dict mapping chrom name to :class:`ChromProgress`
            - ``library_stats``: a :class:`LibraryStats` instance
            - ``reads_by_pair``: dict mapping canonical ``(chrom_a, chrom_b)``
              to lists of :class:`SVRead`
            - ``pair_density``: dict mapping canonical chrom pair to count
            - ``bin_counts``: dict mapping ``(chrom, bin_idx)`` to evidence counts
        """
        if pysam is None:
            raise RuntimeError(
                "pysam is required for BAM/CRAM extraction but is not installed"
            )

        logger.info("Starting SV extraction from %s", self.bam_path)
        start_time = time.monotonic()
        self._last_telemetry_time = start_time

        aln_file = pysam.AlignmentFile(
            self.bam_path,
            reference_filename=self.reference_path,
        )

        try:
            self._init_chrom_info(aln_file)
            self._scan(aln_file)
        finally:
            aln_file.close()

        # Post-processing
        self._compute_library_stats()
        self._build_clip_pileups()

        elapsed = time.monotonic() - start_time
        logger.info(
            "Extraction complete in %.1fs — %d reads, "
            "%d discordant, %d split, %d clip pileups",
            elapsed,
            self._total_reads,
            len(self._discordant_reads),
            len(self._split_reads),
            len(self._clip_pileups),
        )

        total_reads = sum(p.reads_processed for p in self._chrom_progress.values())
        total_bytes = sum(p.bytes_processed for p in self._chrom_progress.values())

        return {
            "discordant_reads": self._discordant_reads,
            "split_reads": self._split_reads,
            "clip_pileups": self._clip_pileups,
            "chrom_progress": self._chrom_progress,
            "chrom_lengths": self._chrom_lengths,
            "library_stats": self._library_stats,
            "reads_by_pair": dict(self._reads_by_pair),
            "pair_density": dict(self._pair_density),
            "bin_counts": dict(self._bin_counts),
            "total_reads_processed": total_reads,
            "total_bytes_processed": total_bytes,
            "debug_evidence": self._debug_evidence,
        }

    # ------------------------------------------------------------------
    # Parallel extraction
    # ------------------------------------------------------------------

    def extract_parallel(
        self, num_workers: Optional[int] = None, cancel_event: Any = None
    ) -> dict[str, Any]:
        """Run per-chromosome parallel extraction over the BAM/CRAM.

        Each chromosome is processed by a separate worker process using
        ``multiprocessing.Pool``.  Workers open their own pysam file
        handles and use the BAI index to jump directly to their assigned
        chromosome, so there is no lock contention.

        Args:
            num_workers: Number of worker processes.  ``None`` or ``0``
                means auto-detect (``min(cpu_count, 24)``).

        Returns:
            Same dict structure as :meth:`extract`.
        """
        if pysam is None:
            raise RuntimeError(
                "pysam is required for BAM/CRAM extraction but is not installed"
            )

        logger.info(
            "Starting parallel SV extraction from %s", self.bam_path
        )
        start_time = time.monotonic()

        # 1. Read header, init chrom info, build raw-name map
        aln_file = pysam.AlignmentFile(
            self.bam_path, reference_filename=self.reference_path
        )
        try:
            self._init_chrom_info(aln_file)
            raw_chrom_map: dict[str, str] = {}
            for sq in aln_file.header.get("SQ", []):
                name = sq["SN"]
                normalized = _normalize_chrom(
                    name, use_chr_prefix=self._use_chr_prefix
                )
                if normalized in self._chrom_progress:
                    raw_chrom_map[normalized] = name
        finally:
            aln_file.close()

        # 2. Determine worker count
        if num_workers is None or num_workers <= 0:
            num_workers = min(mp.cpu_count() or 4, 24)
        num_workers = min(num_workers, len(self._chrom_order))

        valid_chroms = set(self._chrom_progress.keys())
        insert_size_cap = max(
            100_000, INSERT_SIZE_SAMPLE_CAP // num_workers
        )

        # 3. Build work items
        work_items: list[dict[str, Any]] = []
        for chrom in self._chrom_order:
            work_items.append(
                {
                    "chrom": chrom,
                    "raw_chrom": raw_chrom_map.get(chrom, chrom),
                    "length": self._chrom_lengths.get(chrom, 0),
                    "use_chr_prefix": self._use_chr_prefix,
                    "valid_chroms": valid_chroms,
                    "min_mapq": self.min_mapq,
                    "min_clip_length": self.min_clip_length,
                    "min_split_aligned": self.min_split_aligned,
                    "insert_size_cap": insert_size_cap,
                    "debug_regions": self._debug_regions,
                }
            )

        # 4. Progress queue + consumer thread
        progress_queue: mp.Queue = mp.Queue()
        stop_consumer = threading.Event()

        def _consume_progress() -> None:
            chrom_stats: dict[str, dict] = {}
            last_tp_time = time.monotonic()
            last_topbins_time = time.monotonic()
            prev_reads = 0
            prev_bytes = 0
            prev_disc = 0
            prev_split = 0
            # Track accumulated pair density from worker events
            accumulated_density: dict[tuple[str, str], int] = {}

            while not stop_consumer.is_set():
                try:
                    event = progress_queue.get(timeout=0.25)
                    etype = event.get("type")
                    if etype == "chrom.progress":
                        chrom_stats[event["chrom"]] = event
                    elif etype == "pair.density":
                        # Track accumulated counts for top_bins
                        pk = (event.get("chrom_a", ""), event.get("chrom_b", ""))
                        accumulated_density[pk] = event.get("count", 0)
                    self._emit_event(event)
                except _queue_mod.Empty:
                    pass
                except Exception:
                    break

                # Emit aggregate scan.throughput every ~500ms
                now = time.monotonic()
                elapsed = now - last_tp_time
                if elapsed >= 0.5 and chrom_stats:
                    tot_r = sum(
                        cs.get("reads", 0) for cs in chrom_stats.values()
                    )
                    tot_d = sum(
                        cs.get("discordant", 0)
                        for cs in chrom_stats.values()
                    )
                    tot_s = sum(
                        cs.get("split", 0) for cs in chrom_stats.values()
                    )
                    tot_b = tot_r * 150  # ~150 bytes/read estimate
                    if tot_r > prev_reads:
                        self._emit_event({
                            "type": "scan.throughput",
                            "reads_per_sec": int(
                                (tot_r - prev_reads) / elapsed
                            ),
                            "bytes_per_sec": int(
                                (tot_b - prev_bytes) / elapsed
                            ),
                            "discordant_per_sec": int(
                                (tot_d - prev_disc) / elapsed
                            ),
                            "split_per_sec": int(
                                (tot_s - prev_split) / elapsed
                            ),
                        })
                    prev_reads = tot_r
                    prev_bytes = tot_b
                    prev_disc = tot_d
                    prev_split = tot_s
                    last_tp_time = now

                # Emit provisional.top_bins every ~2s
                if now - last_topbins_time >= 2.0 and accumulated_density:
                    top = sorted(
                        accumulated_density.items(),
                        key=lambda x: x[1],
                        reverse=True,
                    )[:12]
                    bins = [
                        {
                            "chrom_a": pk[0],
                            "chrom_b": pk[1],
                            "count": cnt,
                            "pos_a": 0,
                            "pos_b": 0,
                        }
                        for pk, cnt in top
                        if cnt >= 2
                    ]
                    if bins:
                        self._emit_event({
                            "type": "provisional.top_bins",
                            "bins": bins,
                        })
                    last_topbins_time = now

        consumer = threading.Thread(
            target=_consume_progress, daemon=True, name="parallel-progress"
        )
        consumer.start()

        # 5. Launch pool
        logger.info(
            "Launching %d workers for %d chromosomes",
            num_workers,
            len(work_items),
        )
        for chrom in self._chrom_order:
            self._chrom_progress[chrom].status = "scanning"

        try:
            with mp.Pool(
                num_workers,
                initializer=_init_chrom_worker,
                initargs=(
                    progress_queue, self.bam_path,
                    self.reference_path, cancel_event,
                ),
            ) as pool:
                for result in pool.imap_unordered(
                    _extract_single_chrom, work_items
                ):
                    if cancel_event is not None and cancel_event.is_set():
                        from pipeline import CancelledException
                        raise CancelledException(
                            "Scan cancelled during extraction"
                        )
                    chrom = result["chrom"]

                    # Merge read lists
                    self._discordant_reads.extend(
                        result["discordant_reads"]
                    )
                    self._split_reads.extend(result["split_reads"])
                    self._clip_reads.extend(result["clip_reads"])
                    self._insert_sizes.extend(result["insert_sizes"])

                    # Merge bin counts
                    for key, counts in result["bin_counts"].items():
                        if isinstance(key, list):
                            key = tuple(key)
                        for ev_type, count in counts.items():
                            self._bin_counts[key][ev_type] += count

                    # Merge pair density
                    for pair_key, count in result["pair_density"].items():
                        if isinstance(pair_key, list):
                            pair_key = tuple(pair_key)
                        self._pair_density[pair_key] += count

                    # Merge reads by pair
                    for pair_key, reads in result[
                        "reads_by_pair"
                    ].items():
                        if isinstance(pair_key, list):
                            pair_key = tuple(pair_key)
                        self._reads_by_pair[pair_key].extend(reads)

                    # Merge debug evidence
                    self._debug_evidence.extend(
                        result.get("debug_evidence", [])
                    )

                    # Update chrom progress
                    prog = self._chrom_progress[chrom]
                    prog.reads_processed = result["reads_processed"]
                    prog.bytes_processed = result["bytes_processed"]
                    prog.discordant_count = result["discordant_count"]
                    prog.split_count = result["split_count"]
                    prog.clip_count = result["clip_count"]
                    prog.status = "complete"
                    prog.pct = 100.0

                    self._total_reads += result["reads_processed"]
                    self._total_bytes += result["bytes_processed"]

                    self._emit_chrom_progress(chrom)
                    logger.info(
                        "Chromosome %s complete: %d reads, "
                        "%d disc, %d split, %d clip",
                        chrom,
                        result["reads_processed"],
                        result["discordant_count"],
                        result["split_count"],
                        result["clip_count"],
                    )
        finally:
            # 6. Stop consumer, drain remaining events
            stop_consumer.set()
            consumer.join(timeout=2.0)
            try:
                while not progress_queue.empty():
                    event = progress_queue.get_nowait()
                    self._emit_event(event)
            except Exception:
                pass

        # 7. Subsample insert sizes if needed
        if len(self._insert_sizes) > INSERT_SIZE_SAMPLE_CAP:
            import random

            self._insert_sizes = random.sample(
                self._insert_sizes, INSERT_SIZE_SAMPLE_CAP
            )

        # 8. Post-processing (same as serial)
        self._compute_library_stats()
        self._build_clip_pileups()

        elapsed = time.monotonic() - start_time
        logger.info(
            "Parallel extraction complete in %.1fs — %d reads, "
            "%d discordant, %d split, %d clip pileups (%d workers)",
            elapsed,
            self._total_reads,
            len(self._discordant_reads),
            len(self._split_reads),
            len(self._clip_pileups),
            num_workers,
        )

        total_reads = sum(
            p.reads_processed for p in self._chrom_progress.values()
        )
        total_bytes = sum(
            p.bytes_processed for p in self._chrom_progress.values()
        )

        return {
            "discordant_reads": self._discordant_reads,
            "split_reads": self._split_reads,
            "clip_pileups": self._clip_pileups,
            "chrom_progress": self._chrom_progress,
            "chrom_lengths": self._chrom_lengths,
            "library_stats": self._library_stats,
            "reads_by_pair": dict(self._reads_by_pair),
            "pair_density": dict(self._pair_density),
            "bin_counts": dict(self._bin_counts),
            "total_reads_processed": total_reads,
            "total_bytes_processed": total_bytes,
            "debug_evidence": self._debug_evidence,
        }

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def _init_chrom_info(self, aln_file: Any) -> None:
        """Detect chr-prefix convention and populate chromosome metadata."""
        header = aln_file.header
        sq_entries = header.get("SQ", [])
        if not sq_entries:
            raise ValueError("BAM/CRAM header contains no @SQ entries")

        # Detect prefix convention from the first standard chromosome
        self._use_chr_prefix = any(
            sq["SN"].startswith("chr") for sq in sq_entries
        )

        # Build ordered chromosome list (primary assembly only)
        valid_chroms = set()
        for i in range(1, 23):
            valid_chroms.add(str(i))
            valid_chroms.add(f"chr{i}")
        valid_chroms.update(["X", "Y", "chrX", "chrY"])
        if not self.exclude_chrM:
            valid_chroms.update(["M", "MT", "chrM"])

        for sq in sq_entries:
            name = sq["SN"]
            if name in valid_chroms:
                normalized = _normalize_chrom(
                    name, use_chr_prefix=self._use_chr_prefix
                )
                length = sq.get("LN", 0)
                self._chrom_lengths[normalized] = length
                self._chrom_order.append(normalized)
                self._chrom_progress[normalized] = ChromProgress(
                    chrom=normalized,
                    length=length,
                    status="pending",
                )

        logger.info(
            "Detected %d primary chromosomes, chr_prefix=%s",
            len(self._chrom_order),
            self._use_chr_prefix,
        )

    # ------------------------------------------------------------------
    # Core scan loop
    # ------------------------------------------------------------------

    def _scan(self, aln_file: Any) -> None:
        """Iterate over every read in the BAM, classifying SV-relevant ones."""
        current_chrom: Optional[str] = None

        for read in aln_file.fetch(until_eof=True):
            # Skip unmapped reads entirely
            if read.is_unmapped:
                continue

            raw_chrom = read.reference_name
            if raw_chrom is None:
                continue

            normalized_chrom = _normalize_chrom(
                raw_chrom, use_chr_prefix=self._use_chr_prefix
            )

            # Skip non-primary chromosomes (alt contigs, decoys, etc.)
            if normalized_chrom not in self._chrom_progress:
                continue

            self._total_reads += 1
            # Estimate bytes from query length as approximation
            self._total_bytes += read.query_length or 100

            # Track chromosome transitions
            if normalized_chrom != current_chrom:
                if current_chrom is not None and current_chrom in self._chrom_progress:
                    self._chrom_progress[current_chrom].status = "complete"
                    self._chrom_progress[current_chrom].pct = 100.0
                    self._emit_chrom_progress(current_chrom)
                current_chrom = normalized_chrom
                self._chrom_progress[current_chrom].status = "scanning"

            # Update chrom read count
            prog = self._chrom_progress[normalized_chrom]
            prog.reads_processed += 1
            if prog.length > 0:
                prog.pct = min(
                    100.0,
                    (read.reference_start / prog.length) * 100.0,
                )

            # ---- Collect insert sizes for library stats ----
            if (
                len(self._insert_sizes) < INSERT_SIZE_SAMPLE_CAP
                and read.is_proper_pair
                and not read.is_unmapped
                and not read.mate_is_unmapped
                and read.reference_name == read.next_reference_name
                and (read.flag & BAD_FLAGS) == 0
                and read.template_length > 0
            ):
                self._insert_sizes.append(abs(read.template_length))

            # ---- Classify reads ----
            flag = read.flag
            is_bad = (flag & BAD_FLAGS) != 0

            # Discordant: mate on different chrom, good MAPQ, no bad flags
            if self._is_discordant(read, normalized_chrom, is_bad):
                self._handle_discordant(read, normalized_chrom)

            # Split: has SA tag
            if self._is_split(read, is_bad):
                self._handle_split(read, normalized_chrom)

            # Clipped: significant soft clip
            clip_side, clip_len, clip_seq = _get_clip_info(read)
            if clip_side is not None and clip_len >= self.min_clip_length:
                self._handle_clipped(
                    read, normalized_chrom, clip_side, clip_len, clip_seq
                )

            # ---- Periodic telemetry ----
            now = time.monotonic()
            if now - self._last_telemetry_time >= TELEMETRY_INTERVAL_SEC:
                self._emit_telemetry(now)
                self._last_telemetry_time = now

        # Finalize last chromosome
        if current_chrom is not None and current_chrom in self._chrom_progress:
            self._chrom_progress[current_chrom].status = "complete"
            self._chrom_progress[current_chrom].pct = 100.0
            self._emit_chrom_progress(current_chrom)

        # Final telemetry flush
        self._emit_telemetry(time.monotonic())
        self._flush_dirty_bins()

    # ------------------------------------------------------------------
    # Classification predicates
    # ------------------------------------------------------------------

    def _is_discordant(
        self, read: Any, chrom: str, is_bad: bool
    ) -> bool:
        """Check whether a read qualifies as a discordant pair."""
        if is_bad:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if read.mate_is_unmapped:
            return False
        mate_chrom_raw = read.next_reference_name
        if mate_chrom_raw is None:
            return False
        mate_chrom = _normalize_chrom(
            mate_chrom_raw, use_chr_prefix=self._use_chr_prefix
        )
        # Must be on a different chromosome
        if mate_chrom == chrom:
            return False
        # Skip mates landing on non-primary chromosomes
        if mate_chrom not in self._chrom_progress:
            return False
        # Mate MAPQ check — we can only get this from the mate's own
        # record during the scan, but pysam does not expose it from the
        # current record.  We apply MAPQ >= 20 on the primary read only;
        # mates are filtered at clustering time.
        return True

    def _is_split(self, read: Any, is_bad: bool) -> bool:
        """Check whether a read qualifies as a split read via SA tag."""
        if read.mapping_quality < self.min_mapq:
            return False
        # Only process primary alignments — exclude secondary, duplicate,
        # and supplementary.  Processing only primaries prevents double-
        # counting where both primary and supplementary create SVRead
        # entries for the same breakpoint.
        if (read.flag & 0xD00) != 0:  # secondary | duplicate | supplementary
            return False
        if not read.has_tag("SA"):
            return False
        return True

    # ------------------------------------------------------------------
    # Read handlers
    # ------------------------------------------------------------------

    def _handle_discordant(self, read: Any, chrom: str) -> None:
        """Process a discordant read and store it."""
        mate_chrom = _normalize_chrom(
            read.next_reference_name, use_chr_prefix=self._use_chr_prefix
        )

        # Deduplicate by read name + chrom pair
        rn_hash = _hash_read_name(read.query_name)
        dedup_key = f"{rn_hash}_{chrom}_{mate_chrom}"
        if dedup_key in self._seen_disc:
            self._log_debug_evidence(
                read, chrom, mate_chrom, read.next_reference_start,
                "discordant", accepted=False, rejection_reason="dedup",
            )
            return
        self._seen_disc.add(dedup_key)

        sv_read = SVRead(
            read_name_hash=rn_hash,
            chrom=chrom,
            pos=read.reference_start,
            mapq=read.mapping_quality,
            evidence_type=EvidenceType.DISCORDANT,
            mate_chrom=mate_chrom,
            mate_pos=read.next_reference_start,
            is_reverse=read.is_reverse,
            mate_is_reverse=read.mate_is_reverse,
            insert_size=abs(read.template_length),
            flag=read.flag,
        )

        self._discordant_reads.append(sv_read)

        # Debug logging for accepted discordant reads
        self._log_debug_evidence(
            read, chrom, mate_chrom, read.next_reference_start,
            "discordant", accepted=True,
        )

        # Update progress counters
        self._chrom_progress[chrom].discordant_count += 1

        # Bin tracking
        bin_idx = read.reference_start // BIN_SIZE
        self._bin_counts[(chrom, bin_idx)]["discordant"] += 1
        self._dirty_bins.add((chrom, bin_idx))

        # Pair density
        pair_key = _canonical_chrom_pair(chrom, mate_chrom)
        self._pair_density[pair_key] += 1

        # Group by canonical pair for clustering
        self._reads_by_pair[pair_key].append(sv_read)

        # Emit highlight for new pair clusters or significant events
        if self._pair_density[pair_key] in (1, 5, 10, 25, 50, 100):
            self._emit_event(
                {
                    "type": "evidence.highlight",
                    "evidence_type": EvidenceType.DISCORDANT.value,
                    "chrom_a": pair_key[0],
                    "pos_a": read.reference_start,
                    "chrom_b": pair_key[1],
                    "pos_b": read.next_reference_start,
                }
            )

    def _handle_split(self, read: Any, chrom: str) -> None:
        """Process a split read (SA tag present) and store it."""
        sa_raw = read.get_tag("SA")
        sa_entries = _parse_sa_tag(sa_raw)
        if not sa_entries:
            return

        read_len = read.query_length or read.infer_read_length() or 0
        rn_hash = _hash_read_name(read.query_name)

        for sa in sa_entries:
            sa_chrom = _normalize_chrom(
                sa["chrom"], use_chr_prefix=self._use_chr_prefix
            )

            # Skip same-chromosome supplementary alignments
            if sa_chrom == chrom:
                continue

            # Skip non-primary chromosome targets
            if sa_chrom not in self._chrom_progress:
                self._log_debug_evidence(
                    read, chrom, sa_chrom, sa["pos"],
                    "split", accepted=False, sa_mapq=sa["mapq"],
                    rejection_reason="non_primary_chrom",
                )
                continue

            # Check supplementary MAPQ (require >= 20 to reduce noise)
            if sa["mapq"] < 20:
                self._log_debug_evidence(
                    read, chrom, sa_chrom, sa["pos"],
                    "split", accepted=False, sa_mapq=sa["mapq"],
                    rejection_reason=f"sa_mapq_low({sa['mapq']})",
                )
                continue

            # Check minimum aligned length on each side
            sa_aligned = _cigar_aligned_length(sa["cigar"])
            primary_aligned = (
                read.query_alignment_length
                if read.query_alignment_length
                else 0
            )

            min_required = max(self.min_split_aligned, int(read_len * 0.30)) if read_len > 0 else self.min_split_aligned

            if primary_aligned < min_required or sa_aligned < min_required:
                self._log_debug_evidence(
                    read, chrom, sa_chrom, sa["pos"],
                    "split", accepted=False, sa_mapq=sa["mapq"],
                    rejection_reason=f"aligned_too_short(pri={primary_aligned},sa={sa_aligned},min={min_required})",
                )
                continue

            # Deduplicate by read name + chrom pair
            dedup_key = f"{rn_hash}_{chrom}_{sa_chrom}"
            if dedup_key in self._seen_split:
                self._log_debug_evidence(
                    read, chrom, sa_chrom, sa["pos"],
                    "split", accepted=False, sa_mapq=sa["mapq"],
                    rejection_reason="dedup",
                )
                continue
            self._seen_split.add(dedup_key)

            sv_read = SVRead(
                read_name_hash=rn_hash,
                chrom=chrom,
                pos=read.reference_start,
                mapq=read.mapping_quality,
                evidence_type=EvidenceType.SPLIT,
                mate_chrom=sa_chrom,
                mate_pos=sa["pos"],
                mate_mapq=sa["mapq"],
                sa_tag=sa_raw,
                is_reverse=read.is_reverse,
                mate_is_reverse=(sa["strand"] == "-"),
                flag=read.flag,
            )

            self._split_reads.append(sv_read)

            # Debug logging for accepted split reads
            self._log_debug_evidence(
                read, chrom, sa_chrom, sa["pos"],
                "split", accepted=True, sa_mapq=sa["mapq"],
            )

            # Update progress
            self._chrom_progress[chrom].split_count += 1

            # Bin tracking
            bin_idx = read.reference_start // BIN_SIZE
            self._bin_counts[(chrom, bin_idx)]["split"] += 1
            self._dirty_bins.add((chrom, bin_idx))

            # Pair density
            pair_key = _canonical_chrom_pair(chrom, sa_chrom)
            self._pair_density[pair_key] += 1
            self._reads_by_pair[pair_key].append(sv_read)

            # Split reads are always emitted as highlights (rare signal)
            self._emit_event(
                {
                    "type": "evidence.highlight",
                    "evidence_type": EvidenceType.SPLIT.value,
                    "chrom_a": chrom,
                    "pos_a": read.reference_start,
                    "chrom_b": sa_chrom,
                    "pos_b": sa["pos"],
                }
            )

    def _handle_clipped(
        self,
        read: Any,
        chrom: str,
        clip_side: str,
        clip_len: int,
        clip_seq: Optional[str],
    ) -> None:
        """Process a soft-clipped read and store it."""
        sv_read = SVRead(
            read_name_hash=_hash_read_name(read.query_name),
            chrom=chrom,
            pos=read.reference_start,
            mapq=read.mapping_quality,
            evidence_type=EvidenceType.CLIPPED,
            clip_side=clip_side,
            clip_len=clip_len,
            clip_seq=clip_seq,
            is_reverse=read.is_reverse,
            flag=read.flag,
        )

        self._clip_reads.append(sv_read)

        # Update progress
        self._chrom_progress[chrom].clip_count += 1

        # Bin tracking
        bin_idx = read.reference_start // BIN_SIZE
        self._bin_counts[(chrom, bin_idx)]["clip"] += 1
        self._dirty_bins.add((chrom, bin_idx))

    # ------------------------------------------------------------------
    # Library statistics
    # ------------------------------------------------------------------

    def _compute_library_stats(self) -> None:
        """Compute insert-size distribution from sampled properly-paired reads."""
        n = len(self._insert_sizes)
        if n == 0:
            logger.warning("No properly-paired reads sampled for insert-size estimation")
            return

        if np is not None:
            arr = np.array(self._insert_sizes, dtype=np.float64)
            median = float(np.median(arr))
            mean = float(np.mean(arr))
            std = float(np.std(arr))
            mad = float(np.median(np.abs(arr - median)))
            q25 = float(np.percentile(arr, 25))
            q75 = float(np.percentile(arr, 75))
        else:
            # Pure-Python fallback (less efficient, development only)
            sorted_sizes = sorted(self._insert_sizes)
            median = float(sorted_sizes[n // 2])
            mean = sum(self._insert_sizes) / n
            variance = sum((x - mean) ** 2 for x in self._insert_sizes) / n
            std = variance ** 0.5
            mad = float(
                sorted(abs(x - median) for x in self._insert_sizes)[n // 2]
            )
            q25 = float(sorted_sizes[n // 4])
            q75 = float(sorted_sizes[3 * n // 4])

        self._library_stats = LibraryStats(
            median=median,
            mean=mean,
            std=std,
            mad=mad,
            q25=q25,
            q75=q75,
            n_sampled=n,
        )

        logger.info(
            "Library insert-size stats: median=%.0f, std=%.0f, n=%d",
            median,
            std,
            n,
        )

        # Free memory — no longer needed
        self._insert_sizes.clear()

    # ------------------------------------------------------------------
    # Clip pileup construction
    # ------------------------------------------------------------------

    def _build_clip_pileups(self) -> None:
        """Group clipped reads into pileups within +/-5 bp, requiring depth >= 4.

        Clipped reads at approximately the same genomic position are merged
        into a single pileup entry.  Only pileups meeting the minimum depth
        threshold are retained.
        """
        # Group clips by (chrom, clip_side) and sort by position
        groups: dict[tuple[str, str], list[SVRead]] = defaultdict(list)
        for cr in self._clip_reads:
            if cr.clip_side is not None:
                groups[(cr.chrom, cr.clip_side)].append(cr)

        pileups: list[dict[str, Any]] = []

        for (chrom, side), reads in groups.items():
            reads.sort(key=lambda r: r.pos)

            # Sliding-window merge: reads within +/-pileup_window
            cluster: list[SVRead] = []
            cluster_anchor = 0

            for rd in reads:
                if not cluster:
                    cluster = [rd]
                    cluster_anchor = rd.pos
                elif rd.pos - cluster_anchor <= self.pileup_window * 2:
                    cluster.append(rd)
                else:
                    # Finalize previous cluster
                    if len(cluster) >= self.min_pileup_depth:
                        pileups.append(
                            self._finalize_clip_pileup(chrom, side, cluster)
                        )
                    cluster = [rd]
                    cluster_anchor = rd.pos

            # Finalize last cluster
            if len(cluster) >= self.min_pileup_depth:
                pileups.append(
                    self._finalize_clip_pileup(chrom, side, cluster)
                )

        self._clip_pileups = pileups
        logger.info("Built %d clip pileups from %d clipped reads",
                     len(pileups), len(self._clip_reads))

    def _finalize_clip_pileup(
        self, chrom: str, side: str, cluster: list[SVRead]
    ) -> ClipPileup:
        """Create a ClipPileup from a cluster of nearby clipped reads."""
        positions = [r.pos for r in cluster]
        median_pos = sorted(positions)[len(positions) // 2]
        seqs = [r.clip_seq for r in cluster if r.clip_seq]

        return ClipPileup(
            chrom=chrom,
            pos=median_pos,
            depth=len(cluster),
            clip_seqs=seqs,
            clip_side=side,
        )

    # ------------------------------------------------------------------
    # Telemetry emission
    # ------------------------------------------------------------------

    def _emit_event(self, event: dict[str, Any]) -> None:
        """Send a single event through the callback, if registered."""
        if self.callback is not None:
            try:
                self.callback(event)
            except Exception:
                logger.debug("Callback error for event %s", event.get("type"), exc_info=True)

    def _emit_telemetry(self, now: float) -> None:
        """Emit throughput and progress telemetry."""
        elapsed = now - self._last_telemetry_time
        if elapsed <= 0:
            elapsed = TELEMETRY_INTERVAL_SEC  # avoid division by zero

        reads_delta = self._total_reads - self._telemetry_reads_last
        bytes_delta = self._total_bytes - self._telemetry_bytes_last
        disc_delta = len(self._discordant_reads) - self._telemetry_disc_last
        split_delta = len(self._split_reads) - self._telemetry_split_last

        self._emit_event(
            {
                "type": "scan.throughput",
                "reads_per_sec": int(reads_delta / elapsed),
                "bytes_per_sec": int(bytes_delta / elapsed),
                "discordant_per_sec": int(disc_delta / elapsed),
                "split_per_sec": int(split_delta / elapsed),
            }
        )

        self._telemetry_reads_last = self._total_reads
        self._telemetry_bytes_last = self._total_bytes
        self._telemetry_disc_last = len(self._discordant_reads)
        self._telemetry_split_last = len(self._split_reads)

        # Emit per-chromosome progress for the chromosomes being scanned
        for chrom, prog in self._chrom_progress.items():
            if prog.status == "scanning":
                self._emit_chrom_progress(chrom)

        # Flush dirty bins
        self._flush_dirty_bins()

        # Emit pair density deltas
        self._emit_pair_density()

        # Emit provisional top bins every ~2s
        if now - self._last_topbins_time >= 2.0:
            self._emit_top_bins()
            self._last_topbins_time = now

    def _emit_chrom_progress(self, chrom: str) -> None:
        """Emit a chrom.progress event for a single chromosome."""
        prog = self._chrom_progress.get(chrom)
        if prog is None:
            return
        self._emit_event(
            {
                "type": "chrom.progress",
                "chrom": chrom,
                "pct": round(prog.pct, 1),
                "reads": prog.reads_processed,
                "discordant": prog.discordant_count,
                "split": prog.split_count,
            }
        )

    def _flush_dirty_bins(self) -> None:
        """Emit chrom.bin_update events for all bins that changed since last flush."""
        for chrom, bin_idx in self._dirty_bins:
            counts = self._bin_counts[(chrom, bin_idx)]
            bin_start = bin_idx * BIN_SIZE
            bin_end = bin_start + BIN_SIZE
            self._emit_event(
                {
                    "type": "chrom.bin_update",
                    "chrom": chrom,
                    "bin_start": bin_start,
                    "bin_end": bin_end,
                    "discordant": counts["discordant"],
                    "split": counts["split"],
                    "clip": counts["clip"],
                }
            )
        self._dirty_bins.clear()

    def _emit_pair_density(self) -> None:
        """Emit pair.density events for chromosome pairs that changed."""
        for pair_key, count in self._pair_density.items():
            last = self._pair_density_last.get(pair_key, 0)
            if count != last:
                delta = count - last
                self._emit_event(
                    {
                        "type": "pair.density",
                        "chrom_a": pair_key[0],
                        "chrom_b": pair_key[1],
                        "count": count,
                        "delta": delta,
                    }
                )
                self._pair_density_last[pair_key] = count

    def _emit_top_bins(self) -> None:
        """Emit provisional.top_bins event with the top 12 chromosome pairs."""
        if not self._pair_density:
            return
        top = sorted(
            self._pair_density.items(),
            key=lambda x: x[1],
            reverse=True,
        )[:12]
        bins = [
            {
                "chrom_a": pk[0],
                "chrom_b": pk[1],
                "count": cnt,
                "pos_a": 0,
                "pos_b": 0,
            }
            for pk, cnt in top
            if cnt >= 2
        ]
        if bins:
            self._emit_event({
                "type": "provisional.top_bins",
                "bins": bins,
            })
