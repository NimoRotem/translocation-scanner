"""Data models for the translocation scanner pipeline."""
from __future__ import annotations
import time
import uuid
from dataclasses import dataclass, field
from enum import Enum
from typing import NamedTuple, Optional


class EvidenceType(str, Enum):
    DISCORDANT = "discordant"
    SPLIT = "split"
    CLIPPED = "clipped"


class Tier(str, Enum):
    CONFIRMED = "confirmed"
    VALIDATED = "validated"
    LIKELY = "likely"
    STRONG_CANDIDATE = "strong_candidate"
    CANDIDATE = "candidate"
    FILTERED = "filtered"


class ScanStage(str, Enum):
    QUEUED = "queued"
    LIBRARY_STATS = "library_stats"
    EXTRACTION = "extraction"
    CLUSTERING = "clustering"
    CLIP_REALIGNMENT = "clip_realignment"
    EXTERNAL_CALLERS = "external_callers"
    BACKGROUND_MODEL = "background_model"
    FILTERING = "filtering"
    SCORING = "scoring"
    OUTPUT = "output"
    COMPLETED = "completed"
    FAILED = "failed"


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class SVRead:
    """A single SV-relevant read extracted from BAM."""
    read_name_hash: str
    chrom: str
    pos: int
    mapq: int
    evidence_type: EvidenceType
    mate_chrom: Optional[str] = None
    mate_pos: Optional[int] = None
    mate_mapq: Optional[int] = None
    sa_tag: Optional[str] = None
    clip_seq: Optional[str] = None
    clip_side: Optional[str] = None  # "left" or "right"
    clip_len: int = 0
    is_reverse: bool = False
    mate_is_reverse: bool = False
    insert_size: int = 0
    flag: int = 0


class CompactEvidence(NamedTuple):
    """Lightweight per-read record for filter and audit use.

    Stored in EvidenceCluster.reads instead of full SVRead objects.
    Fields match SVRead attribute names so downstream code (e.g.
    filters.py accessing r.is_reverse) works unchanged.
    """
    pos_a: int
    pos_b: int
    mapq: int
    is_reverse: bool        # strand on side A
    mate_is_reverse: bool   # strand on side B
    evidence_type: int      # 0=DISCORDANT, 1=SPLIT, 2=CLIPPED
    flag: int               # BAM flag (bit 0x400 = duplicate)
    read_name_hash: int     # first 12 hex chars of MD5 → uint64


@dataclass
class ClipPileup:
    """A cluster of soft-clipped reads at approximately the same position."""
    chrom: str
    pos: int
    depth: int
    clip_seqs: list[str] = field(default_factory=list)
    clip_side: str = "right"
    partner_chrom: Optional[str] = None
    partner_pos: Optional[int] = None


@dataclass
class EvidenceCluster:
    """A cluster of SV evidence supporting a potential translocation."""
    cluster_id: str = ""
    chrom_a: str = ""
    chrom_b: str = ""
    pos_a: int = 0
    pos_b: int = 0
    ci_a: tuple[int, int] = (-500, 500)
    ci_b: tuple[int, int] = (-500, 500)
    orientation: str = "++"
    discordant_count: int = 0
    split_count: int = 0
    clipped_count: int = 0
    reciprocal_support: int = 0
    background_p: float = 1.0
    median_mapq: float = 0.0
    tier: Tier = Tier.CANDIDATE
    score: float = 0.0
    reads: list[SVRead] = field(default_factory=list)
    microhomology: str = ""
    inserted_seq: str = ""
    filter_flags: list[str] = field(default_factory=list)
    unique_starts_a: int = 0
    unique_starts_b: int = 0
    score_components: dict = field(default_factory=dict)
    reject_reasons: list[str] = field(default_factory=list)
    evidence_label: str = ""
    merged_subclusters: list[str] = field(default_factory=list)
    # Per-side local NB p-values (Step 3)
    local_nb_pvalue_a: float = 1.0
    local_nb_pvalue_b: float = 1.0
    local_rate_a: float = 0.0
    local_rate_b: float = 0.0
    # Per-side median MAPQ
    median_mapq_a: float = 0.0
    median_mapq_b: float = 0.0
    # Coverage ratios
    local_coverage_ratio_a: float = 0.0
    local_coverage_ratio_b: float = 0.0
    # Duplicate fraction
    duplicate_fraction: float = 0.0
    # Mask overlaps
    mask_overlaps_a: list[str] = field(default_factory=list)
    mask_overlaps_b: list[str] = field(default_factory=list)
    segdup_pct_a: float = 0.0
    segdup_pct_b: float = 0.0
    # Promiscuous hotspot
    promiscuous_hotspot: bool = False
    # Flank remap uniqueness
    both_flanks_remap_uniquely: bool = True
    # External caller agreement
    external_callers: list[str] = field(default_factory=list)
    # Assembly status
    assembly_resolved: bool = False
    assembly_bp_offset: int = 0
    # IGV manual confirmation
    igv_confirmed: bool = False
    # Chromosome pair enrichment
    chrom_pair_enrichment: float = 0.0
    # Orientation distribution
    orientation_distribution: dict = field(default_factory=dict)

    @property
    def event_id(self) -> str:
        return f"TRA_{self.chrom_a}_{self.pos_a}_{self.chrom_b}_{self.pos_b}_{self.orientation}"

    @property
    def total_support(self) -> int:
        return self.discordant_count + self.split_count + self.clipped_count

    def to_dict(self) -> dict:
        return {
            "event_id": self.event_id,
            "cluster_id": self.cluster_id,
            "chrom_a": self.chrom_a,
            "pos_a": self.pos_a,
            "ci_a": list(self.ci_a),
            "chrom_b": self.chrom_b,
            "pos_b": self.pos_b,
            "ci_b": list(self.ci_b),
            "orientation": self.orientation,
            "support": {
                "discordant": self.discordant_count,
                "split": self.split_count,
                "clipped": self.clipped_count,
                "total": self.total_support,
            },
            "unique_starts_a": self.unique_starts_a,
            "unique_starts_b": self.unique_starts_b,
            "duplicate_fraction": self.duplicate_fraction,
            "reciprocal_support": self.reciprocal_support,
            "median_mapq": self.median_mapq,
            "median_mapq_a": self.median_mapq_a,
            "median_mapq_b": self.median_mapq_b,
            "local_nb_pvalue_a": self.local_nb_pvalue_a,
            "local_nb_pvalue_b": self.local_nb_pvalue_b,
            "local_rate_a": self.local_rate_a,
            "local_rate_b": self.local_rate_b,
            "local_coverage_ratio_a": self.local_coverage_ratio_a,
            "local_coverage_ratio_b": self.local_coverage_ratio_b,
            "mask_overlaps_a": self.mask_overlaps_a,
            "mask_overlaps_b": self.mask_overlaps_b,
            "segdup_pct_a": self.segdup_pct_a,
            "segdup_pct_b": self.segdup_pct_b,
            "promiscuous_hotspot": self.promiscuous_hotspot,
            "both_flanks_remap_uniquely": self.both_flanks_remap_uniquely,
            "chrom_pair_enrichment": self.chrom_pair_enrichment,
            "orientation_distribution": self.orientation_distribution,
            "external_callers": self.external_callers,
            "assembly_resolved": self.assembly_resolved,
            "assembly_bp_offset": self.assembly_bp_offset,
            "tier": self.tier.value,
            "score": self.score,
            "score_components": self.score_components,
            "filter_flags": self.filter_flags,
            "reject_reasons": self.reject_reasons,
            "evidence_label": self.evidence_label,
            "merged_subclusters": self.merged_subclusters,
        }


@dataclass
class ChromProgress:
    """Per-chromosome scan progress."""
    chrom: str
    length: int = 0
    reads_processed: int = 0
    bytes_processed: int = 0
    discordant_count: int = 0
    split_count: int = 0
    clip_count: int = 0
    pct: float = 0.0
    status: str = "pending"  # pending, scanning, complete


@dataclass
class ScanJob:
    """A translocation scan job."""
    job_id: str = field(default_factory=lambda: str(uuid.uuid4())[:8])
    file_path: str = ""
    reference_path: str = ""
    reference_build: str = "GRCh38"
    status: JobStatus = JobStatus.QUEUED
    stage: ScanStage = ScanStage.QUEUED
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    error: Optional[str] = None
    # Runtime stats
    total_reads: int = 0
    reads_processed: int = 0
    bytes_processed: int = 0
    discordant_count: int = 0
    split_count: int = 0
    clip_count: int = 0
    chimeric_rate: float = 0.0
    insert_size_median: float = 0.0
    insert_size_std: float = 0.0
    # Settings
    settings: dict = field(default_factory=dict)
    # Results
    validated_calls: list[dict] = field(default_factory=list)
    chrom_progress: dict[str, dict] = field(default_factory=dict)
    results_dir: str = ""

    def to_dict(self) -> dict:
        elapsed = 0
        if self.started_at:
            end = self.completed_at or time.time()
            elapsed = end - self.started_at
        return {
            "job_id": self.job_id,
            "file_path": self.file_path,
            "reference_build": self.reference_build,
            "status": self.status.value,
            "stage": self.stage.value,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "completed_at": self.completed_at,
            "elapsed": round(elapsed, 1),
            "error": self.error,
            "total_reads": self.total_reads,
            "reads_processed": self.reads_processed,
            "bytes_processed": self.bytes_processed,
            "discordant_count": self.discordant_count,
            "split_count": self.split_count,
            "clip_count": self.clip_count,
            "settings": self.settings,
            "chimeric_rate": self.chimeric_rate,
            "insert_size_median": self.insert_size_median,
            "insert_size_std": self.insert_size_std,
            "validated_calls": self.validated_calls,
            "num_calls": len(self.validated_calls),
        }
