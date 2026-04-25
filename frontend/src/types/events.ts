/** SSE event type definitions per spec §16 */

export interface ScanStartedEvent {
  type: 'scan.started';
  job_id: string;
  file_path: string;
}

export interface StageChangedEvent {
  type: 'scan.stage_changed';
  stage: string;
}

export interface ScanProgressEvent {
  type: 'scan.progress';
  stage: string;
  pct: number;
  reads_processed?: number;
  bytes_processed?: number;
  detail?: string;
}

export interface ThroughputEvent {
  type: 'scan.throughput';
  reads_per_sec: number;
  bytes_per_sec: number;
  discordant_per_sec: number;
  split_per_sec: number;
}

export interface ChromProgressEvent {
  type: 'chrom.progress';
  chrom: string;
  pct: number;
  reads: number;
  discordant: number;
  split: number;
}

export interface ChromBinUpdateEvent {
  type: 'chrom.bin_update';
  chrom: string;
  bin_start: number;
  bin_end: number;
  discordant: number;
  split: number;
  clip: number;
}

export interface PairDensityEvent {
  type: 'pair.density';
  chrom_a: string;
  chrom_b: string;
  count: number;
  delta: number;
}

export interface EvidenceHighlightEvent {
  type: 'evidence.highlight';
  evidence_type: 'split' | 'clip_pileup' | 'new_pair_cluster' | 'discordant';
  chrom_a: string;
  pos_a: number;
  chrom_b: string;
  pos_b: number;
  summary?: string;
}

export interface ProvisionalTopBinsEvent {
  type: 'provisional.top_bins';
  bins: Array<{
    chrom_a: string;
    chrom_b: string;
    count: number;
    pos_a: number;
    pos_b: number;
  }>;
}

export interface ScanCompletedEvent {
  type: 'scan.completed';
  job_id?: string;
}

export interface ValidationStartedEvent {
  type: 'validation.started';
  job_id?: string;
}

export interface ValidationCallEmittedEvent {
  type: 'validation.call_emitted';
  call: ValidatedCall;
}

export interface ValidationCompletedEvent {
  type: 'validation.completed';
  job_id?: string;
  num_calls: number;
}

export interface ScanCancelledEvent {
  type: 'scan.cancelled';
  job_id?: string;
}

export interface ErrorEvent {
  type: 'error';
  stage: string;
  message: string;
}

export interface ChromPairMatrixEvent {
  type: 'chrom_pair.matrix';
  job_id: string;
  matrix: number[][];
  chrom_order: string[];
  top_pairs: Array<{ chrom_a: string; chrom_b: string; count: number }>;
}

export interface ClusteringScaleEvent {
  type: 'clustering.scale';
  chrom_a: string;
  chrom_b: string;
  scale_bp: number;
  surviving_reads: number;
  reads_before: number;
}

export type SSEEvent =
  | ScanStartedEvent
  | StageChangedEvent
  | ScanProgressEvent
  | ThroughputEvent
  | ChromProgressEvent
  | ChromBinUpdateEvent
  | PairDensityEvent
  | EvidenceHighlightEvent
  | ProvisionalTopBinsEvent
  | ScanCompletedEvent
  | ValidationStartedEvent
  | ValidationCallEmittedEvent
  | ValidationCompletedEvent
  | ScanCancelledEvent
  | ErrorEvent
  | ChromPairMatrixEvent
  | ClusteringScaleEvent;

export interface ValidatedCall {
  event_id: string;
  cluster_id: string;
  chrom_a: string;
  pos_a: number;
  ci_a: [number, number];
  chrom_b: string;
  pos_b: number;
  ci_b: [number, number];
  orientation: string;
  support: {
    discordant: number;
    split: number;
    clipped: number;
    total: number;
  };
  unique_starts_a: number;
  unique_starts_b: number;
  duplicate_fraction: number;
  reciprocal_support: number;
  background_p: number;
  median_mapq: number;
  median_mapq_a: number;
  median_mapq_b: number;
  local_nb_pvalue_a: number;
  local_nb_pvalue_b: number;
  local_rate_a: number;
  local_rate_b: number;
  local_coverage_ratio_a: number;
  local_coverage_ratio_b: number;
  mask_overlaps_a: string[];
  mask_overlaps_b: string[];
  segdup_pct_a: number;
  segdup_pct_b: number;
  promiscuous_hotspot: boolean;
  both_flanks_remap_uniquely: boolean;
  chrom_pair_enrichment: number;
  orientation_distribution: Record<string, number>;
  external_callers: string[];
  assembly_resolved: boolean;
  assembly_bp_offset: number;
  tier: 'confirmed' | 'validated' | 'likely' | 'strong_candidate' | 'candidate' | 'filtered';
  score: number;
  score_components: Record<string, number>;
  microhomology: string;
  inserted_seq: string;
  filter_flags: string[];
  reject_reasons: string[];
  evidence_label: string;
  merged_subclusters: string[];
}

export interface ServerFile {
  path: string;
  name: string;
  dir: string;
  size: number;
  size_human: string;
  format: string;
  indexed: boolean;
}

export interface Reference {
  name: string;
  path: string;
}

export type ScanMode = 'idle' | 'streaming' | 'validated';

export interface ChromData {
  chrom: string;
  pct: number;
  reads: number;
  discordant: number;
  split: number;
  status: 'pending' | 'scanning' | 'complete';
  bins: Map<number, { discordant: number; split: number; clip: number }>;
}

export interface WaterfallEntry {
  id: number;
  timestamp: number;
  type: string;
  chrom_a: string;
  pos_a: number;
  chrom_b: string;
  pos_b: number;
  detail: string;
}
