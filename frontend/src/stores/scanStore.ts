import { create } from 'zustand';
import type {
  ScanMode, ValidatedCall, ChromData, WaterfallEntry,
  ThroughputEvent, ChromProgressEvent, ChromBinUpdateEvent,
  PairDensityEvent, EvidenceHighlightEvent, SSEEvent,
  ChromPairMatrixEvent, ClusteringScaleEvent,
} from '../types/events';

const CHROMS = [
  'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
  'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
  'chr20','chr21','chr22','chrX','chrY',
];
const CHROM_INDEX: Record<string, number> = {};
CHROMS.forEach((c, i) => { CHROM_INDEX[c] = i; });

// Also handle non-chr-prefixed names
for (let i = 1; i <= 22; i++) {
  CHROM_INDEX[String(i)] = i - 1;
}
CHROM_INDEX['X'] = 22;
CHROM_INDEX['Y'] = 23;

export function normalizeChrom(c: string): string {
  if (c.startsWith('chr')) return c;
  if (c === 'X' || c === 'Y') return 'chr' + c;
  const n = parseInt(c);
  if (!isNaN(n) && n >= 1 && n <= 22) return 'chr' + c;
  return c;
}

export function chromIndex(c: string): number {
  return CHROM_INDEX[c] ?? CHROM_INDEX[normalizeChrom(c)] ?? -1;
}

interface ThroughputSnapshot {
  t: number;
  reads_per_sec: number;
  bytes_per_sec: number;
  discordant_per_sec: number;
  split_per_sec: number;
}

interface SSEDebugState {
  connected: boolean;
  jobId: string | null;
  lastEventType: string;
  lastEventTime: number;
  eventCount: number;
  errors: Array<{ time: number; message: string }>;
}

interface ScanState {
  // Mode
  mode: ScanMode;
  jobId: string | null;
  stage: string;
  error: string | null;
  startedAt: number | null;

  // Throughput
  throughput: ThroughputSnapshot;
  throughputHistory: ThroughputSnapshot[];

  // Chromosome progress
  chromProgress: Record<string, ChromData>;

  // Density matrix: 24x24
  densityMatrix: number[][];

  // Waterfall
  waterfall: WaterfallEntry[];
  waterfallCounter: number;

  // Provisional arcs
  provisionalArcs: Array<{
    chrom_a: string; pos_a: number;
    chrom_b: string; pos_b: number;
    count: number; timestamp: number;
  }>;

  // Provisional top bins
  topBins: Array<{
    chrom_a: string; chrom_b: string;
    count: number; pos_a: number; pos_b: number;
  }>;

  // Validated calls
  validatedCalls: ValidatedCall[];

  // Overall progress
  overallPct: number;
  readsProcessed: number;
  bytesProcessed: number;

  // SSE debug state
  sseDebug: SSEDebugState;

  // Actions
  handleEvent: (event: SSEEvent) => void;
  reset: () => void;
  startScan: (jobId: string) => void;
  setSSEDebug: (patch: Partial<SSEDebugState>) => void;
}

const initialDensity = (): number[][] =>
  Array.from({ length: 24 }, () => new Array(24).fill(0));

const initialChromProgress = (): Record<string, ChromData> => {
  const m: Record<string, ChromData> = {};
  CHROMS.forEach(c => {
    m[c] = { chrom: c, pct: 0, reads: 0, discordant: 0, split: 0, status: 'pending', bins: new Map() };
  });
  return m;
};

export const useScanStore = create<ScanState>((set, get) => ({
  mode: 'idle',
  jobId: null,
  stage: '',
  error: null,
  startedAt: null,
  throughput: { t: 0, reads_per_sec: 0, bytes_per_sec: 0, discordant_per_sec: 0, split_per_sec: 0 },
  throughputHistory: [],
  chromProgress: initialChromProgress(),
  densityMatrix: initialDensity(),
  waterfall: [],
  waterfallCounter: 0,
  provisionalArcs: [],
  topBins: [],
  validatedCalls: [],
  overallPct: 0,
  readsProcessed: 0,
  bytesProcessed: 0,
  sseDebug: { connected: false, jobId: null, lastEventType: '', lastEventTime: 0, eventCount: 0, errors: [] },

  setSSEDebug: (patch: Partial<SSEDebugState>) => set(s => ({
    sseDebug: { ...s.sseDebug, ...patch },
  })),

  startScan: (jobId: string) => set({
    mode: 'streaming',
    jobId,
    stage: 'queued',
    error: null,
    startedAt: Date.now(),
    throughput: { t: 0, reads_per_sec: 0, bytes_per_sec: 0, discordant_per_sec: 0, split_per_sec: 0 },
    throughputHistory: [],
    chromProgress: initialChromProgress(),
    densityMatrix: initialDensity(),
    waterfall: [],
    waterfallCounter: 0,
    provisionalArcs: [],
    topBins: [],
    validatedCalls: [],
    overallPct: 0,
    readsProcessed: 0,
    bytesProcessed: 0,
  }),

  reset: () => set({
    mode: 'idle',
    jobId: null,
    stage: '',
    error: null,
    startedAt: null,
    throughputHistory: [],
    chromProgress: initialChromProgress(),
    densityMatrix: initialDensity(),
    waterfall: [],
    provisionalArcs: [],
    topBins: [],
    validatedCalls: [],
    overallPct: 0,
  }),

  handleEvent: (event: SSEEvent) => {
    const state = get();
    switch (event.type) {
      case 'scan.started':
        set({ mode: 'streaming', stage: 'extraction' });
        break;

      case 'scan.stage_changed':
        set({ stage: event.stage });
        break;

      case 'scan.progress':
        set({
          overallPct: event.pct ?? state.overallPct,
          readsProcessed: event.reads_processed ?? state.readsProcessed,
          bytesProcessed: event.bytes_processed ?? state.bytesProcessed,
        });
        break;

      case 'scan.throughput': {
        const snap: ThroughputSnapshot = {
          t: Date.now(),
          reads_per_sec: event.reads_per_sec,
          bytes_per_sec: event.bytes_per_sec,
          discordant_per_sec: event.discordant_per_sec,
          split_per_sec: event.split_per_sec,
        };
        const hist = [...state.throughputHistory, snap].slice(-60);
        set({ throughput: snap, throughputHistory: hist });
        break;
      }

      case 'chrom.progress': {
        const nc = normalizeChrom(event.chrom);
        const cp = { ...state.chromProgress };
        const prev = cp[nc] || { chrom: nc, pct: 0, reads: 0, discordant: 0, split: 0, status: 'pending' as const, bins: new Map() };
        cp[nc] = {
          ...prev,
          pct: event.pct,
          reads: event.reads,
          discordant: event.discordant,
          split: event.split,
          status: event.pct >= 100 ? 'complete' : 'scanning',
        };
        set({ chromProgress: cp });
        break;
      }

      case 'chrom.bin_update': {
        const nc = normalizeChrom(event.chrom);
        const cp = { ...state.chromProgress };
        const prev = cp[nc];
        if (prev) {
          const newBins = new Map(prev.bins);
          newBins.set(event.bin_start, {
            discordant: event.discordant,
            split: event.split,
            clip: event.clip,
          });
          cp[nc] = { ...prev, bins: newBins };
          set({ chromProgress: cp });
        }
        break;
      }

      case 'pair.density': {
        const ia = chromIndex(event.chrom_a);
        const ib = chromIndex(event.chrom_b);
        if (ia >= 0 && ia < 24 && ib >= 0 && ib < 24) {
          const dm = state.densityMatrix.map(r => [...r]);
          dm[ia][ib] = event.count;
          dm[ib][ia] = event.count; // symmetric
          set({ densityMatrix: dm });
        }

        // Update provisional arcs
        const arcs = [...state.provisionalArcs];
        const existing = arcs.find(a =>
          a.chrom_a === event.chrom_a && a.chrom_b === event.chrom_b
        );
        if (existing) {
          existing.count = event.count;
          existing.timestamp = Date.now();
        } else if (event.count >= 2) {
          arcs.push({
            chrom_a: event.chrom_a, pos_a: 0,
            chrom_b: event.chrom_b, pos_b: 0,
            count: event.count,
            timestamp: Date.now(),
          });
        }
        // Keep max 200 arcs
        set({ provisionalArcs: arcs.slice(-200) });
        break;
      }

      case 'evidence.highlight': {
        const wf: WaterfallEntry = {
          id: state.waterfallCounter + 1,
          timestamp: Date.now(),
          type: event.evidence_type,
          chrom_a: event.chrom_a,
          pos_a: event.pos_a,
          chrom_b: event.chrom_b,
          pos_b: event.pos_b,
          detail: event.summary || '',
        };
        const entries = [...state.waterfall, wf].slice(-500);
        set({ waterfall: entries, waterfallCounter: wf.id });

        // Also update provisional arcs with position
        const arcs2 = [...state.provisionalArcs];
        const arc = arcs2.find(a =>
          a.chrom_a === event.chrom_a && a.chrom_b === event.chrom_b
        );
        if (arc) {
          arc.pos_a = event.pos_a;
          arc.pos_b = event.pos_b;
          arc.timestamp = Date.now();
          set({ provisionalArcs: arcs2 });
        }
        break;
      }

      case 'provisional.top_bins':
        set({ topBins: (event as any).bins || [] });
        break;

      case 'scan.completed':
        set({ stage: 'validation', overallPct: 100 });
        break;

      case 'validation.started':
        set({ stage: 'validation' });
        break;

      case 'validation.call_emitted': {
        const calls = [...state.validatedCalls, event.call];
        set({ validatedCalls: calls });
        break;
      }

      case 'validation.completed':
        set({ mode: 'validated', stage: 'completed' });
        break;

      case 'scan.cancelled':
        set({ mode: 'idle', stage: 'cancelled', error: 'Scan cancelled' });
        break;

      case 'error':
        set({ error: event.message, stage: 'failed' });
        break;

      case 'chrom_pair.matrix': {
        // Full matrix update from Phase 1.2
        const matrixEvt = event as ChromPairMatrixEvent;
        if (matrixEvt.matrix && matrixEvt.matrix.length === 24) {
          set({ densityMatrix: matrixEvt.matrix });
        }
        break;
      }

      case 'clustering.scale': {
        // Multi-scale clustering progress from Phase 1.3
        // Just update the stage detail — the density matrix already shows the data
        break;
      }
    }
  },
}));

export { CHROMS, CHROM_INDEX };
