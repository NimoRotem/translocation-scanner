import { useState, useEffect, useRef } from 'react';
import type { ServerFile, Reference } from '../types/events';

const BASE = '/translocation-scanner';

interface ScanSetupProps {
  onStart: (jobId: string, filePath: string) => void;
}

interface Settings {
  min_mapq: number;
  min_clip_length: number;
  min_split_aligned: number;
  min_pileup_depth: number;
  pileup_window: number;
  merge_distance: number;
  bg_bin_size: number;
  bg_pvalue_threshold: number;
  centromere_margin: number;
  skip_clip_realignment: boolean;
  skip_external_callers: boolean;
  parallel_extraction: boolean;
  num_workers: number;
}

const DEFAULTS: Settings = {
  min_mapq: 20,
  min_clip_length: 20,
  min_split_aligned: 20,
  min_pileup_depth: 4,
  pileup_window: 5,
  merge_distance: 500,
  bg_bin_size: 100_000,
  bg_pvalue_threshold: 0.001,
  centromere_margin: 1_000_000,
  skip_clip_realignment: false,
  skip_external_callers: false,
  parallel_extraction: true,
  num_workers: 0,
};

interface SettingDef {
  key: keyof Settings;
  label: string;
  tooltip: string;
  type: 'number' | 'boolean';
  min?: number;
  max?: number;
  step?: number;
  unit?: string;
  category: 'extraction' | 'clustering' | 'model';
}

const SETTING_DEFS: SettingDef[] = [
  {
    key: 'min_mapq', label: 'Min MAPQ', type: 'number',
    min: 0, max: 60, step: 1,
    tooltip: 'Minimum mapping quality for a read to be considered. Higher values increase confidence but may miss reads in repetitive regions. Default 20 is standard for short-read WGS.',
    category: 'extraction',
  },
  {
    key: 'min_clip_length', label: 'Min Clip Length', type: 'number',
    min: 5, max: 100, step: 5, unit: 'bp',
    tooltip: 'Minimum soft-clipped bases for a read to count as clipped evidence. Lower values find more breakpoints but increase false positives from alignment artifacts.',
    category: 'extraction',
  },
  {
    key: 'min_split_aligned', label: 'Min Split Aligned', type: 'number',
    min: 10, max: 100, step: 5, unit: 'bp',
    tooltip: 'Minimum aligned bases on each side of a split read (SA tag). Ensures both halves of a chimeric alignment are confidently mapped.',
    category: 'extraction',
  },
  {
    key: 'min_pileup_depth', label: 'Min Pileup Depth', type: 'number',
    min: 2, max: 20, step: 1, unit: 'reads',
    tooltip: 'Minimum number of clipped reads at a position to form a clip pileup. Lower values are more sensitive but noisier. Increase for high-coverage samples.',
    category: 'extraction',
  },
  {
    key: 'pileup_window', label: 'Pileup Window', type: 'number',
    min: 1, max: 50, step: 1, unit: 'bp',
    tooltip: 'Radius in base-pairs for grouping nearby clipped reads into a single pileup. Wider windows merge more aggressively; narrower windows preserve precision.',
    category: 'extraction',
  },
  {
    key: 'merge_distance', label: 'Merge Distance', type: 'number',
    min: 100, max: 5000, step: 100, unit: 'bp',
    tooltip: 'Maximum distance between SV reads to merge into the same cluster. Larger values capture more evidence per cluster but may merge distinct breakpoints. Tuned for ~500bp insert size libraries.',
    category: 'clustering',
  },
  {
    key: 'bg_bin_size', label: 'Background Bin Size', type: 'number',
    min: 10_000, max: 1_000_000, step: 10_000, unit: 'bp',
    tooltip: 'Bin size for the Poisson background chimerism model. Smaller bins increase resolution but need more reads per bin for stable estimates. 100kb works well for 30x WGS.',
    category: 'model',
  },
  {
    key: 'bg_pvalue_threshold', label: 'P-value Threshold', type: 'number',
    min: 0.0001, max: 0.05, step: 0.0001,
    tooltip: 'BH-corrected p-value threshold for the background model. Clusters above this threshold are filtered out. Lower values are more stringent; higher values allow more candidates through.',
    category: 'model',
  },
  {
    key: 'centromere_margin', label: 'Centromere Margin', type: 'number',
    min: 100_000, max: 5_000_000, step: 100_000, unit: 'bp',
    tooltip: 'Distance from centromere midpoint within which breakpoints are flagged. Centromeric regions are enriched for alignment artifacts. Wider margins filter more aggressively.',
    category: 'model',
  },
  {
    key: 'skip_clip_realignment', label: 'Skip Clip Realignment', type: 'boolean',
    tooltip: 'Skip the minimap2 clip realignment stage. Disabling saves time but loses the ability to identify partner loci from soft-clipped sequences. Recommended to keep enabled when minimap2 is available.',
    category: 'clustering',
  },
  {
    key: 'skip_external_callers', label: 'Skip DELLY (Fast Mode)', type: 'boolean',
    tooltip: 'Skip the external DELLY caller. DELLY on a full 30x WGS BAM takes ~40 min. Skipping gives results in ~7 min using blind discovery only.',
    category: 'clustering',
  },
  {
    key: 'parallel_extraction', label: 'Parallel Extraction', type: 'boolean',
    tooltip: 'Run each chromosome on a separate CPU core. Dramatically speeds up extraction on multi-core systems (e.g. 97 min → 5 min on 24 cores).',
    category: 'extraction',
  },
  {
    key: 'num_workers', label: 'Worker Count', type: 'number',
    min: 0, max: 32, step: 1,
    tooltip: 'Number of parallel processes for extraction. 0 = auto-detect (one per chromosome, max 24). Only used when Parallel Extraction is enabled.',
    category: 'extraction',
  },
];

function formatValue(def: SettingDef, value: number | boolean): string {
  if (def.type === 'boolean') return value ? 'Yes' : 'No';
  const n = value as number;
  if (def.key === 'bg_pvalue_threshold') return n.toFixed(4);
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + 'M';
  if (n >= 1_000) return (n / 1_000).toFixed(0) + 'k';
  return String(n);
}

function Tooltip({ text, targetRef }: { text: string; targetRef: React.RefObject<HTMLElement | null> }) {
  const [pos, setPos] = useState<{ top: number; left: number } | null>(null);

  useEffect(() => {
    if (targetRef.current) {
      const r = targetRef.current.getBoundingClientRect();
      setPos({ top: r.top - 8, left: r.left + r.width / 2 });
    }
  }, [targetRef]);

  if (!pos) return null;

  return (
    <div
      className="setting-tooltip"
      style={{ top: pos.top, left: pos.left }}
    >
      {text}
    </div>
  );
}

export default function ScanSetup({ onStart }: ScanSetupProps) {
  const [files, setFiles] = useState<ServerFile[]>([]);
  const [refs, setRefs] = useState<Reference[]>([]);
  const [selected, setSelected] = useState('');
  const [refBuild, setRefBuild] = useState('GRCh38');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showSettings, setShowSettings] = useState(false);
  const [settings, setSettings] = useState<Settings>({ ...DEFAULTS });
  const [hoveredSetting, setHoveredSetting] = useState<string | null>(null);
  const hoverRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    fetch(`${BASE}/api/server-files`)
      .then(r => r.json())
      .then(data => {
        setFiles(data.files || []);
        setRefs(data.references || []);
        if (data.files?.length > 0) {
          setSelected(data.files[0].path);
        }
      })
      .catch(() => setError('Failed to load server files'));
  }, []);

  const handleStart = async () => {
    if (!selected) return;
    setLoading(true);
    setError(null);
    try {
      const resp = await fetch(`${BASE}/api/scan`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          file_path: selected,
          reference_build: refBuild,
          settings,
        }),
      });
      if (!resp.ok) {
        const data = await resp.json();
        throw new Error(data.detail || 'Failed to start scan');
      }
      const data = await resp.json();
      onStart(data.job_id, selected);
    } catch (e: any) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  };

  const updateSetting = <K extends keyof Settings>(key: K, value: Settings[K]) => {
    setSettings(prev => ({ ...prev, [key]: value }));
  };

  const resetSettings = () => setSettings({ ...DEFAULTS });

  const isModified = JSON.stringify(settings) !== JSON.stringify(DEFAULTS);

  const categories = [
    { id: 'extraction', label: 'Extraction' },
    { id: 'clustering', label: 'Clustering' },
    { id: 'model', label: 'Background Model' },
  ];

  return (
    <div className="scan-setup">
      <div className="setup-hero">
        <h1>Translocation Scanner</h1>
        <p className="setup-subtitle">
          Interchromosomal breakpoint detection with streaming visualization
        </p>
      </div>

      <div className="setup-form">
        <div className="setup-field">
          <label>BAM/CRAM File</label>
          <select
            value={selected}
            onChange={e => setSelected(e.target.value)}
            className="setup-select"
          >
            {files.length === 0 && <option value="">No files found</option>}
            {files.map(f => (
              <option key={f.path} value={f.path}>
                {f.name} — {f.size_human} {f.format}
                {!f.indexed ? ' (no index!)' : ''}
              </option>
            ))}
          </select>
        </div>

        <div className="setup-field">
          <label>Reference Build</label>
          <select
            value={refBuild}
            onChange={e => setRefBuild(e.target.value)}
            className="setup-select"
          >
            <option value="GRCh38">GRCh38 (hg38)</option>
            <option value="GRCh38_numeric">GRCh38 (numeric chroms)</option>
          </select>
          <div className="setup-hint">
            Available references: {refs.map(r => r.name).join(', ') || 'none detected'}
          </div>
        </div>

        {/* Settings toggle */}
        <div className="settings-toggle-row">
          <button
            className={`settings-toggle-btn ${showSettings ? 'active' : ''}`}
            onClick={() => setShowSettings(!showSettings)}
          >
            <span className="settings-icon">&#9881;</span>
            Advanced Settings
            <span className="settings-chevron">{showSettings ? '▲' : '▼'}</span>
          </button>
          {isModified && !showSettings && (
            <span className="settings-modified-badge">Modified</span>
          )}
        </div>

        {/* Settings panel */}
        {showSettings && (
          <div className="settings-panel">
            {categories.map(cat => {
              const defs = SETTING_DEFS.filter(d => d.category === cat.id);
              return (
                <div key={cat.id} className="settings-category">
                  <h4 className="settings-cat-label">{cat.label}</h4>
                  <div className="settings-grid">
                    {defs.map(def => {
                      const val = settings[def.key];
                      return (
                        <div
                          key={def.key}
                          className={`setting-row ${hoveredSetting === def.key ? 'hovered' : ''}`}
                          onMouseEnter={(e) => {
                            setHoveredSetting(def.key);
                            hoverRef.current = e.currentTarget as HTMLDivElement;
                          }}
                          onMouseLeave={() => setHoveredSetting(null)}
                        >
                          <div className="setting-label-row">
                            <span className="setting-label">
                              {def.label}
                              <span className="setting-info-icon" title={def.tooltip}>?</span>
                            </span>
                            <span className="setting-value-display">
                              {formatValue(def, val)}
                              {def.unit && <span className="setting-unit">{def.unit}</span>}
                            </span>
                          </div>
                          {def.type === 'boolean' ? (
                            <label className="setting-toggle">
                              <input
                                type="checkbox"
                                checked={val as boolean}
                                onChange={e => updateSetting(def.key, e.target.checked as any)}
                              />
                              <span className="toggle-track">
                                <span className="toggle-thumb" />
                              </span>
                            </label>
                          ) : (
                            <input
                              type="range"
                              className="setting-slider"
                              min={def.min}
                              max={def.max}
                              step={def.step}
                              value={val as number}
                              onChange={e => updateSetting(def.key, Number(e.target.value) as any)}
                            />
                          )}
                          {hoveredSetting === def.key && (
                            <div className="setting-tooltip-inline">
                              {def.tooltip}
                            </div>
                          )}
                        </div>
                      );
                    })}
                  </div>
                </div>
              );
            })}
            <div className="settings-actions">
              {isModified && (
                <button className="settings-reset-btn" onClick={resetSettings}>
                  Reset to Defaults
                </button>
              )}
            </div>
          </div>
        )}

        {error && <div className="setup-error">{error}</div>}

        <button
          className="setup-start-btn"
          onClick={handleStart}
          disabled={!selected || loading}
        >
          {loading ? 'Starting...' : 'Start Scan'}
        </button>

        <div className="setup-info">
          <h3>How it works</h3>
          <ol>
            <li>Single-pass BAM extraction identifies discordant pairs, split reads, and clipped pileups</li>
            <li>Breakpoint clustering groups evidence by chromosome pair and orientation</li>
            <li>Background chimerism model separates signal from library noise</li>
            <li>Hard/soft filters remove centromeric, low-MAPQ, and artifact clusters</li>
            <li>Validated calls are scored and tiered (confirmed / likely / candidate)</li>
          </ol>
          <p>The streaming dashboard shows live extraction progress — provisional visualizations that resolve into a calibrated callset.</p>
        </div>
      </div>
    </div>
  );
}
