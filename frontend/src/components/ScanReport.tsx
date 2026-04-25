import { useState, useEffect } from 'react';
import type { ValidatedCall } from '../types/events';
import ValidatedCallset from './ValidatedCallset';

const BASE = '/translocation-scanner';

interface ReportData {
  sample: {
    name: string;
    path: string;
    reference_build: string;
    scan_date: number | null;
    elapsed_seconds: number;
  };
  quality: {
    total_reads: number;
    chimeric_rate: number;
    chimeric_rate_pct: string;
    chimeric_assessment: 'normal' | 'elevated' | 'high';
    insert_size_median: number;
    insert_size_std: number;
    insert_size_assessment: string;
  };
  evidence: {
    discordant: number;
    split: number;
    clip_pileups: number;
  };
  pipeline: {
    clusters_formed: number;
    clusters_passing: number;
    timings: Record<string, number>;
    filter_breakdown: Record<string, number>;
  };
  results: {
    total_calls: number;
    by_tier: Record<string, number>;
    filtered: number;
    calls: ValidatedCall[];
  };
  interpretation: {
    summary: string;
    detail: string;
  };
  mask_manifest_version?: string;
  warnings?: string[];
}

interface ScanReportProps {
  jobId: string | null;
  validatedCalls: ValidatedCall[];
}

const ASSESSMENT_COLORS: Record<string, string> = {
  normal: '#10b981',
  elevated: '#f59e0b',
  high: '#ef4444',
  atypical: '#f59e0b',
};

function formatNum(n: number): string {
  if (n >= 1_000_000_000) return (n / 1_000_000_000).toFixed(1) + 'B';
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + 'M';
  if (n >= 1_000) return (n / 1_000).toFixed(1) + 'k';
  return String(n);
}

function formatDuration(sec: number): string {
  if (!sec) return '--';
  const m = Math.floor(sec / 60);
  const s = Math.floor(sec % 60);
  if (m > 0) return `${m}m ${s}s`;
  return `${s}s`;
}

function formatDate(ts: number | null): string {
  if (!ts) return '--';
  const d = new Date(ts * 1000);
  return d.toLocaleDateString('en-US', {
    year: 'numeric',
    month: 'long',
    day: 'numeric',
    hour: '2-digit',
    minute: '2-digit',
  });
}

const FILTER_LABELS: Record<string, string> = {
  centromere_proximity: 'Near centromere',
  acrocentric_parm: 'Acrocentric p-arm',
  telomere_proximity: 'Near telomere',
  blacklist: 'ENCODE blacklist',
  segdup: 'Segmental duplication',
  low_mapq: 'Low mapping quality',
  high_background_p: 'Background noise',
  insufficient_support: 'Insufficient support',
  single_orientation: 'Single strand orientation',
  hard_exclude_MT: 'Mitochondrial',
  hard_exclude_non_primary: 'Non-primary chromosome',
  hard_exclude_blacklist: 'Blacklist core overlap',
  hard_exclude_centromere: 'Centromere core overlap',
  hard_exclude_telomere: 'Telomere core overlap',
  hard_exclude_orientation: 'Orientation incoherent',
  reject_sr0_low_pr: 'SR=0, PR<10',
  reject_high_dup: 'High duplicate fraction',
  reject_low_unique_starts: 'Low unique starts',
  reject_high_coverage_ratio: 'High coverage ratio',
  reject_promiscuous: 'Promiscuous hotspot',
};

export default function ScanReport({ jobId, validatedCalls }: ScanReportProps) {
  const [report, setReport] = useState<ReportData | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!jobId) return;
    fetch(`${BASE}/api/jobs/${jobId}/report`)
      .then(r => r.json())
      .then(data => {
        setReport(data);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, [jobId]);

  if (loading) {
    return <div className="report-loading">Loading report...</div>;
  }

  if (!report) {
    return (
      <div className="report-fallback">
        <ValidatedCallset calls={validatedCalls} jobId={jobId} artifactWarning={false} />
      </div>
    );
  }

  const { sample, quality, evidence, pipeline, results, interpretation } = report;
  const isArtifactDominated = (report.warnings || []).some(w => w.includes('artifact'));
  const totalFilterBreakdown = Object.values(pipeline.filter_breakdown).reduce((a, b) => a + b, 0);

  return (
    <div className="scan-report">
      {/* ===== Interpretation Banner ===== */}
      <div className={`report-interpretation ${results.total_calls > 0 ? 'has-calls' : 'no-calls'}`}>
        <div className="interp-icon">
          {results.total_calls > 0 ? (
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none" stroke="#f59e0b" strokeWidth="2">
              <path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/>
              <line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/>
            </svg>
          ) : (
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none" stroke="#10b981" strokeWidth="2">
              <path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/>
              <polyline points="22 4 12 14.01 9 11.01"/>
            </svg>
          )}
        </div>
        <div className="interp-text">
          <div className="interp-summary">{interpretation.summary}</div>
          <div className="interp-detail">{interpretation.detail}</div>
        </div>
      </div>

      {/* ===== Sample & Quality ===== */}
      <div className="report-grid">
        <div className="report-card">
          <div className="report-card-title">Sample</div>
          <div className="report-stats">
            <div className="rstat">
              <span className="rstat-label">File</span>
              <span className="rstat-value mono">{sample.name}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Reference</span>
              <span className="rstat-value">{sample.reference_build}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Scan date</span>
              <span className="rstat-value">{formatDate(sample.scan_date)}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Runtime</span>
              <span className="rstat-value">{formatDuration(sample.elapsed_seconds)}</span>
            </div>
          </div>
        </div>

        <div className="report-card">
          <div className="report-card-title">Quality Metrics</div>
          <div className="report-stats">
            <div className="rstat">
              <span className="rstat-label">Total reads</span>
              <span className="rstat-value">{formatNum(quality.total_reads)}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Chimeric rate</span>
              <span className="rstat-value">
                {quality.chimeric_rate_pct}
                <span className="assessment-badge" style={{ color: ASSESSMENT_COLORS[quality.chimeric_assessment] }}>
                  {quality.chimeric_assessment}
                </span>
              </span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Insert size</span>
              <span className="rstat-value">
                {quality.insert_size_median.toFixed(0)} &plusmn; {quality.insert_size_std.toFixed(0)} bp
                <span className="assessment-badge" style={{ color: ASSESSMENT_COLORS[quality.insert_size_assessment] }}>
                  {quality.insert_size_assessment}
                </span>
              </span>
            </div>
          </div>
        </div>

        <div className="report-card">
          <div className="report-card-title">Evidence Extracted</div>
          <div className="report-stats">
            <div className="rstat">
              <span className="rstat-label">Discordant pairs</span>
              <span className="rstat-value" style={{ color: '#ef4444' }}>{formatNum(evidence.discordant)}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Split reads</span>
              <span className="rstat-value" style={{ color: '#10b981' }}>{formatNum(evidence.split)}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Clip pileups</span>
              <span className="rstat-value" style={{ color: '#3b82f6' }}>{formatNum(evidence.clip_pileups)}</span>
            </div>
          </div>
        </div>

        <div className="report-card">
          <div className="report-card-title">Pipeline Summary</div>
          <div className="report-stats">
            <div className="rstat">
              <span className="rstat-label">Clusters formed</span>
              <span className="rstat-value">{formatNum(pipeline.clusters_formed)}</span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Clusters passing</span>
              <span className="rstat-value" style={{ color: pipeline.clusters_passing > 0 ? '#10b981' : '#6b7280' }}>
                {pipeline.clusters_passing}
              </span>
            </div>
            <div className="rstat">
              <span className="rstat-label">Filtered out</span>
              <span className="rstat-value" style={{ color: '#6b7280' }}>
                {formatNum(results.filtered)}
              </span>
            </div>
          </div>
        </div>
      </div>

      {/* ===== Pipeline Timing ===== */}
      {Object.keys(pipeline.timings).length > 0 && (
        <div className="report-card full-width">
          <div className="report-card-title">Pipeline Timing</div>
          <div className="timing-bar-container">
            {(() => {
              const total = pipeline.timings.total || Object.values(pipeline.timings).reduce((a, b) => a + b, 0);
              const stages = [
                { key: 'extraction', label: 'Extraction', color: '#3b82f6' },
                { key: 'clustering', label: 'Clustering', color: '#f59e0b' },
                { key: 'clip_realignment', label: 'Clip realign', color: '#8b5cf6' },
                { key: 'external_callers', label: 'External callers', color: '#22d3ee' },
                { key: 'background_model', label: 'Background', color: '#10b981' },
                { key: 'filtering', label: 'Filtering', color: '#ef4444' },
                { key: 'scoring', label: 'Scoring', color: '#ec4899' },
                { key: 'output', label: 'Output', color: '#6b7280' },
              ];
              return (
                <>
                  <div className="timing-bar">
                    {stages.map(s => {
                      const t = pipeline.timings[s.key] || 0;
                      const pct = total > 0 ? (t / total) * 100 : 0;
                      if (pct < 0.5) return null;
                      return (
                        <div
                          key={s.key}
                          className="timing-segment"
                          style={{ width: `${pct}%`, backgroundColor: s.color }}
                          title={`${s.label}: ${formatDuration(t)}`}
                        />
                      );
                    })}
                  </div>
                  <div className="timing-legend">
                    {stages.map(s => {
                      const t = pipeline.timings[s.key] || 0;
                      if (t < 0.5) return null;
                      return (
                        <span key={s.key} className="timing-item">
                          <span className="timing-dot" style={{ backgroundColor: s.color }} />
                          {s.label}: {formatDuration(t)}
                        </span>
                      );
                    })}
                    <span className="timing-item total">Total: {formatDuration(total)}</span>
                  </div>
                </>
              );
            })()}
          </div>
        </div>
      )}

      {/* ===== Filter Breakdown ===== */}
      {totalFilterBreakdown > 0 && (
        <div className="report-card full-width">
          <div className="report-card-title">Filter Breakdown</div>
          <div className="filter-bars">
            {Object.entries(pipeline.filter_breakdown)
              .sort(([, a], [, b]) => b - a)
              .map(([flag, count]) => {
                const pct = (count / totalFilterBreakdown) * 100;
                return (
                  <div key={flag} className="filter-row">
                    <span className="filter-label">{FILTER_LABELS[flag] || flag}</span>
                    <div className="filter-bar-track">
                      <div className="filter-bar-fill" style={{ width: `${pct}%` }} />
                    </div>
                    <span className="filter-count">{formatNum(count)}</span>
                  </div>
                );
              })}
          </div>
        </div>
      )}

      {/* ===== Results Table ===== */}
      {results.total_calls > 0 && (
        <div className="report-card full-width">
          <ValidatedCallset calls={validatedCalls} jobId={jobId} artifactWarning={isArtifactDominated} />
        </div>
      )}

      {/* ===== Downloads ===== */}
      {jobId && results.total_calls > 0 && (
        <div className="report-downloads">
          <span className="report-downloads-label">Download results:</span>
          <a href={`${BASE}/api/jobs/${jobId}/download/translocations.vcf`} className="export-btn">VCF</a>
          <a href={`${BASE}/api/jobs/${jobId}/download/translocations.bedpe`} className="export-btn">BEDPE</a>
          <a href={`${BASE}/api/jobs/${jobId}/download/translocations.json`} className="export-btn">JSON</a>
        </div>
      )}
    </div>
  );
}
