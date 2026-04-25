import { useState, useEffect } from 'react';
import type { ValidatedCall } from '../types/events';

const BASE = '/translocation-scanner';

interface JobSummary {
  job_id: string;
  file_path: string;
  reference_build: string;
  status: string;
  stage: string;
  created_at: number;
  started_at: number | null;
  completed_at: number | null;
  elapsed: number;
  error: string | null;
  total_reads: number;
  reads_processed: number;
  bytes_processed: number;
  discordant_count: number;
  split_count: number;
  clip_count: number;
  chimeric_rate: number;
  insert_size_median: number;
  insert_size_std: number;
  num_calls: number;
  validated_calls: ValidatedCall[];
}

interface PreviousRunsProps {
  onViewJob: (job: JobSummary) => void;
  onResumeScan: (jobId: string, filePath: string) => void;
}

const TIER_COLORS: Record<string, string> = {
  confirmed: '#10b981',
  likely: '#3b82f6',
  candidate: '#6b7280',
  filtered: '#ef4444',
};

function formatDate(ts: number): string {
  if (!ts) return '—';
  const d = new Date(ts * 1000);
  return d.toLocaleDateString('en-US', {
    month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit',
  });
}

function formatDuration(sec: number): string {
  if (!sec) return '—';
  const h = Math.floor(sec / 3600);
  const m = Math.floor((sec % 3600) / 60);
  const s = Math.floor(sec % 60);
  if (h > 0) return `${h}h ${m}m ${s}s`;
  if (m > 0) return `${m}m ${s}s`;
  return `${s}s`;
}

function formatNum(n: number): string {
  if (n >= 1_000_000_000) return (n / 1_000_000_000).toFixed(1) + 'B';
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + 'M';
  if (n >= 1_000) return (n / 1_000).toFixed(1) + 'k';
  return String(n);
}

export default function PreviousRuns({ onViewJob, onResumeScan }: PreviousRunsProps) {
  const [jobList, setJobList] = useState<JobSummary[]>([]);
  const [loading, setLoading] = useState(true);
  const [expanded, setExpanded] = useState<string | null>(null);

  useEffect(() => {
    fetch(`${BASE}/api/jobs`)
      .then(r => r.json())
      .then(data => {
        setJobList(data.jobs || []);
        setLoading(false);
      })
      .catch(() => setLoading(false));
  }, []);

  if (loading) {
    return <div className="prev-runs-loading">Loading previous runs...</div>;
  }

  if (jobList.length === 0) {
    return (
      <div className="prev-runs">
        <h2 className="prev-runs-title">Previous Runs</h2>
        <div className="prev-runs-empty">No previous scans yet. Start one above.</div>
      </div>
    );
  }

  return (
    <div className="prev-runs">
      <h2 className="prev-runs-title">Previous Runs</h2>
      <div className="prev-runs-list">
        {jobList.map(job => {
          const fileName = job.file_path.split('/').pop() || job.file_path;
          const isExpanded = expanded === job.job_id;
          const tierCounts: Record<string, number> = {};
          (job.validated_calls || []).forEach(c => {
            tierCounts[c.tier] = (tierCounts[c.tier] || 0) + 1;
          });

          const statusColor = job.status === 'completed' ? '#10b981'
            : job.status === 'failed' ? '#ef4444'
            : job.status === 'running' ? '#f59e0b'
            : '#6b7280';

          return (
            <div key={job.job_id} className="prev-run-card">
              <div
                className="prev-run-header"
                onClick={() => setExpanded(isExpanded ? null : job.job_id)}
              >
                <div className="prev-run-main">
                  <span className="prev-run-file">{fileName}</span>
                  <span className="prev-run-status" style={{ color: statusColor }}>
                    {job.status.toUpperCase()}
                  </span>
                </div>
                <div className="prev-run-meta">
                  <span>{formatDate(job.created_at)}</span>
                  <span>{formatDuration(job.elapsed)}</span>
                  <span>{job.reference_build}</span>
                  {job.num_calls > 0 && (
                    <span className="prev-run-calls">{job.num_calls} calls</span>
                  )}
                  <span className="prev-run-expand">{isExpanded ? '▲' : '▼'}</span>
                </div>
              </div>

              {isExpanded && (
                <div className="prev-run-detail">
                  <div className="prev-run-stats">
                    <div className="prev-stat">
                      <span className="prev-stat-label">Job ID</span>
                      <span className="prev-stat-value mono">{job.job_id}</span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Total Reads</span>
                      <span className="prev-stat-value">{formatNum(job.total_reads)}</span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Discordant</span>
                      <span className="prev-stat-value" style={{ color: '#ef4444' }}>
                        {formatNum(job.discordant_count)}
                      </span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Split</span>
                      <span className="prev-stat-value" style={{ color: '#10b981' }}>
                        {formatNum(job.split_count)}
                      </span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Clip Pileups</span>
                      <span className="prev-stat-value">{formatNum(job.clip_count)}</span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Chimeric Rate</span>
                      <span className="prev-stat-value">
                        {job.chimeric_rate ? (job.chimeric_rate * 100).toFixed(3) + '%' : '—'}
                      </span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">Insert Size</span>
                      <span className="prev-stat-value">
                        {job.insert_size_median ? `${job.insert_size_median.toFixed(0)} ± ${job.insert_size_std.toFixed(0)}` : '—'}
                      </span>
                    </div>
                    <div className="prev-stat">
                      <span className="prev-stat-label">File</span>
                      <span className="prev-stat-value mono" style={{ fontSize: 11 }}>
                        {job.file_path}
                      </span>
                    </div>
                  </div>

                  {job.error && (
                    <div className="prev-run-error">
                      Error: {job.error}
                    </div>
                  )}

                  {Object.keys(tierCounts).length > 0 && (
                    <div className="prev-run-tiers">
                      {Object.entries(tierCounts).map(([tier, count]) => (
                        <span key={tier} className="prev-tier-badge" style={{ color: TIER_COLORS[tier] }}>
                          <span className="tier-dot" style={{ backgroundColor: TIER_COLORS[tier] }} />
                          {count} {tier}
                        </span>
                      ))}
                    </div>
                  )}

                  {(job.validated_calls || []).length > 0 && (
                    <div className="prev-run-calls-table">
                      <table>
                        <thead>
                          <tr>
                            <th>Tier</th>
                            <th>Breakpoint A</th>
                            <th>Breakpoint B</th>
                            <th>Orient</th>
                            <th>Support</th>
                            <th>Score</th>
                            <th>P-value</th>
                          </tr>
                        </thead>
                        <tbody>
                          {job.validated_calls.map((call, i) => (
                            <tr key={i}>
                              <td>
                                <span className="tier-dot" style={{ backgroundColor: TIER_COLORS[call.tier] }} />
                                {call.tier}
                              </td>
                              <td className="mono">{call.chrom_a}:{call.pos_a.toLocaleString()}</td>
                              <td className="mono">{call.chrom_b}:{call.pos_b.toLocaleString()}</td>
                              <td className="mono">{call.orientation}</td>
                              <td>
                                {call.support.discordant}d {call.support.split}s {call.support.clipped}c
                              </td>
                              <td className="mono">{call.score.toFixed(1)}</td>
                              <td className="mono">
                                {call.background_p < 1e-10
                                  ? call.background_p.toExponential(1)
                                  : call.background_p < 0.001
                                  ? call.background_p.toExponential(2)
                                  : call.background_p.toFixed(4)}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  )}

                  <div className="prev-run-actions">
                    {job.status === 'completed' && job.num_calls > 0 && (
                      <>
                        <a href={`${BASE}/api/jobs/${job.job_id}/download/calls.vcf`} className="prev-action-btn">
                          Download VCF
                        </a>
                        <a href={`${BASE}/api/jobs/${job.job_id}/download/calls.bedpe`} className="prev-action-btn">
                          Download BEDPE
                        </a>
                        <a href={`${BASE}/api/jobs/${job.job_id}/download/calls.json`} className="prev-action-btn">
                          Download JSON
                        </a>
                      </>
                    )}
                    {job.status === 'running' && (
                      <button
                        className="prev-action-btn accent"
                        onClick={() => onResumeScan(job.job_id, job.file_path)}
                      >
                        Resume Stream
                      </button>
                    )}
                  </div>
                </div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
