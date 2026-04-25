import { useState } from 'react';
import { motion, AnimatePresence } from 'framer-motion';
import type { ValidatedCall } from '../types/events';

interface ValidatedCallsetProps {
  calls: ValidatedCall[];
  jobId: string | null;
  artifactWarning?: boolean;
}

type SortKey = 'tier' | 'score' | 'chrom_a' | 'support' | 'pvalue';

const TIER_ORDER: Record<string, number> = {
  confirmed: 0,
  validated: 1,
  likely: 2,
  strong_candidate: 3,
  candidate: 4,
  filtered: 5,
};

const TIER_COLORS: Record<string, string> = {
  confirmed: '#10b981',
  validated: '#22d3ee',
  likely: '#3b82f6',
  strong_candidate: '#a78bfa',
  candidate: '#6b7280',
  filtered: '#ef4444',
};

const TIER_LABELS: Record<string, string> = {
  confirmed: 'Confirmed',
  validated: 'Validated',
  likely: 'Likely',
  strong_candidate: 'Strong Candidate',
  candidate: 'Candidate',
  filtered: 'Filtered',
};

function formatPos(pos: number): string {
  return pos.toLocaleString();
}

function formatPval(p: number): string {
  if (p === 0) return '0';
  if (p < 1e-10) return p.toExponential(1);
  if (p < 0.001) return p.toExponential(2);
  return p.toFixed(4);
}

function worstPval(call: ValidatedCall): number {
  return Math.max(call.local_nb_pvalue_a ?? 1, call.local_nb_pvalue_b ?? 1);
}

const BASE = '/translocation-scanner';

const HIGH_CONFIDENCE_TIERS = new Set(['confirmed', 'validated', 'likely', 'strong_candidate']);

export default function ValidatedCallset({ calls, jobId, artifactWarning }: ValidatedCallsetProps) {
  const [sortKey, setSortKey] = useState<SortKey>('tier');
  const [sortAsc, setSortAsc] = useState(true);
  const [selected, setSelected] = useState<string | null>(null);
  const [showAll, setShowAll] = useState(false);

  const highConfCalls = calls.filter(c => HIGH_CONFIDENCE_TIERS.has(c.tier));
  const filteredCalls = showAll ? calls : highConfCalls;
  const hiddenCount = calls.length - highConfCalls.length;

  const sorted = [...filteredCalls].sort((a, b) => {
    let cmp = 0;
    switch (sortKey) {
      case 'tier':
        cmp = (TIER_ORDER[a.tier] ?? 9) - (TIER_ORDER[b.tier] ?? 9);
        break;
      case 'score':
        cmp = b.score - a.score;
        break;
      case 'chrom_a':
        cmp = a.chrom_a.localeCompare(b.chrom_a) || a.pos_a - b.pos_a;
        break;
      case 'support':
        cmp = (b.support.total ?? (b.support.discordant + b.support.split + b.support.clipped)) -
              (a.support.total ?? (a.support.discordant + a.support.split + a.support.clipped));
        break;
      case 'pvalue':
        cmp = worstPval(a) - worstPval(b);
        break;
    }
    return sortAsc ? cmp : -cmp;
  });

  const handleSort = (key: SortKey) => {
    if (sortKey === key) setSortAsc(!sortAsc);
    else { setSortKey(key); setSortAsc(true); }
  };

  const tierCounts = calls.reduce((acc, c) => {
    acc[c.tier] = (acc[c.tier] || 0) + 1;
    return acc;
  }, {} as Record<string, number>);

  return (
    <div className="validated-callset">
      <div className="validated-header">
        <h2>Validated Translocation Calls</h2>
        <div className="validated-summary">
          {Object.entries(tierCounts)
            .sort(([a], [b]) => (TIER_ORDER[a] ?? 9) - (TIER_ORDER[b] ?? 9))
            .map(([tier, count]) => (
              <span key={tier} className="tier-badge" style={{ color: TIER_COLORS[tier] }}>
                {count} {TIER_LABELS[tier] || tier}
              </span>
            ))}
        </div>
        {jobId && (
          <div className="validated-export">
            <a href={`${BASE}/api/jobs/${jobId}/download/calls.vcf`} className="export-btn">VCF</a>
            <a href={`${BASE}/api/jobs/${jobId}/download/calls.bedpe`} className="export-btn">BEDPE</a>
            <a href={`${BASE}/api/jobs/${jobId}/download/calls.json`} className="export-btn">JSON</a>
          </div>
        )}
      </div>

      {artifactWarning && (
        <div className="artifact-warning" style={{
          background: '#44220a', border: '1px solid #f59e0b', borderRadius: 6,
          padding: '8px 14px', marginBottom: 12, color: '#fbbf24', fontSize: '0.85rem'
        }}>
          Warning: High number of calls detected. Results may be artifact-dominated.
          Review filter flags and evidence support carefully.
        </div>
      )}

      {!showAll && hiddenCount > 0 && (
        <div style={{ marginBottom: 10, fontSize: '0.85rem', color: '#9ca3af' }}>
          Showing {highConfCalls.length} high-confidence call{highConfCalls.length !== 1 ? 's' : ''}.{' '}
          <button
            onClick={() => setShowAll(true)}
            style={{ background: 'none', border: 'none', color: '#60a5fa', cursor: 'pointer', textDecoration: 'underline', padding: 0, fontSize: 'inherit' }}
          >
            Show {hiddenCount} candidate{hiddenCount !== 1 ? 's' : ''}
          </button>
        </div>
      )}
      {showAll && hiddenCount > 0 && (
        <div style={{ marginBottom: 10, fontSize: '0.85rem', color: '#9ca3af' }}>
          Showing all {calls.length} calls.{' '}
          <button
            onClick={() => setShowAll(false)}
            style={{ background: 'none', border: 'none', color: '#60a5fa', cursor: 'pointer', textDecoration: 'underline', padding: 0, fontSize: 'inherit' }}
          >
            Show high-confidence only
          </button>
        </div>
      )}

      <table className="validated-table">
        <thead>
          <tr>
            <th onClick={() => handleSort('tier')} className="sortable">
              Tier {sortKey === 'tier' ? (sortAsc ? '▲' : '▼') : ''}
            </th>
            <th onClick={() => handleSort('chrom_a')} className="sortable">
              Breakpoint A {sortKey === 'chrom_a' ? (sortAsc ? '▲' : '▼') : ''}
            </th>
            <th>Breakpoint B</th>
            <th>Orient</th>
            <th onClick={() => handleSort('support')} className="sortable">
              Support {sortKey === 'support' ? (sortAsc ? '▲' : '▼') : ''}
            </th>
            <th onClick={() => handleSort('score')} className="sortable">
              Score {sortKey === 'score' ? (sortAsc ? '▲' : '▼') : ''}
            </th>
            <th onClick={() => handleSort('pvalue')} className="sortable" title="Worst (max) of per-side local negative binomial p-values">
              NB P-value {sortKey === 'pvalue' ? (sortAsc ? '▲' : '▼') : ''}
            </th>
            <th>External</th>
          </tr>
        </thead>
        <tbody>
          <AnimatePresence>
            {sorted.map((call, i) => (
              <motion.tr
                key={call.event_id}
                initial={{ opacity: 0, y: 10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: i * 0.05, duration: 0.3 }}
                className={`call-row tier-${call.tier} ${selected === call.event_id ? 'selected' : ''}`}
                onClick={() => setSelected(selected === call.event_id ? null : call.event_id)}
              >
                <td>
                  <span className="tier-dot" style={{ backgroundColor: TIER_COLORS[call.tier] }} />
                  {TIER_LABELS[call.tier] || call.tier}
                </td>
                <td className="mono">
                  {call.chrom_a}:{formatPos(call.pos_a)}
                </td>
                <td className="mono">
                  {call.chrom_b}:{formatPos(call.pos_b)}
                </td>
                <td className="mono">{call.orientation}</td>
                <td>
                  <span title="Discordant pairs">{call.support.discordant}d</span>{' '}
                  <span title="Split reads" style={{ color: '#10b981' }}>{call.support.split}s</span>{' '}
                  <span title="Clipped reads">{call.support.clipped}c</span>
                </td>
                <td className="mono">{call.score.toFixed(1)}</td>
                <td className="mono">{formatPval(worstPval(call))}</td>
                <td>
                  {(call.external_callers?.length > 0) ? (
                    <span className="external-badge" title={call.external_callers.join(', ')}>
                      {call.external_callers.length}
                    </span>
                  ) : (
                    <span style={{ color: '#4b5563' }}>--</span>
                  )}
                </td>
              </motion.tr>
            ))}
          </AnimatePresence>
        </tbody>
      </table>

      {sorted.length === 0 && (
        <div className="validated-empty">
          No translocation calls passed validation filters.
        </div>
      )}

      {/* Detail panel for selected call */}
      <AnimatePresence>
        {selected && (() => {
          const call = calls.find(c => c.event_id === selected);
          if (!call) return null;
          return (
            <motion.div
              key="detail"
              initial={{ opacity: 0, height: 0 }}
              animate={{ opacity: 1, height: 'auto' }}
              exit={{ opacity: 0, height: 0 }}
              className="call-detail"
            >
              <h3>{call.event_id}</h3>
              <div className="detail-grid">
                {/* Breakpoints */}
                <div className="detail-section">
                  <div className="detail-section-title">Breakpoints</div>
                  <div>
                    <strong>Side A:</strong> {call.chrom_a}:{formatPos(call.pos_a)}
                    <span className="ci"> CI: [{call.ci_a[0]}, {call.ci_a[1]}]</span>
                    <span className="mapq-tag"> MAPQ: {call.median_mapq_a?.toFixed(0) ?? '--'}</span>
                  </div>
                  <div>
                    <strong>Side B:</strong> {call.chrom_b}:{formatPos(call.pos_b)}
                    <span className="ci"> CI: [{call.ci_b[0]}, {call.ci_b[1]}]</span>
                    <span className="mapq-tag"> MAPQ: {call.median_mapq_b?.toFixed(0) ?? '--'}</span>
                  </div>
                  <div><strong>Orientation:</strong> {call.orientation}</div>
                  <div><strong>Reciprocal support:</strong> {call.reciprocal_support}</div>
                  {call.evidence_label && <div><strong>Evidence:</strong> {call.evidence_label}</div>}
                </div>

                {/* Statistical evidence */}
                <div className="detail-section">
                  <div className="detail-section-title">Statistical Evidence</div>
                  <div>
                    <strong>NB p-value A:</strong> {formatPval(call.local_nb_pvalue_a ?? 1)}
                    <span className="rate-tag"> (rate: {call.local_rate_a?.toFixed(2) ?? '--'})</span>
                  </div>
                  <div>
                    <strong>NB p-value B:</strong> {formatPval(call.local_nb_pvalue_b ?? 1)}
                    <span className="rate-tag"> (rate: {call.local_rate_b?.toFixed(2) ?? '--'})</span>
                  </div>
                  <div>
                    <strong>Coverage ratio A:</strong> {call.local_coverage_ratio_a?.toFixed(2) ?? '--'}
                    {' / '}
                    <strong>B:</strong> {call.local_coverage_ratio_b?.toFixed(2) ?? '--'}
                  </div>
                  <div>
                    <strong>Chrom pair enrichment:</strong> {call.chrom_pair_enrichment?.toFixed(2) ?? '--'}
                  </div>
                  <div>
                    <strong>Duplicate fraction:</strong> {((call.duplicate_fraction ?? 0) * 100).toFixed(1)}%
                  </div>
                </div>

                {/* Uniqueness */}
                <div className="detail-section">
                  <div className="detail-section-title">Uniqueness</div>
                  <div>
                    <strong>Unique starts:</strong> A={call.unique_starts_a ?? '--'}, B={call.unique_starts_b ?? '--'}
                  </div>
                  <div>
                    <strong>Flanks remap uniquely:</strong>{' '}
                    <span style={{ color: call.both_flanks_remap_uniquely ? '#10b981' : '#ef4444' }}>
                      {call.both_flanks_remap_uniquely ? 'Yes' : 'No'}
                    </span>
                  </div>
                  {call.promiscuous_hotspot && (
                    <div style={{ color: '#f59e0b' }}><strong>Promiscuous hotspot</strong></div>
                  )}
                </div>

                {/* Masks */}
                <div className="detail-section">
                  <div className="detail-section-title">Mask Overlaps</div>
                  <div>
                    <strong>Side A:</strong>{' '}
                    {call.mask_overlaps_a?.length > 0 ? call.mask_overlaps_a.join(', ') : 'none'}
                    {call.segdup_pct_a > 0 && <span className="segdup-tag"> (segdup {call.segdup_pct_a.toFixed(0)}%)</span>}
                  </div>
                  <div>
                    <strong>Side B:</strong>{' '}
                    {call.mask_overlaps_b?.length > 0 ? call.mask_overlaps_b.join(', ') : 'none'}
                    {call.segdup_pct_b > 0 && <span className="segdup-tag"> (segdup {call.segdup_pct_b.toFixed(0)}%)</span>}
                  </div>
                </div>

                {/* External callers */}
                {call.external_callers?.length > 0 && (
                  <div className="detail-section">
                    <div className="detail-section-title">External Caller Agreement</div>
                    <div>{call.external_callers.join(', ')}</div>
                  </div>
                )}

                {/* Assembly / sequence features */}
                {(call.microhomology || call.inserted_seq || call.assembly_resolved) && (
                  <div className="detail-section">
                    <div className="detail-section-title">Sequence Features</div>
                    {call.assembly_resolved && (
                      <div><strong>Assembly resolved</strong> (offset: {call.assembly_bp_offset}bp)</div>
                    )}
                    {call.microhomology && <div><strong>Microhomology:</strong> {call.microhomology}</div>}
                    {call.inserted_seq && <div><strong>Inserted sequence:</strong> {call.inserted_seq}</div>}
                  </div>
                )}

                {/* Filter / reject info */}
                {(call.filter_flags?.length > 0 || call.reject_reasons?.length > 0) && (
                  <div className="detail-section">
                    <div className="detail-section-title">Filters</div>
                    {call.filter_flags?.length > 0 && (
                      <div><strong>Filter flags:</strong> {call.filter_flags.join(', ')}</div>
                    )}
                    {call.reject_reasons?.length > 0 && (
                      <div style={{ color: '#ef4444' }}>
                        <strong>Reject reasons:</strong> {call.reject_reasons.join(', ')}
                      </div>
                    )}
                  </div>
                )}

                {/* Score breakdown */}
                {call.score_components && Object.keys(call.score_components).length > 0 && (
                  <div className="detail-section">
                    <div className="detail-section-title">Score Components</div>
                    <div className="score-breakdown">
                      {Object.entries(call.score_components).map(([key, val]) => (
                        <span key={key} className="score-chip">
                          {key}: {typeof val === 'number' ? val.toFixed(1) : val}
                        </span>
                      ))}
                    </div>
                  </div>
                )}
              </div>
            </motion.div>
          );
        })()}
      </AnimatePresence>
    </div>
  );
}
