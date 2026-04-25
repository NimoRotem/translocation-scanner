import { useRef, useEffect, useState } from 'react';
import type { WaterfallEntry } from '../types/events';

interface WaterfallProps {
  entries: WaterfallEntry[];
}

const TYPE_COLORS: Record<string, string> = {
  split: '#10b981',
  clip_pileup: '#3b82f6',
  new_pair_cluster: '#f59e0b',
  discordant: '#6b7280',
};

const TYPE_ICONS: Record<string, string> = {
  split: '★',
  clip_pileup: '◆',
  new_pair_cluster: '●',
  discordant: '+',
};

function formatPos(pos: number): string {
  if (pos >= 1_000_000) return (pos / 1_000_000).toFixed(1) + 'Mb';
  if (pos >= 1_000) return (pos / 1_000).toFixed(1) + 'kb';
  return String(pos);
}

function formatTime(ts: number, base: number): string {
  const s = (ts - base) / 1000;
  const m = Math.floor(s / 60);
  const sec = (s % 60).toFixed(1);
  return `${String(m).padStart(2, '0')}:${sec.padStart(4, '0')}`;
}

export default function Waterfall({ entries }: WaterfallProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const [autoScroll, setAutoScroll] = useState(true);
  const [filter, setFilter] = useState<string | null>(null);
  const baseTime = entries.length > 0 ? entries[0].timestamp : Date.now();

  useEffect(() => {
    if (autoScroll && containerRef.current) {
      containerRef.current.scrollTop = containerRef.current.scrollHeight;
    }
  }, [entries.length, autoScroll]);

  const handleScroll = () => {
    const el = containerRef.current;
    if (!el) return;
    const atBottom = el.scrollHeight - el.scrollTop - el.clientHeight < 40;
    setAutoScroll(atBottom);
  };

  const filtered = filter
    ? entries.filter(e => e.type === filter)
    : entries;

  // Show max 200 visible entries
  const visible = filtered.slice(-200);

  return (
    <div className="waterfall">
      <div className="waterfall-header">
        <span className="waterfall-title">Evidence Waterfall</span>
        <div className="waterfall-filters">
          {['split', 'clip_pileup', 'new_pair_cluster', 'discordant'].map(t => (
            <button
              key={t}
              className={`waterfall-filter-btn ${filter === t ? 'active' : ''}`}
              style={{ color: TYPE_COLORS[t] }}
              onClick={() => setFilter(filter === t ? null : t)}
            >
              {TYPE_ICONS[t]} {t.replace('_', ' ')}
            </button>
          ))}
        </div>
      </div>
      <div
        ref={containerRef}
        className="waterfall-scroll"
        onScroll={handleScroll}
      >
        {visible.map(entry => (
          <div key={entry.id} className="waterfall-row">
            <span className="wf-time">{formatTime(entry.timestamp, baseTime)}</span>
            <span className="wf-locus">
              {entry.chrom_a}:{formatPos(entry.pos_a)} ↔ {entry.chrom_b}:{formatPos(entry.pos_b)}
            </span>
            <span className="wf-type" style={{ color: TYPE_COLORS[entry.type] || '#9ca3af' }}>
              {TYPE_ICONS[entry.type] || '·'} {entry.type.replace('_', ' ')}
            </span>
            {entry.detail && (
              <span className="wf-detail">{entry.detail}</span>
            )}
          </div>
        ))}
        {visible.length === 0 && (
          <div className="waterfall-empty">Waiting for evidence events...</div>
        )}
      </div>
      {!autoScroll && (
        <button
          className="waterfall-jump"
          onClick={() => {
            setAutoScroll(true);
            if (containerRef.current) {
              containerRef.current.scrollTop = containerRef.current.scrollHeight;
            }
          }}
        >
          ↓ Jump to live
        </button>
      )}
    </div>
  );
}
