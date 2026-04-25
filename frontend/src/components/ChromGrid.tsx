import { CHROMS } from '../stores/scanStore';
import type { ChromData } from '../types/events';

interface ChromGridProps {
  chromProgress: Record<string, ChromData>;
}

function formatNum(n: number): string {
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + 'M';
  if (n >= 1_000) return (n / 1_000).toFixed(1) + 'k';
  return String(n);
}

export default function ChromGrid({ chromProgress }: ChromGridProps) {
  return (
    <div className="chrom-grid">
      <div className="chrom-grid-header">
        <span>Chr</span>
        <span>Progress</span>
        <span>Reads</span>
        <span>Disc</span>
        <span>Split</span>
      </div>
      {CHROMS.map(chrom => {
        const data = chromProgress[chrom];
        if (!data) return null;
        const pct = Math.min(100, data.pct);
        const statusClass = data.status === 'complete' ? 'chrom-complete'
          : data.status === 'scanning' ? 'chrom-scanning'
          : 'chrom-pending';

        // Anomalous density detection
        const anomalous = data.discordant > 0 &&
          data.reads > 0 &&
          (data.discordant / data.reads) > 0.05;

        const label = chrom.replace('chr', '');

        return (
          <div
            key={chrom}
            className={`chrom-row ${statusClass} ${anomalous ? 'chrom-anomalous' : ''}`}
          >
            <span className="chrom-label">{label}</span>
            <div className="chrom-bar-container">
              <div
                className="chrom-bar-fill"
                style={{
                  width: `${pct}%`,
                  backgroundColor: data.status === 'complete' ? '#10b981'
                    : data.status === 'scanning' ? '#f59e0b'
                    : '#374151',
                }}
              />
              <span className="chrom-bar-pct">{pct.toFixed(0)}%</span>
            </div>
            <span className="chrom-stat">{formatNum(data.reads)}</span>
            <span className="chrom-stat chrom-disc">{formatNum(data.discordant)}</span>
            <span className="chrom-stat chrom-split">{formatNum(data.split)}</span>
          </div>
        );
      })}
    </div>
  );
}
