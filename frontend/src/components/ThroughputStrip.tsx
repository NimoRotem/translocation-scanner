import Sparkline from './Sparkline';

interface ThroughputSnapshot {
  t: number;
  reads_per_sec: number;
  bytes_per_sec: number;
  discordant_per_sec: number;
  split_per_sec: number;
}

interface ThroughputStripProps {
  current: ThroughputSnapshot;
  history: ThroughputSnapshot[];
}

function formatNum(n: number): string {
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + 'M';
  if (n >= 1_000) return (n / 1_000).toFixed(1) + 'k';
  return n.toFixed(0);
}

function formatBytes(n: number): string {
  if (n >= 1_073_741_824) return (n / 1_073_741_824).toFixed(1) + ' GB/s';
  if (n >= 1_048_576) return (n / 1_048_576).toFixed(0) + ' MB/s';
  if (n >= 1024) return (n / 1024).toFixed(0) + ' KB/s';
  return n.toFixed(0) + ' B/s';
}

const rows: Array<{
  key: keyof ThroughputSnapshot;
  label: string;
  format: (n: number) => string;
  color: string;
}> = [
  { key: 'reads_per_sec', label: 'reads/sec', format: formatNum, color: '#f59e0b' },
  { key: 'bytes_per_sec', label: 'bytes/sec', format: formatBytes, color: '#3b82f6' },
  { key: 'discordant_per_sec', label: 'discordant/sec', format: formatNum, color: '#ef4444' },
  { key: 'split_per_sec', label: 'splits/sec', format: formatNum, color: '#10b981' },
];

export default function ThroughputStrip({ current, history }: ThroughputStripProps) {
  return (
    <div className="throughput-strip">
      {rows.map(row => {
        const val = current[row.key] as number;
        const sparkData = history.map(h => h[row.key] as number);

        // Detect stall: if reads/sec drops below 20% of median for >5s
        let warning = false;
        if (row.key === 'reads_per_sec' && sparkData.length > 5) {
          const sorted = [...sparkData].sort((a, b) => a - b);
          const median = sorted[Math.floor(sorted.length / 2)];
          const recent = sparkData.slice(-5);
          if (median > 0 && recent.every(v => v < median * 0.2)) {
            warning = true;
          }
        }

        return (
          <div key={row.key} className="throughput-row">
            <span className="throughput-label">{row.label}</span>
            <span className="throughput-value" style={{ color: row.color }}>
              {row.format(val)}
              {warning && <span className="throughput-warn" title="I/O stall detected"> ⚠</span>}
            </span>
            <Sparkline data={sparkData} width={100} height={16} color={row.color} />
          </div>
        );
      })}
    </div>
  );
}
