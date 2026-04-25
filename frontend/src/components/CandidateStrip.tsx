interface CandidateStripProps {
  topBins: Array<{
    chrom_a: string;
    chrom_b: string;
    count: number;
    pos_a: number;
    pos_b: number;
  }>;
  mode: 'streaming' | 'validated' | 'idle';
}

function formatPos(pos: number): string {
  if (pos === 0) return '—';
  if (pos >= 1_000_000) return '~' + (pos / 1_000_000).toFixed(1) + 'M';
  if (pos >= 1_000) return '~' + (pos / 1_000).toFixed(0) + 'k';
  return String(pos);
}

export default function CandidateStrip({ topBins, mode }: CandidateStripProps) {
  if (topBins.length === 0) return null;

  const sorted = [...topBins].sort((a, b) => b.count - a.count);
  const maxCount = sorted[0]?.count || 1;

  return (
    <div className="candidate-strip">
      <div className="candidate-strip-label">
        {mode === 'streaming' ? 'PROVISIONAL · Top bins by raw pair count' : 'Top evidence bins'}
      </div>
      <div className="candidate-strip-scroll">
        {sorted.slice(0, 12).map((bin, i) => {
          const intensity = Math.min(1, bin.count / maxCount);
          return (
            <div
              key={`${bin.chrom_a}-${bin.chrom_b}-${i}`}
              className="candidate-card"
              style={{
                borderColor: `rgba(245, 158, 11, ${0.3 + intensity * 0.7})`,
                backgroundColor: `rgba(245, 158, 11, ${0.05 + intensity * 0.1})`,
              }}
            >
              <div className="candidate-pair">
                {bin.chrom_a.replace('chr', '')} ↔ {bin.chrom_b.replace('chr', '')}
              </div>
              <div className="candidate-count">{bin.count} pairs</div>
              <div className="candidate-pos">
                {formatPos(bin.pos_a)} / {formatPos(bin.pos_b)}
              </div>
              {mode === 'streaming' && (
                <div className="candidate-provisional">PROVISIONAL</div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
