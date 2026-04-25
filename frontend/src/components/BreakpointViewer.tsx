import type { ValidatedCall } from '../types/events';

interface BreakpointViewerProps {
  call: ValidatedCall;
  jobId: string;
}

const BASE = '/translocation-scanner';

export default function BreakpointViewer({ call, jobId }: BreakpointViewerProps) {
  return (
    <div className="breakpoint-viewer" style={{
      background: 'var(--bg-card)',
      border: '1px solid var(--border)',
      borderRadius: 8,
      padding: 20,
    }}>
      <h3 style={{ color: 'var(--accent-green)', fontFamily: 'var(--font-mono)', marginBottom: 12 }}>
        Breakpoint Pileup — {call.event_id}
      </h3>
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 16 }}>
        <div style={{
          background: 'var(--bg)',
          border: '1px solid var(--border)',
          borderRadius: 4,
          padding: 16,
          minHeight: 120,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontFamily: 'var(--font-mono)',
          color: 'var(--text-muted)',
          fontSize: 12,
        }}>
          {call.chrom_a}:{call.pos_a.toLocaleString()}
          <br />
          {call.support.discordant}d {call.support.split}s {call.support.clipped}c
          <br />
          <span style={{ fontSize: 10, marginTop: 8, display: 'block' }}>
            Pileup viewer coming in M2
          </span>
        </div>
        <div style={{
          background: 'var(--bg)',
          border: '1px solid var(--border)',
          borderRadius: 4,
          padding: 16,
          minHeight: 120,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          fontFamily: 'var(--font-mono)',
          color: 'var(--text-muted)',
          fontSize: 12,
        }}>
          {call.chrom_b}:{call.pos_b.toLocaleString()}
          <br />
          Reciprocal: {call.reciprocal_support}
          <br />
          <span style={{ fontSize: 10, marginTop: 8, display: 'block' }}>
            Pileup viewer coming in M2
          </span>
        </div>
      </div>
    </div>
  );
}
