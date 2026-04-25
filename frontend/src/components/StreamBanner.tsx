import { motion, AnimatePresence } from 'framer-motion';
import type { ScanMode } from '../types/events';

interface StreamBannerProps {
  mode: ScanMode;
  stage: string;
  elapsed: number;
  fileName: string;
  onCancel?: () => void;
}

function formatTime(ms: number): string {
  const s = Math.floor(ms / 1000);
  const h = Math.floor(s / 3600);
  const m = Math.floor((s % 3600) / 60);
  const sec = s % 60;
  return `${String(h).padStart(2, '0')}:${String(m).padStart(2, '0')}:${String(sec).padStart(2, '0')}`;
}

export default function StreamBanner({ mode, stage, elapsed, fileName, onCancel }: StreamBannerProps) {
  if (mode === 'idle') return null;

  const isStreaming = mode === 'streaming';
  const bgColor = isStreaming ? '#78350f' : '#064e3b';
  const accentColor = isStreaming ? '#f59e0b' : '#10b981';
  const label = isStreaming ? 'STREAMING' : 'VALIDATED';
  const subtitle = isStreaming
    ? 'Live extraction in progress — provisional visuals below reflect raw evidence, not final calls'
    : 'Scan complete — calibrated, artifact-filtered callset below';

  const shortName = fileName.split('/').pop() || fileName;

  return (
    <motion.div
      className="stream-banner"
      animate={{ backgroundColor: bgColor }}
      transition={{ duration: 1.5 }}
      style={{
        padding: '12px 20px',
        borderBottom: `2px solid ${accentColor}`,
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center',
        fontFamily: 'var(--font-mono)',
      }}
    >
      <div style={{ display: 'flex', alignItems: 'center', gap: 16 }}>
        <span style={{ fontSize: 13, color: 'var(--text-muted)' }}>
          Sample: <strong style={{ color: 'var(--text)' }}>{shortName}</strong>
        </span>
        <span style={{ fontSize: 13, color: 'var(--text-muted)' }}>
          Runtime: <strong style={{ color: 'var(--text)' }}>{formatTime(elapsed)}</strong>
        </span>
        {stage && stage !== 'completed' && (
          <span style={{ fontSize: 12, color: 'var(--text-muted)', textTransform: 'capitalize' }}>
            Stage: {stage.replace(/_/g, ' ')}
          </span>
        )}
      </div>
      <div style={{ display: 'flex', alignItems: 'center', gap: 10 }}>
        {isStreaming && onCancel && (
          <button
            onClick={onCancel}
            style={{
              background: 'rgba(239, 68, 68, 0.15)',
              border: '1px solid rgba(239, 68, 68, 0.4)',
              color: '#ef4444',
              padding: '4px 12px',
              borderRadius: 4,
              fontSize: 12,
              fontFamily: 'var(--font-mono)',
              fontWeight: 600,
              cursor: 'pointer',
              letterSpacing: '0.05em',
            }}
          >
            Stop Scan
          </button>
        )}
        <motion.span
          animate={{ color: accentColor }}
          transition={{ duration: 1.5 }}
          style={{ fontSize: 14, fontWeight: 700, letterSpacing: '0.1em' }}
        >
          {label}
        </motion.span>
        <motion.div
          animate={{
            backgroundColor: accentColor,
            boxShadow: isStreaming
              ? ['0 0 4px #f59e0b', '0 0 12px #f59e0b', '0 0 4px #f59e0b']
              : '0 0 8px #10b981',
          }}
          transition={{ duration: isStreaming ? 1 : 0, repeat: isStreaming ? Infinity : 0 }}
          style={{ width: 10, height: 10, borderRadius: '50%' }}
        />
      </div>
      <div style={{
        position: 'absolute',
        bottom: -1,
        left: 0,
        right: 0,
        fontSize: 11,
        color: accentColor,
        opacity: 0.7,
        textAlign: 'center',
        padding: '2px 0',
        pointerEvents: 'none',
      }}>
        {subtitle}
      </div>
    </motion.div>
  );
}
