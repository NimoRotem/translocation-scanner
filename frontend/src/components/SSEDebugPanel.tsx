import { useEffect, useState } from 'react';
import { useScanStore } from '../stores/scanStore';

export function SSEDebugPanel() {
  const [visible, setVisible] = useState(false);
  const debug = useScanStore(s => s.sseDebug);
  const stage = useScanStore(s => s.stage);
  const jobId = useScanStore(s => s.jobId);

  // Toggle with Ctrl+Shift+D
  useEffect(() => {
    const handler = (e: KeyboardEvent) => {
      if (e.ctrlKey && e.shiftKey && e.key === 'D') {
        e.preventDefault();
        setVisible(v => !v);
      }
    };
    window.addEventListener('keydown', handler);
    return () => window.removeEventListener('keydown', handler);
  }, []);

  if (!visible) return null;

  const eventsPerSec = debug.eventCount > 0 && debug.lastEventTime > 0
    ? Math.round(debug.eventCount / ((Date.now() - (debug.lastEventTime - debug.eventCount * 200)) / 1000) * 10) / 10
    : 0;

  const timeSince = debug.lastEventTime
    ? `${((Date.now() - debug.lastEventTime) / 1000).toFixed(1)}s ago`
    : 'never';

  return (
    <div style={{
      position: 'fixed',
      bottom: 16,
      right: 16,
      width: 320,
      background: 'rgba(0, 0, 0, 0.9)',
      border: '1px solid #444',
      borderRadius: 8,
      padding: 12,
      fontFamily: 'monospace',
      fontSize: 11,
      color: '#ccc',
      zIndex: 9999,
    }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 8 }}>
        <strong style={{ color: '#fff' }}>SSE Debug</strong>
        <button
          onClick={() => setVisible(false)}
          style={{ background: 'none', border: 'none', color: '#888', cursor: 'pointer', fontSize: 14 }}
        >
          ×
        </button>
      </div>

      <div style={{ marginBottom: 4 }}>
        <span style={{ color: debug.connected ? '#4ade80' : '#ef4444' }}>●</span>
        {' '}{debug.connected ? 'Connected' : 'Disconnected'}
      </div>

      <div style={{ marginBottom: 4 }}>
        Job: <span style={{ color: '#93c5fd' }}>{jobId || 'none'}</span>
        {' | '}Stage: <span style={{ color: '#fbbf24' }}>{stage || 'none'}</span>
      </div>

      <div style={{ marginBottom: 4 }}>
        Last event: <span style={{ color: '#c4b5fd' }}>{debug.lastEventType || 'none'}</span>
        {' '}({timeSince})
      </div>

      <div style={{ marginBottom: 4 }}>
        Total events: <strong>{debug.eventCount}</strong>
        {' | ~'}{eventsPerSec}/sec
      </div>

      {debug.errors.length > 0 && (
        <div style={{ marginTop: 8, borderTop: '1px solid #333', paddingTop: 4 }}>
          <div style={{ color: '#ef4444', marginBottom: 4 }}>Errors ({debug.errors.length}):</div>
          {debug.errors.slice(-5).map((err, i) => (
            <div key={i} style={{ color: '#fca5a5', fontSize: 10 }}>
              {new Date(err.time).toLocaleTimeString()}: {err.message}
            </div>
          ))}
        </div>
      )}

      <div style={{ marginTop: 8, color: '#666', fontSize: 10 }}>
        Ctrl+Shift+D to toggle
      </div>
    </div>
  );
}
