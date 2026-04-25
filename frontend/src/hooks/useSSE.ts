import { useEffect, useRef } from 'react';
import { useScanStore } from '../stores/scanStore';
import type { SSEEvent } from '../types/events';

const BASE = '/translocation-scanner';

export function useSSE(jobId: string | null) {
  const handleEvent = useScanStore(s => s.handleEvent);
  const setSSEDebug = useScanStore(s => s.setSSEDebug);
  const esRef = useRef<EventSource | null>(null);

  useEffect(() => {
    if (!jobId) return;

    const url = `${BASE}/api/jobs/${jobId}/stream`;
    const es = new EventSource(url);
    esRef.current = es;

    const eventTypes = [
      'scan.started', 'scan.stage_changed', 'scan.progress', 'scan.throughput',
      'chrom.progress', 'chrom.bin_update', 'pair.density', 'evidence.highlight',
      'provisional.top_bins', 'scan.completed', 'validation.started',
      'validation.call_emitted', 'validation.completed', 'scan.cancelled', 'error',
      'job.state', 'stream.end',
      // New events from rebuild
      'chrom_pair.matrix', 'clustering.scale',
    ];

    eventTypes.forEach(type => {
      es.addEventListener(type, (e: MessageEvent) => {
        try {
          const data = JSON.parse(e.data);
          if (type === 'stream.end') {
            es.close();
            setSSEDebug({ connected: false });
            return;
          }
          if (type === 'job.state') {
            return;
          }
          handleEvent({ ...data, type } as SSEEvent);
          setSSEDebug({
            lastEventType: type,
            lastEventTime: Date.now(),
            eventCount: (useScanStore.getState().sseDebug?.eventCount ?? 0) + 1,
          });
        } catch (err) {
          console.error('SSE parse error:', type, err);
          setSSEDebug({
            errors: [
              ...((useScanStore.getState().sseDebug?.errors ?? []).slice(-9)),
              { time: Date.now(), message: `Parse error: ${type}` },
            ],
          });
        }
      });
    });

    es.onopen = () => {
      setSSEDebug({ connected: true, jobId });
    };

    es.onerror = () => {
      console.warn('SSE connection error, will retry...');
      setSSEDebug({
        connected: false,
        errors: [
          ...((useScanStore.getState().sseDebug?.errors ?? []).slice(-9)),
          { time: Date.now(), message: 'Connection error' },
        ],
      });
    };

    return () => {
      es.close();
      esRef.current = null;
      setSSEDebug({ connected: false });
    };
  }, [jobId, handleEvent, setSSEDebug]);
}
