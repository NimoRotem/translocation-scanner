import { useEffect, useRef } from 'react';
import { useScanStore } from '../stores/scanStore';
import type { SSEEvent } from '../types/events';

const BASE = '/translocation-scanner';

export function useSSE(jobId: string | null) {
  const handleEvent = useScanStore(s => s.handleEvent);
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
    ];

    eventTypes.forEach(type => {
      es.addEventListener(type, (e: MessageEvent) => {
        try {
          const data = JSON.parse(e.data);
          if (type === 'stream.end') {
            es.close();
            return;
          }
          if (type === 'job.state') {
            // Initial state — could restore from it
            return;
          }
          handleEvent({ ...data, type } as SSEEvent);
        } catch (err) {
          console.error('SSE parse error:', type, err);
        }
      });
    });

    es.onerror = () => {
      console.warn('SSE connection error, will retry...');
    };

    return () => {
      es.close();
      esRef.current = null;
    };
  }, [jobId, handleEvent]);
}
