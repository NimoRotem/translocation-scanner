import { useState, useEffect } from 'react';
import { useScanStore } from './stores/scanStore';
import { useSSE } from './hooks/useSSE';
import StreamBanner from './components/StreamBanner';
import ScanSetup from './components/ScanSetup';
import CircosRing from './components/CircosRing';
import ThroughputStrip from './components/ThroughputStrip';
import ChromGrid from './components/ChromGrid';
import DensityHeatmap from './components/DensityHeatmap';
import Waterfall from './components/Waterfall';
import CandidateStrip from './components/CandidateStrip';
import ValidatedCallset from './components/ValidatedCallset';
import ScanReport from './components/ScanReport';
import PreviousRuns from './components/PreviousRuns';

export default function App() {
  const {
    mode, jobId, stage, error, startedAt,
    throughput, throughputHistory,
    chromProgress, densityMatrix,
    waterfall, provisionalArcs, topBins,
    validatedCalls,
    startScan, reset,
  } = useScanStore();

  const [filePath, setFilePath] = useState('');

  useSSE(jobId);

  // Elapsed time counter
  const [elapsed, setElapsed] = useState(0);
  useEffect(() => {
    if (mode === 'idle' || !startedAt) return;
    const iv = setInterval(() => {
      setElapsed(Date.now() - startedAt);
    }, 1000);
    return () => clearInterval(iv);
  }, [mode, startedAt]);

  const handleStart = (newJobId: string, path: string) => {
    startScan(newJobId);
    setFilePath(path);
  };

  const handleNewScan = () => {
    reset();
    setFilePath('');
    setElapsed(0);
  };

  const handleResumeScan = (resumeJobId: string, path: string) => {
    startScan(resumeJobId);
    setFilePath(path);
  };

  const handleCancel = async () => {
    if (!jobId) return;
    try {
      await fetch(`/translocation-scanner/api/jobs/${jobId}/cancel`, { method: 'POST' });
    } catch (err) {
      console.error('Failed to cancel scan:', err);
    }
  };

  return (
    <div className={`scanner-app mode-${mode}`}>
      {mode !== 'idle' && (
        <StreamBanner
          mode={mode}
          stage={stage}
          elapsed={elapsed}
          fileName={filePath}
          onCancel={handleCancel}
        />
      )}

      {error && (
        <div className="error-bar">
          Error: {error}
          <button onClick={handleNewScan} className="error-dismiss">New Scan</button>
        </div>
      )}

      {mode === 'idle' && (
        <>
          <ScanSetup onStart={handleStart} />
          <div style={{ maxWidth: 700, margin: '0 auto', padding: '0 24px 40px' }}>
            <PreviousRuns
              onViewJob={() => {}}
              onResumeScan={handleResumeScan}
            />
          </div>
        </>
      )}

      {mode === 'streaming' && (
        <div className="streaming-dashboard">
          <div className="dash-top">
            <div className="dash-circos">
              <CircosRing
                chromProgress={chromProgress}
                provisionalArcs={provisionalArcs}
                validatedCalls={[]}
                mode={mode}
                size={460}
              />
            </div>
            <div className="dash-right-top">
              <ThroughputStrip current={throughput} history={throughputHistory} />
            </div>
          </div>
          <div className="dash-middle">
            <div className="dash-chrom">
              <ChromGrid chromProgress={chromProgress} />
            </div>
            <div className="dash-density">
              <DensityHeatmap densityMatrix={densityMatrix} size={340} />
            </div>
          </div>
          <div className="dash-bottom">
            <Waterfall entries={waterfall} />
          </div>
          <CandidateStrip topBins={topBins} mode="streaming" />
        </div>
      )}

      {mode === 'validated' && (
        <div className="validated-dashboard">
          <div className="validated-top">
            <div className="validated-circos">
              <CircosRing
                chromProgress={chromProgress}
                provisionalArcs={[]}
                validatedCalls={validatedCalls}
                mode={mode}
                size={400}
              />
            </div>
            <div className="validated-summary-side">
              <div className="validated-calls-count">
                {validatedCalls.length === 0 ? (
                  <span className="no-calls-badge">No translocations detected</span>
                ) : (
                  <>
                    <span className="calls-number">{validatedCalls.length}</span>
                    <span className="calls-label">
                      translocation{validatedCalls.length !== 1 ? 's' : ''} detected
                    </span>
                  </>
                )}
              </div>
              {validatedCalls.length > 0 && (
                <div className="validated-tier-summary">
                  {['confirmed', 'validated', 'likely', 'strong_candidate', 'candidate'].map(tier => {
                    const n = validatedCalls.filter(c => c.tier === tier).length;
                    if (n === 0) return null;
                    const labels: Record<string, string> = {
                      confirmed: 'confirmed', validated: 'validated', likely: 'likely',
                      strong_candidate: 'strong candidate', candidate: 'candidate',
                    };
                    return (
                      <div key={tier} className={`tier-badge tier-${tier}`}>
                        {n} {labels[tier] || tier}
                      </div>
                    );
                  })}
                </div>
              )}
            </div>
          </div>
          <ScanReport jobId={jobId} validatedCalls={validatedCalls} />
          <div className="validated-actions">
            <button onClick={handleNewScan} className="new-scan-btn">
              New Scan
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
