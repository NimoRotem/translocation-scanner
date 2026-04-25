import { test, expect } from '@playwright/test';

const API_BASE = 'https://23andclaude.com/translocation-scanner/api';
const NANO_BAM = '/data/scan_archive/test_corpus/nano.bam';

test.describe('SSE Streaming', () => {
  test('full scan lifecycle via browser', async ({ page }) => {
    // Step 1: Start a scan via API
    const startRes = await page.request.post(`${API_BASE}/scan`, {
      data: {
        file_path: NANO_BAM,
        reference_build: 'GRCh38',
        settings: {
          skip_external_callers: true,
          skip_clip_realignment: true,
        },
      },
    });
    expect(startRes.ok()).toBeTruthy();
    const { job_id } = await startRes.json();
    expect(job_id).toBeTruthy();

    // Step 2: Navigate to scanner page (will auto-connect SSE)
    await page.goto('https://23andclaude.com/translocation-scanner/');
    await page.waitForLoadState('networkidle');

    // The scanner starts in idle mode; we need to trigger scan start from the UI
    // or navigate to a URL that auto-starts SSE. For now, we poll the API.

    // Step 3: Poll until job completes (max 120s)
    let jobStatus = 'running';
    let attempts = 0;
    while (jobStatus !== 'completed' && jobStatus !== 'failed' && attempts < 120) {
      const statusRes = await page.request.get(`${API_BASE}/jobs/${job_id}`);
      const statusData = await statusRes.json();
      jobStatus = statusData.status;
      if (jobStatus === 'completed' || jobStatus === 'failed') break;
      await page.waitForTimeout(1000);
      attempts++;
    }

    expect(jobStatus).toBe('completed');

    // Step 4: Verify report endpoint returns valid data
    const reportRes = await page.request.get(`${API_BASE}/jobs/${job_id}/report`);
    expect(reportRes.ok()).toBeTruthy();
    const report = await reportRes.json();

    // Verify report structure
    expect(report.sample).toBeTruthy();
    expect(report.quality).toBeTruthy();
    expect(report.evidence).toBeTruthy();
    expect(report.pipeline).toBeTruthy();
    expect(report.results).toBeTruthy();
    expect(report.interpretation).toBeTruthy();

    // Verify timings exist
    expect(report.pipeline.timings).toBeTruthy();
    expect(report.pipeline.timings.extraction).toBeGreaterThan(0);

    // Step 5: Check no JS console errors during page load
    const consoleErrors: string[] = [];
    page.on('console', msg => {
      if (msg.type() === 'error') {
        consoleErrors.push(msg.text());
      }
    });
    await page.reload();
    await page.waitForLoadState('networkidle');

    // Filter out expected warnings
    const realErrors = consoleErrors.filter(e =>
      !e.includes('SSE') && !e.includes('favicon')
    );
    expect(realErrors).toHaveLength(0);
  });

  test('SSE event stream delivers events', async ({ page }) => {
    // Start a scan
    const startRes = await page.request.post(`${API_BASE}/scan`, {
      data: {
        file_path: NANO_BAM,
        reference_build: 'GRCh38',
        settings: {
          skip_external_callers: true,
          skip_clip_realignment: true,
        },
      },
    });
    const { job_id } = await startRes.json();

    // Connect to SSE and collect events
    const events: string[] = [];
    const sseUrl = `${API_BASE}/jobs/${job_id}/stream`;

    // Use page.evaluate to test EventSource in-browser
    const collectedEvents = await page.evaluate(async (url) => {
      return new Promise<string[]>((resolve) => {
        const collected: string[] = [];
        const es = new EventSource(url);
        const types = [
          'scan.started', 'scan.stage_changed', 'scan.progress',
          'scan.completed', 'validation.completed', 'stream.end',
          'chrom_pair.matrix', 'clustering.scale',
        ];
        types.forEach(type => {
          es.addEventListener(type, () => {
            collected.push(type);
            if (type === 'stream.end' || type === 'validation.completed') {
              es.close();
              resolve(collected);
            }
          });
        });
        // Timeout after 120s
        setTimeout(() => {
          es.close();
          resolve(collected);
        }, 120_000);
      });
    }, sseUrl);

    expect(collectedEvents.length).toBeGreaterThan(0);
    expect(collectedEvents).toContain('scan.stage_changed');
  });
});
