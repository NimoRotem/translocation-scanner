import { test, expect } from '@playwright/test';

const BASE_URL = 'https://23andclaude.com/translocation-scanner';
const API_BASE = `${BASE_URL}/api`;

test.describe('Report Rendering', () => {
  test('scanner page loads without errors', async ({ page }) => {
    const consoleErrors: string[] = [];
    page.on('console', msg => {
      if (msg.type() === 'error') {
        consoleErrors.push(msg.text());
      }
    });

    const response = await page.goto(BASE_URL);
    expect(response?.status()).toBe(200);
    await page.waitForLoadState('networkidle');

    // Page should not be blank
    const bodyText = await page.textContent('body');
    expect(bodyText?.length).toBeGreaterThan(0);

    // Filter out non-fatal errors
    const realErrors = consoleErrors.filter(e =>
      !e.includes('favicon') && !e.includes('SSE')
    );
    expect(realErrors).toHaveLength(0);
  });

  test('health endpoint responds', async ({ page }) => {
    const res = await page.request.get(`${API_BASE}/health`);
    expect(res.ok()).toBeTruthy();
    const data = await res.json();
    expect(data.status).toBe('ok');
  });

  test('server-files endpoint responds', async ({ page }) => {
    const res = await page.request.get(`${API_BASE}/server-files`);
    expect(res.ok()).toBeTruthy();
    const data = await res.json();
    expect(data.files).toBeDefined();
    expect(data.references).toBeDefined();
  });

  test('report endpoint for completed job has valid structure', async ({ page }) => {
    // List jobs and find a completed one
    const jobsRes = await page.request.get(`${API_BASE}/jobs`);
    if (!jobsRes.ok()) {
      test.skip();
      return;
    }
    const { jobs } = await jobsRes.json();
    const completed = jobs.find((j: any) => j.status === 'completed');
    if (!completed) {
      test.skip();
      return;
    }

    const reportRes = await page.request.get(`${API_BASE}/jobs/${completed.job_id}/report`);
    expect(reportRes.ok()).toBeTruthy();
    const report = await reportRes.json();

    // Required sections
    expect(report.sample).toBeTruthy();
    expect(report.quality).toBeTruthy();
    expect(report.evidence).toBeTruthy();
    expect(report.pipeline).toBeTruthy();
    expect(report.results).toBeTruthy();
    expect(report.interpretation).toBeTruthy();
    expect(report.interpretation.summary).toBeTruthy();
    expect(report.interpretation.detail).toBeTruthy();

    // Tier counts should sum correctly
    const byTier = report.results.by_tier;
    const tierSum = Object.values(byTier).reduce((a: number, b: any) => a + b, 0);
    expect(tierSum).toBe(report.results.total_calls);
  });
});
