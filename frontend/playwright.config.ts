import { defineConfig } from '@playwright/test';

export default defineConfig({
  testDir: './tests',
  timeout: 300_000,  // 5 min per test (nano BAM takes ~17s + startup)
  retries: 0,
  use: {
    baseURL: 'https://23andclaude.com/translocation-scanner',
    ignoreHTTPSErrors: true,
  },
  projects: [
    {
      name: 'chromium',
      use: { browserName: 'chromium', headless: true },
    },
  ],
});
