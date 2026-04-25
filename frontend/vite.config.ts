import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  base: '/translocation-scanner/',
  server: {
    proxy: {
      '/translocation-scanner/api': {
        target: 'http://localhost:8750',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/translocation-scanner/, ''),
      },
    },
  },
  build: {
    outDir: 'dist',
    sourcemap: true,
  },
})
