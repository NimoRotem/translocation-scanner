import React, { useRef, useEffect, useCallback } from 'react';
import { CHROMS } from '../stores/scanStore';

interface DensityHeatmapProps {
  densityMatrix: number[][]; // 24x24 array from store
  size?: number; // default 350
}

/** Abbreviated chromosome labels: chr1-chr22 -> 1-22, chrX -> X, chrY -> Y */
const LABELS: string[] = CHROMS.map((c) =>
  c.replace('chr', ''),
);

/** Custom log-scaled color ramp based on count thresholds */
function countToColor(count: number): string {
  if (count <= 0) return '#0f0f0f';
  if (count <= 5) return '#1e3a5f';
  if (count <= 20) return '#2563eb';
  if (count <= 50) return '#f59e0b';
  if (count <= 100) return '#ef4444';
  return '#fef3c7';
}

/**
 * Interpolate between two hex colors by factor t (0..1).
 * Used for the pulse brightening effect.
 */
function lerpColor(a: string, b: string, t: number): string {
  const parseHex = (hex: string) => {
    const v = parseInt(hex.slice(1), 16);
    return [(v >> 16) & 0xff, (v >> 8) & 0xff, v & 0xff] as const;
  };
  const [ar, ag, ab] = parseHex(a);
  const [br, bg, bb] = parseHex(b);
  const r = Math.round(ar + (br - ar) * t);
  const g = Math.round(ag + (bg - ag) * t);
  const bv = Math.round(ab + (bb - ab) * t);
  return `rgb(${r},${g},${bv})`;
}

const PULSE_DURATION = 500; // ms
const LABEL_AREA = 28; // px reserved for labels on left and top
const GRID_SIZE = 24;
const TOOLTIP_BG = '#1a1a2e';
const TOOLTIP_BORDER = '#2563eb';

const DensityHeatmap: React.FC<DensityHeatmapProps> = ({
  densityMatrix,
  size = 350,
}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const tooltipRef = useRef<HTMLDivElement>(null);

  // Track previous matrix values to detect changes (for pulse effect)
  const prevMatrixRef = useRef<number[][]>(
    Array.from({ length: GRID_SIZE }, () => new Array(GRID_SIZE).fill(0)),
  );

  // Track pulse timestamps: when a cell last changed
  const pulseTimesRef = useRef<number[][]>(
    Array.from({ length: GRID_SIZE }, () => new Array(GRID_SIZE).fill(0)),
  );

  // Animation frame id for cleanup
  const animFrameRef = useRef<number>(0);

  // Tooltip state (avoid re-renders by using refs)
  const hoverCellRef = useRef<{ row: number; col: number } | null>(null);

  /** Compute cell size from available space */
  const cellSize = (size - LABEL_AREA) / GRID_SIZE;

  /**
   * Detect which cells changed since last frame and record pulse timestamps.
   */
  const detectChanges = useCallback(() => {
    const now = Date.now();
    const prev = prevMatrixRef.current;
    const pulses = pulseTimesRef.current;

    for (let r = 0; r < GRID_SIZE; r++) {
      for (let c = 0; c < GRID_SIZE; c++) {
        const current = densityMatrix[r]?.[c] ?? 0;
        const previous = prev[r][c];
        if (current !== previous) {
          pulses[r][c] = now;
          prev[r][c] = current;
        }
      }
    }
  }, [densityMatrix]);

  /**
   * Main render loop using Canvas + requestAnimationFrame.
   */
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const dpr = window.devicePixelRatio || 1;
    const w = size;
    const h = size;

    // Set canvas resolution for crisp rendering
    if (canvas.width !== w * dpr || canvas.height !== h * dpr) {
      canvas.width = w * dpr;
      canvas.height = h * dpr;
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
      ctx.scale(dpr, dpr);
    }

    // Detect cell changes for pulse animation
    detectChanges();

    const now = Date.now();
    const pulses = pulseTimesRef.current;
    const cs = cellSize;

    // Clear
    ctx.fillStyle = '#0a0a0a';
    ctx.fillRect(0, 0, w, h);

    // Draw chromosome labels - top row
    ctx.font = `${Math.max(8, Math.min(10, cs * 0.7))}px "JetBrains Mono", "Fira Code", monospace`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillStyle = '#6b7280';

    for (let c = 0; c < GRID_SIZE; c++) {
      const x = LABEL_AREA + c * cs + cs / 2;
      ctx.fillText(LABELS[c], x, LABEL_AREA / 2);
    }

    // Draw chromosome labels - left column
    ctx.textAlign = 'right';
    for (let r = 0; r < GRID_SIZE; r++) {
      const y = LABEL_AREA + r * cs + cs / 2;
      ctx.fillText(LABELS[r], LABEL_AREA - 3, y);
    }

    // Draw grid cells
    for (let r = 0; r < GRID_SIZE; r++) {
      for (let c = 0; c < GRID_SIZE; c++) {
        const x = LABEL_AREA + c * cs;
        const y = LABEL_AREA + r * cs;
        const count = densityMatrix[r]?.[c] ?? 0;

        // Diagonal cells (same chromosome) are dimmed/unused
        if (r === c) {
          ctx.fillStyle = '#0a0a0a';
          ctx.fillRect(x, y, cs, cs);
          // Draw a subtle diagonal marker
          ctx.fillStyle = '#1a1a1a';
          ctx.fillRect(x + 1, y + 1, cs - 2, cs - 2);
          continue;
        }

        // Base color from count
        let color = countToColor(count);

        // Pulse effect: brighten cells that changed recently
        const pulseAge = now - pulses[r][c];
        if (pulses[r][c] > 0 && pulseAge < PULSE_DURATION) {
          const pulseFactor = 1 - pulseAge / PULSE_DURATION;
          color = lerpColor(color, '#ffffff', pulseFactor * 0.6);
        }

        ctx.fillStyle = color;
        ctx.fillRect(x + 0.5, y + 0.5, cs - 1, cs - 1);
      }
    }

    // Draw hover highlight
    const hc = hoverCellRef.current;
    if (hc && hc.row !== hc.col) {
      const hx = LABEL_AREA + hc.col * cs;
      const hy = LABEL_AREA + hc.row * cs;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 1.5;
      ctx.strokeRect(hx, hy, cs, cs);

      // Also highlight the symmetric cell
      if (hc.row !== hc.col) {
        const sx = LABEL_AREA + hc.row * cs;
        const sy = LABEL_AREA + hc.col * cs;
        ctx.strokeStyle = 'rgba(255,255,255,0.4)';
        ctx.lineWidth = 1;
        ctx.strokeRect(sx, sy, cs, cs);
      }
    }

    // Continue animation loop
    animFrameRef.current = requestAnimationFrame(draw);
  }, [size, cellSize, densityMatrix, detectChanges]);

  /** Start the animation loop on mount, clean up on unmount */
  useEffect(() => {
    animFrameRef.current = requestAnimationFrame(draw);
    return () => {
      if (animFrameRef.current) {
        cancelAnimationFrame(animFrameRef.current);
      }
    };
  }, [draw]);

  /** Convert mouse position to grid cell */
  const getCellFromEvent = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>): { row: number; col: number } | null => {
      const canvas = canvasRef.current;
      if (!canvas) return null;
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      const col = Math.floor((mx - LABEL_AREA) / cellSize);
      const row = Math.floor((my - LABEL_AREA) / cellSize);

      if (row < 0 || row >= GRID_SIZE || col < 0 || col >= GRID_SIZE) {
        return null;
      }
      return { row, col };
    },
    [cellSize],
  );

  /** Handle mouse move for tooltip */
  const handleMouseMove = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const cell = getCellFromEvent(e);
      const tooltip = tooltipRef.current;
      hoverCellRef.current = cell;

      if (!cell || !tooltip || cell.row === cell.col) {
        if (tooltip) tooltip.style.display = 'none';
        hoverCellRef.current = null;
        return;
      }

      const count = densityMatrix[cell.row]?.[cell.col] ?? 0;
      const chromA = CHROMS[cell.row];
      const chromB = CHROMS[cell.col];
      const labelA = LABELS[cell.row];
      const labelB = LABELS[cell.col];

      // Build tooltip content
      let html = `<div style="font-weight:600;margin-bottom:4px;color:#e5e7eb">chr${labelA} &#8596; chr${labelB}: ${count} pairs</div>`;

      // No additional hotspot data is available from densityMatrix alone,
      // but we show the count prominently
      if (count === 0) {
        html += `<div style="color:#6b7280;font-size:11px">No discordant pairs detected</div>`;
      } else {
        const intensity =
          count <= 5
            ? 'Low'
            : count <= 20
              ? 'Moderate'
              : count <= 50
                ? 'High'
                : count <= 100
                  ? 'Very High'
                  : 'Extreme';
        html += `<div style="color:#9ca3af;font-size:11px">Intensity: ${intensity}</div>`;
      }

      tooltip.innerHTML = html;
      tooltip.style.display = 'block';

      // Position tooltip near cursor, but keep within container
      const canvas = canvasRef.current;
      if (canvas) {
        const rect = canvas.getBoundingClientRect();
        const containerRect = canvas.parentElement?.getBoundingClientRect();
        if (containerRect) {
          let tx = e.clientX - containerRect.left + 12;
          let ty = e.clientY - containerRect.top - 10;

          // Clamp to stay within container
          const tw = tooltip.offsetWidth || 180;
          const th = tooltip.offsetHeight || 60;
          if (tx + tw > containerRect.width) {
            tx = e.clientX - containerRect.left - tw - 8;
          }
          if (ty + th > containerRect.height) {
            ty = containerRect.height - th - 4;
          }
          if (ty < 0) ty = 4;

          tooltip.style.left = `${tx}px`;
          tooltip.style.top = `${ty}px`;
        }
      }
    },
    [densityMatrix, getCellFromEvent],
  );

  const handleMouseLeave = useCallback(() => {
    hoverCellRef.current = null;
    const tooltip = tooltipRef.current;
    if (tooltip) tooltip.style.display = 'none';
  }, []);

  return (
    <div style={{ position: 'relative', width: size, height: size }}>
      <canvas
        ref={canvasRef}
        style={{
          width: size,
          height: size,
          cursor: 'crosshair',
          borderRadius: 4,
        }}
        onMouseMove={handleMouseMove}
        onMouseLeave={handleMouseLeave}
      />
      <div
        ref={tooltipRef}
        style={{
          display: 'none',
          position: 'absolute',
          pointerEvents: 'none',
          background: TOOLTIP_BG,
          border: `1px solid ${TOOLTIP_BORDER}`,
          borderRadius: 6,
          padding: '8px 12px',
          fontSize: 12,
          fontFamily: '"JetBrains Mono", "Fira Code", monospace',
          color: '#e5e7eb',
          zIndex: 100,
          whiteSpace: 'nowrap',
          boxShadow: '0 4px 12px rgba(0,0,0,0.5)',
          maxWidth: 240,
        }}
      />
    </div>
  );
};

export default DensityHeatmap;
