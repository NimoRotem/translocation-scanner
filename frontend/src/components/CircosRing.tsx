import { useRef, useEffect, useCallback } from 'react';
import type { ValidatedCall, ChromData } from '../types/events';
import { CHROMS, chromIndex } from '../stores/scanStore';

// ---------------------------------------------------------------------------
// GRCh38 chromosome lengths (Mb, approximate)
// ---------------------------------------------------------------------------
const CHROM_LENGTHS: Record<string, number> = {
  chr1: 249, chr2: 243, chr3: 198, chr4: 191, chr5: 182, chr6: 171,
  chr7: 159, chr8: 146, chr9: 138, chr10: 134, chr11: 135, chr12: 133,
  chr13: 114, chr14: 107, chr15: 102, chr16: 90, chr17: 84, chr18: 80,
  chr19: 59, chr20: 64, chr21: 47, chr22: 51, chrX: 156, chrY: 57,
};

const TOTAL_GENOME_MB = CHROMS.reduce((s, c) => s + (CHROM_LENGTHS[c] ?? 0), 0);

// ---------------------------------------------------------------------------
// Layout constants
// ---------------------------------------------------------------------------
const GAP_DEG = 2;                       // degrees between chromosomes
const TOTAL_GAP_DEG = GAP_DEG * CHROMS.length;
const USABLE_DEG = 360 - TOTAL_GAP_DEG;
const DEG = Math.PI / 180;

// Tier colors for validated calls
const TIER_COLORS: Record<string, string> = {
  confirmed: '#10b981',
  likely: '#3b82f6',
  candidate: '#6b7280',
  filtered: '#4b5563',
};

// ---------------------------------------------------------------------------
// Pre-compute arc layout: start/end angle for each chromosome
// ---------------------------------------------------------------------------
interface ArcLayout {
  chrom: string;
  startAngle: number;   // radians
  endAngle: number;     // radians
  lengthMb: number;
}

function computeLayout(): ArcLayout[] {
  const layouts: ArcLayout[] = [];
  let cursor = -Math.PI / 2; // start at 12 o'clock

  for (const chrom of CHROMS) {
    const mb = CHROM_LENGTHS[chrom] ?? 0;
    const fraction = mb / TOTAL_GENOME_MB;
    const sweep = fraction * USABLE_DEG * DEG;
    layouts.push({
      chrom,
      startAngle: cursor,
      endAngle: cursor + sweep,
      lengthMb: mb,
    });
    cursor += sweep + GAP_DEG * DEG;
  }
  return layouts;
}

const LAYOUT = computeLayout();
const LAYOUT_MAP: Record<string, ArcLayout> = {};
LAYOUT.forEach(l => { LAYOUT_MAP[l.chrom] = l; });

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Map a (chrom, posMb) to an angle on the ring. */
function chromPosToAngle(chrom: string, posMb: number): number {
  const lay = LAYOUT_MAP[chrom];
  if (!lay) return 0;
  const frac = Math.min(Math.max(posMb / lay.lengthMb, 0), 1);
  return lay.startAngle + frac * (lay.endAngle - lay.startAngle);
}

/** Get the midpoint angle of a chromosome arc. */
function chromMidAngle(chrom: string): number {
  const lay = LAYOUT_MAP[chrom];
  if (!lay) return 0;
  return (lay.startAngle + lay.endAngle) / 2;
}

/** Lerp between colors based on t in [0,1]. */
function lerpColor(a: [number, number, number], b: [number, number, number], t: number): string {
  const r = Math.round(a[0] + (b[0] - a[0]) * t);
  const g = Math.round(a[1] + (b[1] - a[1]) * t);
  const bl = Math.round(a[2] + (b[2] - a[2]) * t);
  return `rgb(${r},${g},${bl})`;
}

// Color ramp stops: yellow -> orange -> red
const YELLOW: [number, number, number] = [251, 191, 36];   // #fbbf24
const ORANGE: [number, number, number] = [245, 158, 11];   // #f59e0b
const RED: [number, number, number]    = [239, 68, 68];    // #ef4444

function countToColor(count: number, maxCount: number): string {
  if (maxCount <= 0) return lerpColor(YELLOW, YELLOW, 0);
  const t = Math.min(count / maxCount, 1);
  if (t < 0.5) {
    return lerpColor(YELLOW, ORANGE, t * 2);
  }
  return lerpColor(ORANGE, RED, (t - 0.5) * 2);
}

// ---------------------------------------------------------------------------
// Props
// ---------------------------------------------------------------------------
interface CircosRingProps {
  chromProgress: Record<string, ChromData>;
  provisionalArcs: Array<{
    chrom_a: string; pos_a: number;
    chrom_b: string; pos_b: number;
    count: number; timestamp: number;
  }>;
  validatedCalls: ValidatedCall[];
  mode: 'idle' | 'streaming' | 'validated';
  size?: number;
}

// ---------------------------------------------------------------------------
// Component
// ---------------------------------------------------------------------------
function CircosRing({
  chromProgress,
  provisionalArcs,
  validatedCalls,
  mode,
  size = 500,
}: CircosRingProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const rafRef = useRef<number>(0);
  const dprRef = useRef<number>(1);

  // Stable refs so the animation loop always sees latest props without
  // re-registering the effect.
  const propsRef = useRef({ chromProgress, provisionalArcs, validatedCalls, mode, size });
  propsRef.current = { chromProgress, provisionalArcs, validatedCalls, mode, size };

  // ------------------------------------------------------------------
  // Drawing
  // ------------------------------------------------------------------
  const draw = useCallback((now: number) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const { chromProgress: cp, provisionalArcs: arcs, validatedCalls: calls, mode: m, size: sz } = propsRef.current;
    const dpr = dprRef.current;
    const w = sz;
    const cx = w / 2;
    const cy = w / 2;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    ctx.save();
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, w);

    // Radii
    const outerR = cx * 0.88;
    const ringWidth = cx * 0.09;
    const innerR = outerR - ringWidth;
    const labelR = outerR + cx * 0.06;
    const arcR = innerR - cx * 0.03; // radius for inner bezier endpoints

    // ----------------------------------------------------------------
    // 1. Outer chromosome ring
    // ----------------------------------------------------------------
    for (const lay of LAYOUT) {
      const chromData = cp[lay.chrom];
      const pct = chromData?.pct ?? 0;

      // Background arc (dark)
      ctx.beginPath();
      ctx.arc(cx, cy, outerR - ringWidth / 2, lay.startAngle, lay.endAngle);
      ctx.lineWidth = ringWidth;
      ctx.strokeStyle = 'rgba(55, 65, 81, 0.6)'; // gray-700
      ctx.lineCap = 'butt';
      ctx.stroke();

      // Progress fill (amber gradient)
      if (pct > 0) {
        const fillEnd = lay.startAngle + (pct / 100) * (lay.endAngle - lay.startAngle);
        const grad = ctx.createConicGradient(lay.startAngle, cx, cy);
        // Conic gradient is relative to startAngle; we approximate with linear
        ctx.beginPath();
        ctx.arc(cx, cy, outerR - ringWidth / 2, lay.startAngle, fillEnd);
        ctx.lineWidth = ringWidth - 2;
        ctx.strokeStyle = pct >= 100
          ? 'rgba(16, 185, 129, 0.85)'   // green when complete
          : 'rgba(245, 158, 11, 0.75)';  // amber while scanning
        ctx.lineCap = 'butt';
        ctx.stroke();
      }

      // Label
      const midA = (lay.startAngle + lay.endAngle) / 2;
      const lx = cx + labelR * Math.cos(midA);
      const ly = cy + labelR * Math.sin(midA);
      ctx.save();
      ctx.translate(lx, ly);
      // Rotate text to be readable
      let textAngle = midA + Math.PI / 2;
      // Flip if on bottom half so text isn't upside-down
      if (midA > Math.PI / 2 && midA < Math.PI * 1.5) {
        textAngle += Math.PI;
      }
      // Handle wrapping around: angles can be negative
      const normMid = ((midA % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
      if (normMid > Math.PI * 0.5 && normMid < Math.PI * 1.5) {
        textAngle = midA - Math.PI / 2;
      }
      ctx.rotate(textAngle);
      ctx.fillStyle = 'rgba(209, 213, 219, 0.8)'; // gray-300
      ctx.font = `${Math.max(8, w * 0.018)}px ui-monospace, monospace`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      // Shorten label: chr1 -> 1, chrX -> X
      const shortLabel = lay.chrom.replace('chr', '');
      ctx.fillText(shortLabel, 0, 0);
      ctx.restore();
    }

    // ----------------------------------------------------------------
    // 2. Provisional arcs (during streaming)
    // ----------------------------------------------------------------
    if (m === 'streaming' && arcs.length > 0) {
      // Find max count for color ramp and top 1% threshold
      const counts = arcs.map(a => a.count);
      const maxCount = Math.max(...counts, 1);
      const sorted = [...counts].sort((a, b) => b - a);
      const top1pctThreshold = sorted[Math.max(0, Math.floor(sorted.length * 0.01))] ?? maxCount;

      ctx.save();
      ctx.globalCompositeOperation = 'lighter';

      for (const arc of arcs) {
        const age = (now - arc.timestamp) / 1000; // seconds
        // Fade arcs older than 2 seconds
        let fadeFactor = 1;
        if (age > 2) {
          fadeFactor = Math.max(0, 1 - (age - 2) / 3); // fully gone after 5s total
        }
        if (fadeFactor <= 0) continue;

        const baseOpacity = 0.15 + 0.65 * Math.min(arc.count / maxCount, 1);
        const opacity = baseOpacity * fadeFactor;

        // Endpoints
        const angleA = arc.pos_a > 0
          ? chromPosToAngle(arc.chrom_a, arc.pos_a / 1e6)
          : chromMidAngle(arc.chrom_a);
        const angleB = arc.pos_b > 0
          ? chromPosToAngle(arc.chrom_b, arc.pos_b / 1e6)
          : chromMidAngle(arc.chrom_b);

        const x1 = cx + arcR * Math.cos(angleA);
        const y1 = cy + arcR * Math.sin(angleA);
        const x2 = cx + arcR * Math.cos(angleB);
        const y2 = cy + arcR * Math.sin(angleB);

        // Control point: toward center, weighted by angular distance
        const cpX = cx;
        const cpY = cy;

        // Pulsing for top 1%
        let lineW = 1 + 2 * Math.min(arc.count / maxCount, 1);
        if (arc.count >= top1pctThreshold) {
          const pulse = 0.5 + 0.5 * Math.sin(now / 1000 * 2 * Math.PI); // 1 Hz
          lineW += 2 * pulse;
        }

        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.quadraticCurveTo(cpX, cpY, x2, y2);
        ctx.strokeStyle = countToColor(arc.count, maxCount);
        ctx.globalAlpha = opacity;
        ctx.lineWidth = lineW;
        ctx.stroke();
      }

      ctx.globalAlpha = 1;
      ctx.globalCompositeOperation = 'source-over';
      ctx.restore();
    }

    // ----------------------------------------------------------------
    // 3. Validated arcs (post-scan)
    // ----------------------------------------------------------------
    if (m === 'validated' && calls.length > 0) {
      for (const call of calls) {
        const color = TIER_COLORS[call.tier] ?? TIER_COLORS.candidate;

        const angleA = chromPosToAngle(call.chrom_a, call.pos_a / 1e6);
        const angleB = chromPosToAngle(call.chrom_b, call.pos_b / 1e6);

        const x1 = cx + arcR * Math.cos(angleA);
        const y1 = cy + arcR * Math.sin(angleA);
        const x2 = cx + arcR * Math.cos(angleB);
        const y2 = cy + arcR * Math.sin(angleB);

        // Solid bezier through center
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.quadraticCurveTo(cx, cy, x2, y2);
        ctx.strokeStyle = color;
        ctx.lineWidth = 2.5;
        ctx.globalAlpha = 0.9;
        ctx.stroke();
        ctx.globalAlpha = 1;

        // Label with event_id at midpoint of the bezier
        // Midpoint of quadratic bezier: P(0.5) = 0.25*P0 + 0.5*Cp + 0.25*P1
        const mx = 0.25 * x1 + 0.5 * cx + 0.25 * x2;
        const my = 0.25 * y1 + 0.5 * cy + 0.25 * y2;

        ctx.fillStyle = color;
        ctx.font = `bold ${Math.max(8, w * 0.016)}px ui-monospace, monospace`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(call.event_id, mx, my);
      }
    }

    // ----------------------------------------------------------------
    // 4. Center info (idle / streaming)
    // ----------------------------------------------------------------
    if (m === 'idle') {
      ctx.fillStyle = 'rgba(156, 163, 175, 0.4)';
      ctx.font = `${w * 0.028}px ui-monospace, monospace`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('AWAITING SCAN', cx, cy);
    }

    ctx.restore();
  }, []);

  // ------------------------------------------------------------------
  // Animation loop
  // ------------------------------------------------------------------
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    // Set up DPR scaling
    const dpr = window.devicePixelRatio || 1;
    dprRef.current = dpr;
    const sz = propsRef.current.size;
    canvas.width = sz * dpr;
    canvas.height = sz * dpr;
    canvas.style.width = `${sz}px`;
    canvas.style.height = `${sz}px`;

    let running = true;

    function loop(ts: number) {
      if (!running) return;
      draw(ts);
      rafRef.current = requestAnimationFrame(loop);
    }

    rafRef.current = requestAnimationFrame(loop);

    return () => {
      running = false;
      cancelAnimationFrame(rafRef.current);
    };
  }, [draw]);

  // Resize canvas when size prop changes
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const dpr = window.devicePixelRatio || 1;
    dprRef.current = dpr;
    canvas.width = size * dpr;
    canvas.height = size * dpr;
    canvas.style.width = `${size}px`;
    canvas.style.height = `${size}px`;
  }, [size]);

  return (
    <canvas
      ref={canvasRef}
      style={{
        width: size,
        height: size,
        maxWidth: '100%',
        aspectRatio: '1 / 1',
      }}
    />
  );
}

export default CircosRing;
