import { useRef, useEffect } from 'react';

interface SparklineProps {
  data: number[];
  width?: number;
  height?: number;
  color?: string;
  filled?: boolean;
}

export default function Sparkline({
  data,
  width = 120,
  height = 20,
  color = '#f59e0b',
  filled = false,
}: SparklineProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || data.length < 2) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, width, height);

    const max = Math.max(...data, 1);
    const step = width / (data.length - 1);

    ctx.beginPath();
    ctx.moveTo(0, height - (data[0] / max) * (height - 2));
    for (let i = 1; i < data.length; i++) {
      ctx.lineTo(i * step, height - (data[i] / max) * (height - 2));
    }

    if (filled) {
      ctx.lineTo((data.length - 1) * step, height);
      ctx.lineTo(0, height);
      ctx.closePath();
      ctx.fillStyle = color + '33';
      ctx.fill();
      // Re-draw line on top
      ctx.beginPath();
      ctx.moveTo(0, height - (data[0] / max) * (height - 2));
      for (let i = 1; i < data.length; i++) {
        ctx.lineTo(i * step, height - (data[i] / max) * (height - 2));
      }
    }

    ctx.strokeStyle = color;
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }, [data, width, height, color, filled]);

  return (
    <canvas
      ref={canvasRef}
      style={{ width, height, display: 'block' }}
    />
  );
}
