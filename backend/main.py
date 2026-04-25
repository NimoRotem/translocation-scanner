"""Translocation Scanner — FastAPI server.

Routes:
  GET  /api/server-files          — list available BAM/CRAM files
  POST /api/scan                  — start a new scan job
  GET  /api/jobs/{job_id}         — poll job status
  GET  /api/jobs/{job_id}/stream  — SSE event stream
  GET  /api/jobs/{job_id}/results — download results
  POST /api/jobs/{job_id}/cancel  — cancel a running scan
  GET  /api/jobs/{job_id}/log     — retrieve per-job log
  GET  /{path}                    — SPA catch-all
"""
from __future__ import annotations

import asyncio
import json
import logging
import multiprocessing as mp
import os
import sys
import threading
import time
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import (
    FileResponse, HTMLResponse, JSONResponse, StreamingResponse,
)
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

# Add backend dir to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from models import ScanJob, JobStatus, ScanStage
from aggregator import EventAggregator
from pipeline_v2 import PipelineV2

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("translocation-scanner")

# --- Configuration ---
ROOT_PATH = os.environ.get("ROOT_PATH", "/translocation-scanner")
PORT = int(os.environ.get("PORT", "8750"))
SAMPLE_DIRS = [
    "/data/aligned_bams",
    "/data/ancestry_app/uploads",
    "/scratch",
]
REFERENCE_PATHS = {
    "GRCh38": "/data/refs/hs38DH.fa",
    "GRCh38_numeric": "/data/genom-nimo/reference.fasta",
}
FRONTEND_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "frontend", "dist")
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
JOBS_FILE = os.path.join(DATA_DIR, "jobs.json")

# --- Global state ---
jobs: dict[str, ScanJob] = {}
aggregators: dict[str, EventAggregator] = {}
cancel_events: dict[str, mp.Event] = {}
_flush_tasks: dict[str, asyncio.Task] = {}


# --- Job persistence ---
def _save_jobs():
    """Persist all jobs to disk."""
    os.makedirs(DATA_DIR, exist_ok=True)
    data = {}
    for jid, job in jobs.items():
        d = job.to_dict()
        d["reference_path"] = job.reference_path
        d["results_dir"] = job.results_dir
        d["validated_calls"] = job.validated_calls
        d["chrom_progress"] = {}
        data[jid] = d
    with open(JOBS_FILE, "w") as f:
        json.dump(data, f, default=str, indent=2)


def _load_jobs():
    """Restore jobs from disk on startup."""
    if not os.path.isfile(JOBS_FILE):
        return
    try:
        with open(JOBS_FILE) as f:
            data = json.load(f)
        for jid, d in data.items():
            job = ScanJob(job_id=jid)
            job.file_path = d.get("file_path", "")
            job.reference_path = d.get("reference_path", "")
            job.reference_build = d.get("reference_build", "GRCh38")
            job.status = JobStatus(d.get("status", "completed"))
            job.stage = ScanStage(d.get("stage", "completed"))
            job.created_at = d.get("created_at", 0)
            job.started_at = d.get("started_at")
            job.completed_at = d.get("completed_at")
            job.error = d.get("error")
            job.total_reads = d.get("total_reads", 0)
            job.reads_processed = d.get("reads_processed", 0)
            job.bytes_processed = d.get("bytes_processed", 0)
            job.discordant_count = d.get("discordant_count", 0)
            job.split_count = d.get("split_count", 0)
            job.clip_count = d.get("clip_count", 0)
            job.chimeric_rate = d.get("chimeric_rate", 0)
            job.insert_size_median = d.get("insert_size_median", 0)
            job.insert_size_std = d.get("insert_size_std", 0)
            job.validated_calls = d.get("validated_calls", [])
            job.results_dir = d.get("results_dir", "")
            jobs[jid] = job
        logger.info("Loaded %d jobs from disk", len(data))
    except Exception:
        logger.exception("Failed to load jobs from disk")


# --- Lifespan ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Translocation Scanner starting on port %d with root_path %s", PORT, ROOT_PATH)
    _load_jobs()
    yield
    _save_jobs()
    for task in _flush_tasks.values():
        task.cancel()
    logger.info("Translocation Scanner shutting down")


app = FastAPI(
    title="Translocation Scanner",
    root_path=ROOT_PATH,
    lifespan=lifespan,
)


# --- Request models ---
class ScanSettings(BaseModel):
    min_mapq: int = 20
    min_clip_length: int = 20
    min_split_aligned: int = 20
    min_pileup_depth: int = 4
    pileup_window: int = 5
    merge_distance: int = 500
    bg_bin_size: int = 100_000
    bg_pvalue_threshold: float = 0.001
    centromere_margin: int = 1_000_000
    skip_clip_realignment: bool = False
    skip_external_callers: bool = False
    parallel_extraction: bool = True
    num_workers: int = 0
    exclude_chrM: bool = True
    top_n_candidates: int = 100
    breakpoint_merge_window: int = 5000
    bg_window_size: int = 10_000
    promiscuous_threshold: int = 5
    min_cluster_support: int = 3
    debug_region_a: Optional[str] = None  # e.g. "chr9:100000000"
    debug_region_b: Optional[str] = None  # e.g. "chr22:23000000"
    debug_margin: int = 2_000_000  # 2MB window around each debug region


class ScanRequest(BaseModel):
    file_path: str
    reference_path: Optional[str] = None
    reference_build: str = "GRCh38"
    settings: Optional[ScanSettings] = None


# --- API Routes ---

@app.get("/api/server-files")
async def list_server_files():
    """List available BAM/CRAM files from configured directories."""
    files = []
    for d in SAMPLE_DIRS:
        if not os.path.isdir(d):
            continue
        try:
            for entry in sorted(os.listdir(d)):
                full = os.path.join(d, entry)
                if not os.path.isfile(full):
                    continue
                ext = entry.lower().rsplit(".", 1)[-1] if "." in entry else ""
                if ext in ("bam", "cram"):
                    # Check for index
                    has_index = (
                        os.path.exists(full + ".bai")
                        or os.path.exists(full + ".crai")
                        or os.path.exists(full.rsplit(".", 1)[0] + ".bai")
                    )
                    size = os.path.getsize(full)
                    files.append({
                        "path": full,
                        "name": entry,
                        "dir": d,
                        "size": size,
                        "size_human": _human_size(size),
                        "format": ext.upper(),
                        "indexed": has_index,
                    })
        except PermissionError:
            continue

    # Also list available references
    refs = []
    for name, path in REFERENCE_PATHS.items():
        if os.path.exists(path):
            refs.append({"name": name, "path": path})

    return {"files": files, "references": refs}


@app.post("/api/scan")
async def start_scan(req: ScanRequest):
    """Start a new translocation scan job."""
    if not os.path.isfile(req.file_path):
        raise HTTPException(400, f"File not found: {req.file_path}")

    # Auto-detect reference if not provided
    ref_path = req.reference_path
    if not ref_path:
        ref_path = REFERENCE_PATHS.get(req.reference_build)
        if not ref_path or not os.path.exists(ref_path):
            # Try to detect from BAM header
            ref_path = _detect_reference(req.file_path)

    if not ref_path or not os.path.exists(ref_path):
        raise HTTPException(400, "No reference FASTA found. Provide reference_path.")

    settings = (req.settings or ScanSettings()).model_dump()
    job = ScanJob(
        file_path=req.file_path,
        reference_path=ref_path,
        reference_build=req.reference_build,
        settings=settings,
    )
    jobs[job.job_id] = job

    # Create aggregator, cancel event, and start flush loop
    agg = EventAggregator()
    aggregators[job.job_id] = agg
    ce = mp.Event()
    cancel_events[job.job_id] = ce
    loop = asyncio.get_event_loop()
    _flush_tasks[job.job_id] = loop.create_task(agg.run_flush_loop())

    # Start pipeline in background thread
    def _run_pipeline():
        orch = PipelineV2(job, event_callback=agg.push_event, cancel_event=ce)
        orch.run()
        # Persist jobs after pipeline completes or fails
        try:
            _save_jobs()
        except Exception:
            logger.exception("Failed to save jobs after pipeline run")

    thread = threading.Thread(target=_run_pipeline, daemon=True, name=f"pipeline-{job.job_id}")
    thread.start()

    logger.info("Started scan job %s for %s", job.job_id, req.file_path)
    return {"job_id": job.job_id, "status": "queued"}


@app.get("/api/jobs")
async def list_jobs():
    """List all jobs, most recent first."""
    job_list = sorted(jobs.values(), key=lambda j: j.created_at, reverse=True)
    return {"jobs": [j.to_dict() for j in job_list]}


@app.get("/api/jobs/{job_id}")
async def get_job(job_id: str):
    """Get job status."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    return job.to_dict()


@app.get("/api/jobs/{job_id}/stream")
async def stream_job(job_id: str, request: Request):
    """SSE event stream for a job."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")

    agg = aggregators.get(job_id)
    if not agg:
        raise HTTPException(404, "No event stream for this job")

    client_id = str(uuid.uuid4())[:8]
    queue = agg.subscribe(client_id)

    async def event_generator():
        try:
            # Send initial state
            yield _sse_format("job.state", job.to_dict())

            while True:
                # Check if client disconnected
                if await request.is_disconnected():
                    break

                try:
                    event = await asyncio.wait_for(queue.get(), timeout=1.0)
                    yield _sse_format(event.get("type", "unknown"), event)
                except asyncio.TimeoutError:
                    # Send keepalive
                    yield ": keepalive\n\n"

                    # Stop streaming after job completes (give 5s for final events)
                    if job.status in (JobStatus.COMPLETED, JobStatus.FAILED):
                        if queue.empty():
                            yield _sse_format("stream.end", {"reason": job.status.value})
                            break
        finally:
            agg.unsubscribe(client_id)

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )


@app.get("/api/jobs/{job_id}/results")
async def get_results(job_id: str):
    """Get results for a completed job."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if job.status != JobStatus.COMPLETED:
        raise HTTPException(400, f"Job not completed (status: {job.status.value})")

    results = {"calls": job.validated_calls, "summary": job.to_dict()}

    # Include file paths if results were written
    if job.results_dir and os.path.isdir(job.results_dir):
        result_files = {}
        for f in os.listdir(job.results_dir):
            result_files[f] = os.path.join(job.results_dir, f)
        results["files"] = {k: f"/api/jobs/{job_id}/download/{k}" for k in result_files}

    return results


@app.get("/api/jobs/{job_id}/report")
async def get_report(job_id: str):
    """Generate a full scan report with interpretation."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")

    report_data = (job.settings or {}).get("_report", {})
    timings = report_data.get("timings", {})
    filter_breakdown = report_data.get("filter_breakdown", {})
    tier_counts = report_data.get("tier_counts", {})
    cluster_counts = report_data.get("cluster_counts", {})
    warnings = report_data.get("warnings", [])
    near_miss_count = report_data.get("near_miss_count", 0)

    raw_clusters = cluster_counts.get("raw_clusters_formed", 0)
    clusters_passing = cluster_counts.get("candidates_retained", 0)

    elapsed = 0.0
    if job.started_at:
        end = job.completed_at or 0
        elapsed = end - job.started_at if end else 0

    # Quality assessment
    chimeric_pct = job.chimeric_rate * 100
    if chimeric_pct < 1.0:
        chimeric_assessment = "normal"
    elif chimeric_pct < 5.0:
        chimeric_assessment = "elevated"
    else:
        chimeric_assessment = "high"

    insert_ok = 200 <= job.insert_size_median <= 600
    insert_assessment = "normal" if insert_ok else "atypical"

    file_name = os.path.basename(job.file_path)

    confirmed = tier_counts.get("confirmed", 0)
    validated = tier_counts.get("validated", 0)
    likely = tier_counts.get("likely", 0)
    strong_candidate = tier_counts.get("strong_candidate", 0)
    candidate = tier_counts.get("candidate", 0)
    filtered = tier_counts.get("filtered", 0)
    high_conf = confirmed + validated + likely
    n_calls = high_conf + strong_candidate + candidate

    # Build interpretation
    if high_conf == 0:
        interp_summary = "No high-confidence translocations detected."
        interp_detail = (
            f"The chimeric rate of {chimeric_pct:.2f}% is "
            f"{'within the normal range (< 1%) expected for standard WGS library preparation' if chimeric_pct < 1 else 'elevated, which may indicate library preparation issues or sample contamination'}. "
            f"All {raw_clusters:,} raw candidate junctions were attributed to background chimerism, "
            f"mapping artifacts, or insufficient evidence after scoring and filtering. "
        )
        if candidate > 0:
            interp_detail += (
                f"{candidate} lower-confidence candidate{'s' if candidate != 1 else ''} "
                f"remain{'s' if candidate == 1 else ''} for review. "
            )
        if near_miss_count > 0:
            interp_detail += (
                f"{near_miss_count} near-miss junction{'s' if near_miss_count != 1 else ''} "
                f"were close to passing but lacked sufficient evidence. "
            )
        interp_detail += (
            "This result is consistent with a normal germline genome without large-scale "
            "chromosomal rearrangements."
        )
    else:
        tier_parts = []
        if confirmed: tier_parts.append(f"{confirmed} confirmed")
        if validated: tier_parts.append(f"{validated} validated")
        if likely: tier_parts.append(f"{likely} likely")
        if candidate: tier_parts.append(f"{candidate} candidate")
        interp_summary = f"{n_calls} translocation call{'s' if n_calls != 1 else ''}: {', '.join(tier_parts)}."
        interp_detail = (
            f"The scan identified {high_conf} high-confidence interchromosomal "
            f"rearrangement{'s' if high_conf != 1 else ''} that passed quality scoring. "
            f"{'Confirmed calls have multi-evidence support (split + discordant) with stringent QC. ' if confirmed else ''}"
            f"{'Validated calls have strong support from 2+ evidence types. ' if validated else ''}"
            f"Review the breakpoint coordinates, score breakdown, and supporting evidence below. "
            f"Filtered junctions ({filtered:,}) were excluded for insufficient support or "
            f"very low mapping quality."
        )

    # Near-miss data for frontend
    near_misses = (job.settings or {}).get("_near_misses", [])

    return {
        "sample": {
            "name": file_name,
            "path": job.file_path,
            "reference_build": job.reference_build,
            "scan_date": job.started_at,
            "elapsed_seconds": round(elapsed, 1),
        },
        "quality": {
            "total_reads": job.total_reads,
            "chimeric_rate": job.chimeric_rate,
            "chimeric_rate_pct": f"{chimeric_pct:.3f}%",
            "chimeric_assessment": chimeric_assessment,
            "insert_size_median": job.insert_size_median,
            "insert_size_std": job.insert_size_std,
            "insert_size_assessment": insert_assessment,
        },
        "evidence": {
            "discordant": job.discordant_count,
            "split": job.split_count,
            "clip_pileups": job.clip_count,
        },
        "pipeline": {
            "clusters_formed": raw_clusters,
            "clusters_passing": clusters_passing,
            "timings": timings,
            "filter_breakdown": filter_breakdown,
        },
        "results": {
            "total_calls": n_calls,
            "by_tier": {
                "confirmed": confirmed,
                "validated": validated,
                "likely": likely,
                "strong_candidate": strong_candidate,
                "candidate": candidate,
            },
            "filtered": filtered,
            "calls": job.validated_calls,
            "near_misses": near_misses[:20],
            "near_miss_count": near_miss_count,
            "mask_manifest_version": report_data.get("mask_manifest_version", ""),
            "reference_build": job.reference_build,
        },
        "interpretation": {
            "summary": interp_summary,
            "detail": interp_detail,
        },
        "warnings": warnings,
    }


@app.post("/api/jobs/{job_id}/cancel")
async def cancel_job(job_id: str):
    """Cancel a running scan job."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if job.status not in (JobStatus.QUEUED, JobStatus.RUNNING):
        raise HTTPException(400, f"Job cannot be cancelled (status: {job.status.value})")
    ce = cancel_events.get(job_id)
    if ce:
        ce.set()
    logger.info("Cancel requested for job %s", job_id)
    return {"status": "cancelling", "job_id": job_id}


@app.get("/api/jobs/{job_id}/log")
async def get_job_log(job_id: str):
    """Retrieve the per-job log file."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    log_path = getattr(job, 'log_path', None)
    if not log_path:
        log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'logs')
        log_path = os.path.join(log_dir, f'{job_id}.log')
    if not os.path.isfile(log_path):
        raise HTTPException(404, "Log file not found")
    return FileResponse(log_path, filename=f'{job_id}.log', media_type='text/plain')


@app.get("/api/jobs/{job_id}/download/{filename}")
async def download_result(job_id: str, filename: str):
    """Download a specific result file."""
    job = jobs.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if not job.results_dir:
        raise HTTPException(404, "No results available")

    # Prevent path traversal
    safe_name = os.path.basename(filename)
    filepath = os.path.join(job.results_dir, safe_name)
    if not os.path.isfile(filepath):
        raise HTTPException(404, f"File not found: {safe_name}")

    return FileResponse(filepath, filename=safe_name)


# --- SPA serving ---

@app.get("/api/health")
async def health():
    return {"status": "ok", "jobs": len(jobs)}


# Mount static files if frontend is built
if os.path.isdir(FRONTEND_DIR):
    app.mount("/assets", StaticFiles(directory=os.path.join(FRONTEND_DIR, "assets")), name="assets")

    @app.get("/{path:path}")
    async def spa_catchall(path: str):
        # Serve index.html for all non-API, non-asset routes
        index = os.path.join(FRONTEND_DIR, "index.html")
        if os.path.isfile(index):
            return HTMLResponse(open(index).read())
        raise HTTPException(404, "Frontend not built")
else:
    @app.get("/")
    async def root():
        return HTMLResponse(
            "<h1>Translocation Scanner</h1>"
            "<p>Frontend not built. Run <code>cd frontend && npm run build</code></p>"
        )


# --- Helpers ---

def _human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} PB"


def _sse_format(event_type: str, data: dict) -> str:
    return f"event: {event_type}\ndata: {json.dumps(data, default=str)}\n\n"


def _detect_reference(bam_path: str) -> Optional[str]:
    """Try to detect correct reference from BAM header chromosome naming."""
    try:
        import pysam
        with pysam.AlignmentFile(bam_path, check_sq=False) as bam:
            refs = bam.references
            if refs and refs[0].startswith("chr"):
                return REFERENCE_PATHS.get("GRCh38")
            else:
                return REFERENCE_PATHS.get("GRCh38_numeric")
    except Exception:
        return REFERENCE_PATHS.get("GRCh38")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=PORT,
        root_path=ROOT_PATH,
        log_level="info",
    )
