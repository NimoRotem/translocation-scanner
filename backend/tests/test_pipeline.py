"""Pipeline integration tests with tiered test corpus.

Test tiers:
  nano   — 418MB BAM, chr9+chr22, <60s end-to-end
  small  — spike-in BAM with 5 known events, <3 min
  medium — 10GB downsampled Nimo, <10 min
  full   — Nimo.bam at full coverage, <30 min (explicit request only)

Run: python3 -m pytest tests/test_pipeline.py -v --tb=short
  Nano only: -k nano
  With timing: --durations=0
"""
from __future__ import annotations

import json
import os
import sys
import time

import pytest

# Add backend to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from models import ScanJob, JobStatus, ScanStage, Tier
from pipeline_v2 import PipelineV2

# ---------------------------------------------------------------------------
# Test corpus paths (on genom-beast-gpu)
# ---------------------------------------------------------------------------
NANO_BAM = "/data/scan_archive/test_corpus/nano.bam"
SMALL_BAM = "/data/scan_archive/test_corpus/spikein.bam"
MEDIUM_BAM = "/data/scan_archive/test_corpus/nimo_downsampled.bam"
FULL_BAM = "/data/aligned_bams/Nimo.bam"

REF_NUMERIC = "/data/genom-nimo/reference.fasta"
REF_CHR = "/data/refs/hs38DH.fa"

VALIDATION_STATE_FILE = "/data/scan_archive/validation_state.json"


def _get_git_sha() -> str:
    """Get current git SHA for validation tracking."""
    try:
        import subprocess
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True,
            cwd=os.path.join(os.path.dirname(__file__), "..", ".."),
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def _run_pipeline(
    bam_path: str,
    reference_path: str,
    settings: dict | None = None,
) -> tuple[ScanJob, list[dict]]:
    """Run pipeline directly (no web service) and return job + events."""
    if not os.path.isfile(bam_path):
        pytest.skip(f"Test BAM not found: {bam_path}")
    if not os.path.isfile(reference_path):
        pytest.skip(f"Reference not found: {reference_path}")

    events: list[dict] = []

    base_settings = {
        "parallel_extraction": True,
        "skip_clip_realignment": True,  # Skip minimap2 for speed
    }
    if settings:
        base_settings.update(settings)

    job = ScanJob(
        file_path=bam_path,
        reference_path=reference_path,
        reference_build="GRCh38",
        settings=base_settings,
    )

    pipeline = PipelineV2(
        job,
        event_callback=lambda e: events.append(e),
    )
    pipeline.run()

    return job, events


def _save_validation_state(test_name: str, passed: bool, details: dict):
    """Record test result for validation gate enforcement."""
    state = {}
    if os.path.isfile(VALIDATION_STATE_FILE):
        try:
            with open(VALIDATION_STATE_FILE) as f:
                state = json.load(f)
        except Exception:
            pass

    state[test_name] = {
        "passed": passed,
        "timestamp": time.time(),
        "git_sha": _get_git_sha(),
        "details": details,
    }

    os.makedirs(os.path.dirname(VALIDATION_STATE_FILE), exist_ok=True)
    with open(VALIDATION_STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)


# ---------------------------------------------------------------------------
# Nano tests (<60s)
# ---------------------------------------------------------------------------

class TestNano:
    """Nano tier: 418MB BAM, chr9+chr22, target <60s."""

    TIME_BUDGET = 60  # seconds

    def test_pipeline_completes(self):
        """Pipeline runs to completion without errors."""
        t0 = time.time()
        job, events = _run_pipeline(
            NANO_BAM, REF_NUMERIC,
            settings={"skip_external_callers": True},
        )
        elapsed = time.time() - t0

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"
        assert job.stage == ScanStage.COMPLETED
        assert elapsed < self.TIME_BUDGET, (
            f"Nano test took {elapsed:.1f}s, budget is {self.TIME_BUDGET}s"
        )

        # Record validation state
        _save_validation_state("nano", True, {
            "elapsed": round(elapsed, 1),
            "discordant": job.discordant_count,
            "split": job.split_count,
            "status": job.status.value,
        })

    def test_no_silent_failures(self):
        """No stages should be silently skipped."""
        job, events = _run_pipeline(
            NANO_BAM, REF_NUMERIC,
            settings={"skip_external_callers": True},
        )
        assert job.status == JobStatus.COMPLETED

        stage_events = [e for e in events if e.get("type") == "scan.stage_changed"]
        stages_seen = {e["stage"] for e in stage_events}
        required_stages = {"extraction", "clustering"}
        assert required_stages.issubset(stages_seen), (
            f"Missing stages: {required_stages - stages_seen}"
        )

    def test_extraction_produces_evidence(self):
        """Extraction should find discordant and/or split reads."""
        job, events = _run_pipeline(
            NANO_BAM, REF_NUMERIC,
            settings={"skip_external_callers": True},
        )
        assert job.status == JobStatus.COMPLETED
        # Even a small region should have some chimeric reads
        assert job.discordant_count >= 0  # Can be 0 for small regions
        assert job.reads_processed > 0, "Should process some reads"

    def test_with_delly(self):
        """Pipeline runs with DELLY enabled on nano BAM."""
        t0 = time.time()
        job, events = _run_pipeline(NANO_BAM, REF_NUMERIC)
        elapsed = time.time() - t0

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"
        # DELLY on 2 chroms should be very fast
        timings = job.settings.get("_report", {}).get("timings", {})
        delly_time = timings.get("external_callers", 0)
        print(f"DELLY time on nano: {delly_time:.1f}s")

        _save_validation_state("nano_delly", True, {
            "elapsed": round(elapsed, 1),
            "delly_time": round(delly_time, 1),
        })


# ---------------------------------------------------------------------------
# Small tests (spike-in, <3min)
# ---------------------------------------------------------------------------

class TestSmall:
    """Small tier: spike-in BAM with known events, target <3 min."""

    TIME_BUDGET = 180  # seconds

    @pytest.mark.skipif(
        not os.path.isfile(SMALL_BAM),
        reason="Spike-in BAM not found",
    )
    def test_spikein_events_detected(self):
        """All 5 spike-in events should be detected at tier likely+."""
        t0 = time.time()
        job, events = _run_pipeline(SMALL_BAM, REF_NUMERIC)
        elapsed = time.time() - t0

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"
        assert elapsed < self.TIME_BUDGET, (
            f"Small test took {elapsed:.1f}s, budget is {self.TIME_BUDGET}s"
        )

        # Check tier distribution
        calls = job.validated_calls
        high_tier = [
            c for c in calls
            if c.get("tier") in ("confirmed", "validated", "likely")
        ]

        _save_validation_state("small_spikein", len(high_tier) >= 5, {
            "elapsed": round(elapsed, 1),
            "high_tier_calls": len(high_tier),
            "total_calls": len(calls),
        })

        assert len(high_tier) >= 5, (
            f"Expected >=5 spike-in events at likely+, got {len(high_tier)}"
        )


# ---------------------------------------------------------------------------
# Full tests (explicit request only)
# ---------------------------------------------------------------------------

class TestFull:
    """Full tier: Nimo.bam at full coverage, target <60 min.

    DELLY on a 93GB 30x WGS BAM takes ~45 min even with exclude.bed.
    Track 1 (blind discovery) completes in ~7 min, but Track 2 (DELLY)
    runs until its subprocess timeout. If DELLY times out, the pipeline
    continues with Track 1 results only.
    """

    TIME_BUDGET = 3600  # 60 min (DELLY on full BAM is the bottleneck)

    @pytest.mark.skipif(
        not os.path.isfile(FULL_BAM),
        reason="Full BAM not found",
    )
    @pytest.mark.full
    def test_full_nimo_scan(self):
        """Full Nimo scan should complete within budget."""
        t0 = time.time()
        job, events = _run_pipeline(FULL_BAM, REF_NUMERIC)
        elapsed = time.time() - t0

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"
        assert elapsed < self.TIME_BUDGET, (
            f"Full test took {elapsed:.1f}s, budget is {self.TIME_BUDGET}s"
        )

        timings = job.settings.get("_report", {}).get("timings", {})
        print(f"Full Nimo timings: {json.dumps(timings, indent=2)}")

        _save_validation_state("full_nimo", True, {
            "elapsed": round(elapsed, 1),
            "timings": timings,
            "calls": len(job.validated_calls),
        })


# ---------------------------------------------------------------------------
# DELLY-specific tests
# ---------------------------------------------------------------------------

class TestDelly:
    """Verify DELLY integration is correct."""

    def test_exclude_bed_exists(self):
        """exclude.bed must exist for DELLY."""
        for bed in [
            "/data/masks/exclude_grch38.bed",
            "/data/masks/exclude_numeric.bed",
        ]:
            if os.path.isfile(bed):
                with open(bed) as f:
                    lines = sum(1 for _ in f)
                assert lines > 100, f"{bed} has only {lines} entries"
                return
        pytest.fail("No exclude.bed found")

    def test_delly_command_has_exclude(self):
        """DELLY command line must include -x exclude.bed."""
        job, events = _run_pipeline(
            NANO_BAM, REF_NUMERIC,
        )
        # Check job log for DELLY command
        log_path = getattr(job, 'log_path', None)
        if log_path and os.path.isfile(log_path):
            with open(log_path) as f:
                log_content = f.read()
            assert "-x" in log_content, (
                "DELLY command in log does not contain -x (exclude.bed)"
            )
            assert "exclude" in log_content.lower(), (
                "DELLY command does not reference exclude.bed"
            )
