"""Negative control — run pipeline on a sample with no known translocations.

Asserts:
  - Pipeline completes without error.
  - Zero confirmed/likely tier calls (noise-floor candidates are acceptable).

Requires nano.bam on the test server. Skipped if not available.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

NANO_BAM = "/data/scan_archive/test_corpus/nano.bam"

# Reference path — not strictly needed when skipping external callers,
# but the pipeline may use it for header validation.
REFERENCE = ""


@pytest.mark.skipif(
    not os.path.isfile(NANO_BAM),
    reason=f"Negative control BAM not found: {NANO_BAM}",
)
class TestNegativeControl:
    """Run pipeline on nano.bam and verify no high-confidence calls."""

    def test_no_confirmed_or_likely_calls(self):
        """Pipeline should produce zero confirmed/likely calls on a normal sample."""
        from models import ScanJob
        from pipeline_v2 import PipelineV2

        job = ScanJob(
            file_path=NANO_BAM,
            reference_path=REFERENCE,
            reference_build="GRCh38",
            settings={
                "skip_external_callers": True,
                "skip_clip_realignment": True,
                "parallel_extraction": False,
                "min_cluster_support": 2,
            },
        )

        events: list[dict] = []
        pipeline = PipelineV2(
            job,
            event_callback=lambda e: events.append(e),
        )
        pipeline.run()

        from models import JobStatus

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"

        # Check calls — confirmed and likely tiers should be absent
        HIGH_TIERS = {"confirmed", "likely", "CONFIRMED", "LIKELY"}
        high_confidence = [
            c for c in (job.validated_calls or [])
            if c.get("tier", "").lower() in {"confirmed", "likely"}
        ]

        if high_confidence:
            for c in high_confidence:
                print(
                    f"  Unexpected high-confidence call: "
                    f"{c.get('chrom_a')}:{c.get('pos_a')} -> "
                    f"{c.get('chrom_b')}:{c.get('pos_b')} "
                    f"(tier={c.get('tier')}, score={c.get('score')})"
                )

        assert len(high_confidence) == 0, (
            f"Expected 0 confirmed/likely calls on negative control, "
            f"got {len(high_confidence)}"
        )

    def test_candidate_scores_below_threshold(self):
        """Any candidate-tier calls should have low scores (noise floor)."""
        from models import ScanJob
        from pipeline_v2 import PipelineV2

        job = ScanJob(
            file_path=NANO_BAM,
            reference_path=REFERENCE,
            reference_build="GRCh38",
            settings={
                "skip_external_callers": True,
                "skip_clip_realignment": True,
                "parallel_extraction": False,
                "min_cluster_support": 2,
            },
        )

        pipeline = PipelineV2(
            job,
            event_callback=lambda _: None,
        )
        pipeline.run()

        from models import JobStatus

        assert job.status == JobStatus.COMPLETED, f"Pipeline failed: {job.error}"

        # Candidate calls are acceptable but should have modest scores
        MAX_NOISE_SCORE = 50.0
        candidates = [
            c for c in (job.validated_calls or [])
            if c.get("tier", "").lower() in {"candidate", "strong_candidate"}
        ]

        high_scoring = [c for c in candidates if c.get("score", 0) > MAX_NOISE_SCORE]
        if high_scoring:
            for c in high_scoring:
                print(
                    f"  High-scoring candidate on negative control: "
                    f"{c.get('chrom_a')}:{c.get('pos_a')} -> "
                    f"{c.get('chrom_b')}:{c.get('pos_b')} "
                    f"(tier={c.get('tier')}, score={c.get('score')})"
                )

        assert len(high_scoring) == 0, (
            f"Expected candidate scores <= {MAX_NOISE_SCORE} on negative control, "
            f"got {len(high_scoring)} above threshold"
        )
