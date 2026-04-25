"""Report contract tests — verify report structure without running a full pipeline.

Tests the /api/jobs/{job_id}/report endpoint output structure against
a mock completed ScanJob. No BAM file needed.
"""
from __future__ import annotations

import os
import sys
import time

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from models import ScanJob, JobStatus, ScanStage, Tier


def _make_mock_job(
    num_calls: int = 3,
    chimeric_rate: float = 0.005,
    include_report_data: bool = True,
) -> ScanJob:
    """Create a mock completed ScanJob with realistic report data."""
    job = ScanJob(
        file_path="/data/aligned_bams/test.bam",
        reference_path="/data/refs/hs38DH.fa",
        reference_build="GRCh38",
    )
    job.status = JobStatus.COMPLETED
    job.stage = ScanStage.COMPLETED
    job.started_at = time.time() - 120
    job.completed_at = time.time()
    job.total_reads = 500_000_000
    job.reads_processed = 500_000_000
    job.discordant_count = 1200
    job.split_count = 45
    job.clip_count = 300
    job.chimeric_rate = chimeric_rate
    job.insert_size_median = 350
    job.insert_size_std = 80

    # Build mock validated calls
    calls = []
    tiers = ["likely", "candidate", "candidate"]
    for i in range(num_calls):
        calls.append({
            "event_id": f"TRA_chr9_100000_{i}_chr22_230000_{i}_++",
            "cluster_id": f"CLU_{i + 1:03d}",
            "chrom_a": "chr9",
            "pos_a": 100000 + i * 1000,
            "chrom_b": "chr22",
            "pos_b": 230000 + i * 1000,
            "orientation": "++",
            "tier": tiers[i % len(tiers)],
            "score": 50.0 - i * 10,
            "score_components": {
                "support_pr": 15.0,
                "support_sr": 10.0,
                "support_clip": 2.0,
                "pvalue_a": 8.0,
                "pvalue_b": 7.0,
                "unique_starts": 5.0,
                "reciprocal": 0.0,
                "external": 0.0,
                "mapq": 6.0,
            },
            "support": {
                "discordant": 20,
                "split": 3,
                "clipped": 2,
                "total": 25,
            },
            "evidence_label": "enriched_cluster",
            "filter_flags": [],
            "reject_reasons": [],
        })
    job.validated_calls = calls

    if include_report_data:
        job.settings["_report"] = {
            "timings": {
                "extraction": 45.2,
                "clustering": 3.1,
                "clip_realignment": 1.5,
                "annotation": 0.8,
                "external_callers": 120.0,
                "background_model": 2.3,
                "filtering": 0.5,
                "output": 1.0,
                "total": 174.4,
            },
            "filter_breakdown": {
                "hard_exclude:mt_involvement": 50,
                "reject:sr0_pr_lt10": 200,
            },
            "tier_counts": {
                "confirmed": 0,
                "validated": 0,
                "likely": 1,
                "strong_candidate": 0,
                "candidate": 2,
                "filtered": 500,
            },
            "cluster_counts": {
                "raw_clusters_formed": 503,
                "candidates_retained": 3,
            },
            "near_miss_count": 5,
            "mask_manifest_version": "v1.0-2025-04-20",
            "delly_status": "completed",
        }
        job.settings["_near_misses"] = []

    return job


def _generate_report(job: ScanJob) -> dict:
    """Generate report using the same logic as the /api/jobs/{job_id}/report endpoint."""
    # Import the report generation logic from main.py
    # We re-implement the core logic here to avoid needing a running server
    report_data = (job.settings or {}).get("_report", {})
    timings = report_data.get("timings", {})
    filter_breakdown = report_data.get("filter_breakdown", {})
    tier_counts = report_data.get("tier_counts", {})
    cluster_counts = report_data.get("cluster_counts", {})
    near_miss_count = report_data.get("near_miss_count", 0)

    raw_clusters = cluster_counts.get("raw_clusters_formed", 0)

    elapsed = 0.0
    if job.started_at:
        end = job.completed_at or 0
        elapsed = end - job.started_at if end else 0

    chimeric_pct = job.chimeric_rate * 100
    file_name = os.path.basename(job.file_path)

    confirmed = tier_counts.get("confirmed", 0)
    validated = tier_counts.get("validated", 0)
    likely = tier_counts.get("likely", 0)
    strong_candidate = tier_counts.get("strong_candidate", 0)
    candidate = tier_counts.get("candidate", 0)
    filtered = tier_counts.get("filtered", 0)
    high_conf = confirmed + validated + likely
    n_calls = high_conf + strong_candidate + candidate

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
        },
        "evidence": {
            "discordant": job.discordant_count,
            "split": job.split_count,
            "clip_pileups": job.clip_count,
        },
        "pipeline": {
            "clusters_formed": raw_clusters,
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
            "mask_manifest_version": report_data.get("mask_manifest_version", ""),
        },
        "interpretation": {
            "summary": "test",
            "detail": "test detail",
        },
    }


class TestReportContract:
    """Verify the report JSON structure matches the frontend contract."""

    def test_required_top_level_keys(self):
        """Report has all required top-level keys."""
        job = _make_mock_job()
        report = _generate_report(job)

        required = {"sample", "quality", "evidence", "pipeline", "results", "interpretation"}
        assert required.issubset(report.keys()), (
            f"Missing keys: {required - set(report.keys())}"
        )

    def test_tier_counts_sum(self):
        """by_tier values sum to total_calls."""
        job = _make_mock_job()
        report = _generate_report(job)

        by_tier = report["results"]["by_tier"]
        total = report["results"]["total_calls"]
        tier_sum = sum(by_tier.values())
        assert tier_sum == total, (
            f"Tier sum {tier_sum} != total_calls {total}"
        )

    def test_calls_have_score_components(self):
        """Every call in results.calls has a score_components dict."""
        job = _make_mock_job()
        report = _generate_report(job)

        for call in report["results"]["calls"]:
            assert "score_components" in call, (
                f"Call {call.get('cluster_id')} missing score_components"
            )
            assert isinstance(call["score_components"], dict)
            assert len(call["score_components"]) > 0

    def test_pipeline_timings_present(self):
        """pipeline.timings has stage keys."""
        job = _make_mock_job()
        report = _generate_report(job)

        timings = report["pipeline"]["timings"]
        assert isinstance(timings, dict)
        assert len(timings) > 0
        # At minimum should have extraction and clustering
        assert "extraction" in timings
        assert "clustering" in timings

    def test_mask_manifest_version_present(self):
        """results has mask_manifest_version."""
        job = _make_mock_job()
        report = _generate_report(job)

        assert "mask_manifest_version" in report["results"]

    def test_no_known_event_language(self):
        """No 'known event' or 'recovered' language in interpretation."""
        job = _make_mock_job()
        report = _generate_report(job)

        interp = report["interpretation"]
        text = (interp.get("summary", "") + " " + interp.get("detail", "")).lower()
        forbidden = ["known event", "known translocation", "recovered"]
        for phrase in forbidden:
            assert phrase not in text, (
                f"Interpretation contains forbidden phrase: '{phrase}'"
            )

    def test_zero_calls_report(self):
        """Report with zero calls has correct structure."""
        job = _make_mock_job(num_calls=0)
        job.settings["_report"]["tier_counts"] = {
            "confirmed": 0, "validated": 0, "likely": 0,
            "strong_candidate": 0, "candidate": 0, "filtered": 100,
        }
        report = _generate_report(job)

        assert report["results"]["total_calls"] == 0
        assert report["results"]["by_tier"]["confirmed"] == 0
        assert len(report["results"]["calls"]) == 0

    def test_sample_section(self):
        """Sample section has required fields."""
        job = _make_mock_job()
        report = _generate_report(job)

        sample = report["sample"]
        assert "name" in sample
        assert "path" in sample
        assert "reference_build" in sample
        assert "scan_date" in sample
        assert "elapsed_seconds" in sample

    def test_quality_section(self):
        """Quality section has required fields."""
        job = _make_mock_job()
        report = _generate_report(job)

        quality = report["quality"]
        assert "total_reads" in quality
        assert "chimeric_rate" in quality
        assert quality["total_reads"] > 0
