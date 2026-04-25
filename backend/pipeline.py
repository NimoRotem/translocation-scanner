"""Pipeline orchestrator — runs all stages in a background thread."""
from __future__ import annotations
import logging
import os
import subprocess
import tempfile
import threading
import time
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Optional

from models import (
    ScanJob, JobStatus, ScanStage, EvidenceCluster, SVRead, ClipPileup, Tier
)

logger = logging.getLogger(__name__)


class CancelledException(Exception):
    """Raised when a scan is cancelled by the user."""
    pass


class PipelineOrchestrator:
    """Runs the full translocation scanning pipeline in a background thread.

    Stages:
    1. SV-read extraction with telemetry
    2. Clustering
    3. Clip realignment (minimap2)
    4. Background model scoring
    5. Filtering (mostly soft flags)
    6. Scoring and tier assignment
    7. Breakpoint window aggregation
    8. Reciprocal consolidation
    9. Top-N selection + near-miss assembly
    10. Output writing
    """

    def __init__(self, job: ScanJob, event_callback=None, cancel_event=None):
        self.job = job
        self.emit = event_callback or (lambda e: None)
        self._chrom_lengths: dict[str, int] = {}
        self._all_reads: list[SVRead] = []
        self.settings = job.settings or {}
        self._cancel_event = cancel_event
        self._job_log: Optional[logging.Logger] = None
        self._job_log_handler: Optional[logging.Handler] = None
        self._report_data: dict = {"timings": {}, "filter_breakdown": {}}
        self._debug_evidence: list[dict] = []

    def _check_cancel(self):
        if self._cancel_event is not None and self._cancel_event.is_set():
            raise CancelledException("Scan cancelled by user")

    def _setup_job_log(self):
        log_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..', 'logs'
        )
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, f'{self.job.job_id}.log')
        handler = logging.FileHandler(log_path)
        handler.setFormatter(logging.Formatter(
            '%(asctime)s [%(levelname)s] %(message)s'
        ))
        self._job_log = logging.getLogger(f'pipeline.{self.job.job_id}')
        self._job_log.addHandler(handler)
        self._job_log.setLevel(logging.INFO)
        self._job_log_handler = handler
        self.job.log_path = log_path
        return log_path

    def _log(self, msg: str, *args, level: int = logging.INFO):
        logger.log(level, msg, *args)
        if self._job_log:
            self._job_log.log(level, msg, *args)

    def run(self):
        """Execute the full pipeline. Called from a thread."""
        log_path = self._setup_job_log()
        try:
            self.job.status = JobStatus.RUNNING
            self.job.started_at = time.time()
            self._log("Pipeline started for job %s", self.job.job_id)
            self._log("Input: %s", self.job.file_path)
            self._log("Reference: %s (%s)", self.job.reference_path, self.job.reference_build)
            self._log("Settings: %s", self.settings)
            self.emit({"type": "scan.started", "job_id": self.job.job_id,
                        "file_path": self.job.file_path})

            # Stage 1+2: Extract SV reads
            t0 = time.time()
            self._stage_extraction()
            dt = time.time() - t0
            self._report_data["timings"]["extraction"] = round(dt, 1)
            self._log("Extraction: %.1fs — %d disc, %d split, %d clip pileups",
                       dt, self.job.discordant_count,
                       self.job.split_count, len(self._clip_pileups))
            self._check_cancel()

            # Stage 3: Clustering
            t0 = time.time()
            clusters = self._stage_clustering()
            dt = time.time() - t0
            self._report_data["timings"]["clustering"] = round(dt, 1)
            raw_cluster_count = len(clusters)
            self._report_data["raw_clusters_formed"] = raw_cluster_count
            self._log("Clustering: %.1fs — %d raw clusters", dt, raw_cluster_count)
            self._check_cancel()

            # Stage 4: Clip realignment
            t0 = time.time()
            clusters = self._stage_clip_realignment(clusters)
            self._report_data["timings"]["clip_realignment"] = round(time.time() - t0, 1)
            self._check_cancel()

            # Stage 5: Background model
            t0 = time.time()
            clusters = self._stage_background_model(clusters)
            dt = time.time() - t0
            self._report_data["timings"]["background_model"] = round(dt, 1)
            self._log("Background model: %.1fs — chimeric rate %.4f%%",
                       dt, self.job.chimeric_rate * 100)
            self._check_cancel()

            # Stage 6: Filtering (mostly soft flags now)
            t0 = time.time()
            clusters = self._stage_filtering(clusters)
            n_hard_filtered = sum(1 for c in clusters if c.tier == Tier.FILTERED)
            dt = time.time() - t0
            self._report_data["timings"]["filtering"] = round(dt, 1)
            self._log("Filtering: %.1fs — %d hard-filtered, %d retained",
                       dt, n_hard_filtered, len(clusters) - n_hard_filtered)

            # Compute filter breakdown
            fb: dict[str, int] = {}
            for c in clusters:
                for flag in c.filter_flags:
                    fb[flag] = fb.get(flag, 0) + 1
            self._report_data["filter_breakdown"] = fb
            self._check_cancel()

            # Stage 7: Scoring + tier assignment
            t0 = time.time()
            clusters = self._stage_scoring(clusters)
            self._report_data["timings"]["scoring"] = round(time.time() - t0, 1)

            # Stage 8: Breakpoint window aggregation
            t0 = time.time()
            merge_window = self.settings.get("breakpoint_merge_window", 5000)
            clusters = self._breakpoint_window_aggregation(clusters, merge_window)
            self._report_data["timings"]["window_aggregation"] = round(time.time() - t0, 1)

            # Stage 9: Reciprocal consolidation
            clusters = self._dedup_reciprocal(clusters)

            # Re-score after merge/consolidation updated counts
            from scoring import ScoringEngine
            engine = ScoringEngine()
            clusters = engine.score_and_tier(clusters)

            # Stage 10: Top-N selection + near-miss assembly
            top_n = self.settings.get("top_n_candidates", 100)
            all_clusters = clusters  # keep reference for near-miss
            clusters, near_misses = self._select_top_n(clusters, top_n)

            # Tier breakdown
            tier_counts: dict[str, int] = {}
            for tier_name in ["confirmed", "validated", "likely", "candidate", "filtered"]:
                n = sum(1 for c in all_clusters if c.tier.value == tier_name)
                tier_counts[tier_name] = n
                if n:
                    self._log("  %s: %d", tier_name, n)
            self._report_data["tier_counts"] = tier_counts

            # Build report counts
            self._report_data["cluster_counts"] = {
                "raw_clusters_formed": raw_cluster_count,
                "hard_filtered": n_hard_filtered,
                "candidates_retained": len(clusters),
                "confirmed": tier_counts.get("confirmed", 0),
                "validated": tier_counts.get("validated", 0),
                "likely": tier_counts.get("likely", 0),
                "candidate": tier_counts.get("candidate", 0),
                "rejected": tier_counts.get("filtered", 0),
            }
            self._report_data["near_miss_count"] = len(near_misses)

            # Artifact-dominated run sanity gate
            confirmed_count = tier_counts.get("confirmed", 0)
            validated_count = tier_counts.get("validated", 0)
            high_confidence = confirmed_count + validated_count
            if high_confidence > 50:
                self._log(
                    "WARNING: %d high-confidence calls — artifact-dominated run",
                    high_confidence,
                )
                self.emit({
                    "type": "scan.warning",
                    "message": (
                        f"Artifact-dominated run: {high_confidence} high-confidence "
                        f"calls detected. A normal germline should have 0-5."
                    ),
                    "warning_type": "artifact_dominated",
                })
                self._report_data.setdefault("warnings", []).append({
                    "type": "artifact_dominated",
                    "high_confidence_count": high_confidence,
                })

            # Notify scan complete
            self.emit({"type": "scan.completed", "job_id": self.job.job_id})
            self.emit({"type": "validation.started", "job_id": self.job.job_id})

            # Emit validated calls
            for c in clusters:
                self.emit({
                    "type": "validation.call_emitted",
                    "call": c.to_dict(),
                })

            # Stage 11: Output
            t0 = time.time()
            self._stage_output(all_clusters)
            self._report_data["timings"]["output"] = round(time.time() - t0, 1)

            self.emit({"type": "validation.completed", "job_id": self.job.job_id,
                        "num_calls": len(clusters)})

            # Store results: top-N candidates + near-misses
            self.job.validated_calls = [c.to_dict() for c in clusters]
            self.job.settings["_near_misses"] = [c.to_dict() for c in near_misses[:50]]
            self.job.status = JobStatus.COMPLETED
            self.job.stage = ScanStage.COMPLETED
            self.job.completed_at = time.time()

            self.job.settings["_report"] = self._report_data

            elapsed = self.job.completed_at - self.job.started_at
            self._report_data["timings"]["total"] = round(elapsed, 1)
            self._log(
                "Pipeline completed in %.1fs — %d confirmed, %d validated, "
                "%d likely, %d candidate (top %d shown), %d near-misses",
                elapsed,
                tier_counts.get("confirmed", 0),
                tier_counts.get("validated", 0),
                tier_counts.get("likely", 0),
                tier_counts.get("candidate", 0),
                top_n,
                len(near_misses),
            )
            self._log("Results: %s", self.job.results_dir)

        except CancelledException:
            self._log("Pipeline cancelled at stage %s", self.job.stage.value)
            self.job.status = JobStatus.CANCELLED
            self.job.stage = ScanStage.FAILED
            self.job.error = "Cancelled by user"
            self.job.completed_at = time.time()
            self.emit({"type": "scan.cancelled", "job_id": self.job.job_id})

        except Exception as e:
            self._log("Pipeline FAILED: %s", e, level=logging.ERROR)
            self._log("Traceback:\n%s", traceback.format_exc(), level=logging.ERROR)
            self.job.status = JobStatus.FAILED
            self.job.stage = ScanStage.FAILED
            self.job.error = str(e)
            self.emit({"type": "error", "stage": self.job.stage.value,
                        "message": f"Pipeline failed: {e}"})

        finally:
            if self._job_log_handler:
                self._job_log.removeHandler(self._job_log_handler)
                self._job_log_handler.close()

    # ------------------------------------------------------------------
    # Pipeline stages
    # ------------------------------------------------------------------

    def _stage_extraction(self):
        self.job.stage = ScanStage.EXTRACTION
        self.emit({"type": "scan.stage_changed", "stage": "extraction"})

        from extractor import SVExtractor
        extractor = SVExtractor(
            self.job.file_path,
            self.job.reference_path,
            callback=self.emit,
            min_mapq=self.settings.get("min_mapq", 20),
            min_clip_length=self.settings.get("min_clip_length", 20),
            min_split_aligned=self.settings.get("min_split_aligned", 20),
            min_pileup_depth=self.settings.get("min_pileup_depth", 4),
            pileup_window=self.settings.get("pileup_window", 5),
            exclude_chrM=self.settings.get("exclude_chrM", True),
        )
        debug_a = self.settings.get("debug_region_a")
        debug_b = self.settings.get("debug_region_b")
        if debug_a or debug_b:
            extractor.set_debug_regions(
                debug_a, debug_b,
                margin=self.settings.get("debug_margin", 2_000_000),
            )
        if self.settings.get("parallel_extraction", True):
            num_workers = self.settings.get("num_workers", 0) or None
            result = extractor.extract_parallel(
                num_workers=num_workers, cancel_event=self._cancel_event
            )
        else:
            result = extractor.extract()

        self._discordant_reads = result["discordant_reads"]
        self._split_reads = result["split_reads"]
        self._clip_pileups = result["clip_pileups"]
        self._chrom_lengths = result.get("chrom_lengths", {})
        if not self._chrom_lengths:
            chrom_prog = result.get("chrom_progress", {})
            self._chrom_lengths = {k: v.length for k, v in chrom_prog.items() if hasattr(v, 'length')}

        lib = result.get("library_stats", None)
        self.job.insert_size_median = getattr(lib, "median", 0)
        self.job.insert_size_std = getattr(lib, "std", 0)
        self.job.reads_processed = result.get("total_reads_processed", 0)
        self.job.bytes_processed = result.get("total_bytes_processed", 0)
        self.job.discordant_count = len(self._discordant_reads)
        self.job.split_count = len(self._split_reads)
        self.job.clip_count = sum(p.depth for p in self._clip_pileups)
        self.job.total_reads = self.job.reads_processed
        self._debug_evidence = result.get("debug_evidence", [])

    def _stage_clustering(self) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.CLUSTERING
        self.emit({"type": "scan.stage_changed", "stage": "clustering"})

        from clustering import ClusterEngine
        engine = ClusterEngine(
            merge_distance=self.settings.get("merge_distance", 500),
        )
        clusters = engine.cluster(
            self._discordant_reads,
            self._split_reads,
            self._clip_pileups,
            callback=self.emit,
            cancel_event=self._cancel_event,
        )
        return clusters

    def _stage_clip_realignment(self, clusters: list[EvidenceCluster]) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.CLIP_REALIGNMENT
        self.emit({"type": "scan.stage_changed", "stage": "clip_realignment"})

        if self.settings.get("skip_clip_realignment", False):
            return clusters
        if not self._clip_pileups or not self.job.reference_path:
            return clusters

        try:
            clips_with_seqs = [p for p in self._clip_pileups if p.clip_seqs]
            if not clips_with_seqs:
                return clusters

            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.fa', delete=False, prefix='clip_'
            ) as fa:
                fa_path = fa.name
                for i, pileup in enumerate(clips_with_seqs):
                    for j, seq in enumerate(pileup.clip_seqs[:10]):
                        fa.write(f">{pileup.chrom}_{pileup.pos}_{i}_{j}\n{seq}\n")

            sam_path = fa_path.replace('.fa', '.sam')
            cmd = [
                'minimap2', '-a', '-x', 'sr', '--secondary=no',
                self.job.reference_path, fa_path
            ]
            with open(sam_path, 'w') as out:
                subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, timeout=300)

            partner_map = {}
            with open(sam_path) as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 11:
                        continue
                    flag = int(fields[1])
                    if flag & 4:
                        continue
                    qname = fields[0]
                    parts = qname.split('_')
                    if len(parts) >= 3:
                        src_chrom = parts[0]
                        src_pos = int(parts[1])
                        target_chrom = fields[2]
                        target_pos = int(fields[3])
                        mapq = int(fields[4])
                        if mapq >= 10 and target_chrom != src_chrom:
                            partner_map[(src_chrom, src_pos)] = (target_chrom, target_pos)

            for pileup in self._clip_pileups:
                key = (pileup.chrom, pileup.pos)
                if key in partner_map:
                    pileup.partner_chrom, pileup.partner_pos = partner_map[key]

            for cluster in clusters:
                for pileup in self._clip_pileups:
                    if (pileup.partner_chrom and
                        pileup.chrom == cluster.chrom_a and
                        abs(pileup.pos - cluster.pos_a) < 1000 and
                        pileup.partner_chrom == cluster.chrom_b and
                        abs(pileup.partner_pos - cluster.pos_b) < 1000):
                        cluster.clipped_count += pileup.depth
                        cluster.ci_a = (-5, 5)
                        cluster.ci_b = (-5, 5)

            for p in [fa_path, sam_path]:
                try:
                    os.unlink(p)
                except OSError:
                    pass

            self.emit({"type": "scan.progress", "stage": "clip_realignment",
                        "pct": 100,
                        "detail": f"Realigned {len(clips_with_seqs)} clip pileups, "
                                  f"{len(partner_map)} found partners"})

        except FileNotFoundError:
            logger.warning("minimap2 not found, skipping clip realignment")
        except subprocess.TimeoutExpired:
            logger.warning("minimap2 timed out during clip realignment")
        except Exception:
            logger.exception("Clip realignment failed (non-fatal)")

        return clusters

    def _stage_background_model(self, clusters: list[EvidenceCluster]) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.BACKGROUND_MODEL
        self.emit({"type": "scan.stage_changed", "stage": "background_model"})

        from background_model import BackgroundModel
        model = BackgroundModel(
            bin_size=self.settings.get("bg_bin_size", 100_000),
        )
        all_reads = self._discordant_reads + self._split_reads
        model.compute_coverage_bins(all_reads, self._chrom_lengths)
        rate = model.estimate_chimeric_rate(
            self._discordant_reads, self._chrom_lengths,
            total_reads=self.job.reads_processed,
        )
        self.job.chimeric_rate = rate
        model.compute_bin_pair_noise(self._discordant_reads)
        clusters = model.score_clusters(clusters, self._chrom_lengths, callback=self.emit)
        return clusters

    def _stage_filtering(self, clusters: list[EvidenceCluster]) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.FILTERING
        self.emit({"type": "scan.stage_changed", "stage": "filtering"})

        blacklist_index = None
        try:
            from regions import load_encode_blacklist_v2
            blacklist_index = load_encode_blacklist_v2()
        except Exception:
            logger.warning("Failed to load ENCODE blacklist", exc_info=True)

        from filters import FilterEngine
        engine = FilterEngine(
            centromere_margin=self.settings.get("centromere_margin", 1_000_000),
            bg_pvalue_threshold=self.settings.get("bg_pvalue_threshold", 0.001),
            blacklist_index=blacklist_index,
        )
        clusters = engine.apply_filters(clusters, self._chrom_lengths, callback=self.emit)
        return clusters

    def _stage_scoring(self, clusters: list[EvidenceCluster]) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.SCORING
        self.emit({"type": "scan.stage_changed", "stage": "scoring"})

        from scoring import ScoringEngine
        engine = ScoringEngine()
        clusters = engine.score_and_tier(clusters, callback=self.emit)
        return clusters

    # ------------------------------------------------------------------
    # Breakpoint window aggregation (merge nearby clusters)
    # ------------------------------------------------------------------

    def _breakpoint_window_aggregation(
        self,
        clusters: list[EvidenceCluster],
        merge_window: int = 5000,
    ) -> list[EvidenceCluster]:
        """Merge clusters within a breakpoint window on same chrom pair + orientation.

        A true translocation may be split into multiple nearby clusters
        due to alignment noise.  This step merges them while preserving
        the subcluster IDs for transparency.
        """
        # Separate filtered from non-filtered
        passing = [c for c in clusters if c.tier != Tier.FILTERED]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        # Group by (chrom_a, chrom_b, orientation)
        groups: dict[tuple[str, str, str], list[EvidenceCluster]] = defaultdict(list)
        for c in passing:
            groups[(c.chrom_a, c.chrom_b, c.orientation)].append(c)

        merged_all: list[EvidenceCluster] = []

        for key, group in groups.items():
            group.sort(key=lambda c: c.pos_a)
            merged: list[EvidenceCluster] = [group[0]]

            for c in group[1:]:
                last = merged[-1]
                if (abs(c.pos_a - last.pos_a) <= merge_window
                        and abs(c.pos_b - last.pos_b) <= merge_window):
                    # Merge c into last
                    last.discordant_count += c.discordant_count
                    last.split_count += c.split_count
                    last.clipped_count += c.clipped_count
                    last.reads.extend(c.reads)
                    last.merged_subclusters.append(
                        c.cluster_id or c.event_id
                    )
                    # Update unique starts
                    if last.reads:
                        last.unique_starts_a = len(set(
                            r.pos_a for r in last.reads
                        ))
                        last.unique_starts_b = len(set(
                            r.pos_b for r in last.reads
                        ))
                    # Take better MAPQ
                    last.median_mapq = max(last.median_mapq, c.median_mapq)
                    # Take better background_p
                    last.background_p = min(last.background_p, c.background_p)
                    # Merge filter flags (deduplicated)
                    for flag in c.filter_flags:
                        if flag not in last.filter_flags:
                            last.filter_flags.append(flag)
                else:
                    merged.append(c)
            merged_all.extend(merged)

        n_before = len(passing)
        n_after = len(merged_all)
        if n_before != n_after:
            self._log(
                "Window aggregation (±%d bp): %d → %d clusters",
                merge_window, n_before, n_after,
            )

        return merged_all + filtered

    # ------------------------------------------------------------------
    # Reciprocal consolidation
    # ------------------------------------------------------------------

    @staticmethod
    def _dedup_reciprocal(clusters: list[EvidenceCluster]) -> list[EvidenceCluster]:
        """Merge reciprocal cluster pairs (A→B ++ with A→B --) into single calls."""
        _reciprocal_orient = {
            "++": "--", "--": "++", "+-": "-+", "-+": "+-",
        }
        compat_distance = 1000

        passing = [c for c in clusters if c.tier != Tier.FILTERED]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        index: dict[tuple[str, str], list[EvidenceCluster]] = defaultdict(list)
        for c in passing:
            index[(c.chrom_a, c.chrom_b)].append(c)

        absorbed: set[str] = set()

        for c in passing:
            if c.event_id in absorbed:
                continue

            expected_orient = _reciprocal_orient.get(c.orientation)
            if expected_orient is None:
                continue

            candidates = index.get((c.chrom_a, c.chrom_b), [])
            for other in candidates:
                if other is c or other.event_id in absorbed:
                    continue
                if other.orientation != expected_orient:
                    continue
                if (abs(other.pos_a - c.pos_a) <= compat_distance
                        and abs(other.pos_b - c.pos_b) <= compat_distance):
                    if c.score >= other.score:
                        keeper, donor = c, other
                    else:
                        keeper, donor = other, c

                    keeper.reciprocal_support = donor.total_support
                    keeper.discordant_count += donor.discordant_count
                    keeper.split_count += donor.split_count
                    keeper.clipped_count += donor.clipped_count
                    keeper.reads.extend(donor.reads)
                    keeper.merged_subclusters.append(
                        donor.cluster_id or donor.event_id
                    )
                    # Update unique starts
                    if keeper.reads:
                        keeper.unique_starts_a = len(set(
                            r.pos_a for r in keeper.reads
                        ))
                        keeper.unique_starts_b = len(set(
                            r.pos_b for r in keeper.reads
                        ))
                    absorbed.add(donor.event_id)
                    break

        deduped = [c for c in passing if c.event_id not in absorbed]
        deduped.sort(key=lambda c: c.score, reverse=True)
        for idx, c in enumerate(deduped, start=1):
            c.cluster_id = f"CLU_{idx:03d}"

        logger.info(
            "Reciprocal dedup: %d → %d non-filtered (%d absorbed)",
            len(passing), len(deduped), len(absorbed),
        )
        return deduped + filtered

    # ------------------------------------------------------------------
    # Top-N selection + near-miss assembly
    # ------------------------------------------------------------------

    @staticmethod
    def _select_top_n(
        clusters: list[EvidenceCluster],
        top_n: int = 100,
    ) -> tuple[list[EvidenceCluster], list[EvidenceCluster]]:
        """Select top-N candidates for display + assemble near-misses.

        All confirmed/validated/likely calls are always included.
        Then the top-N candidates by score fill the remainder.
        Near-misses are candidates that didn't make the cut.

        Returns:
            (selected, near_misses)
        """
        # Always include high-tier calls
        high_tier = {Tier.CONFIRMED, Tier.VALIDATED, Tier.LIKELY}
        always_show = [c for c in clusters if c.tier in high_tier]
        candidates = [c for c in clusters if c.tier == Tier.CANDIDATE]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        # Candidates already sorted by score (from scoring step)
        candidates.sort(key=lambda c: c.score, reverse=True)

        remaining_slots = max(0, top_n - len(always_show))
        selected_candidates = candidates[:remaining_slots]
        near_miss_candidates = candidates[remaining_slots:]

        # Near-misses also include top filtered clusters with
        # total_support >= 3 (they were close but artifact-flagged)
        interesting_filtered = [
            c for c in filtered if c.total_support >= 3
        ]
        interesting_filtered.sort(key=lambda c: c.total_support, reverse=True)

        selected = always_show + selected_candidates
        selected.sort(key=lambda c: c.score, reverse=True)

        # Re-assign sequential IDs
        for idx, c in enumerate(selected, start=1):
            c.cluster_id = f"CLU_{idx:03d}"

        near_misses = near_miss_candidates + interesting_filtered[:20]
        near_misses.sort(key=lambda c: c.score, reverse=True)

        return selected, near_misses

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _stage_output(self, clusters: list[EvidenceCluster]):
        self.job.stage = ScanStage.OUTPUT
        self.emit({"type": "scan.stage_changed", "stage": "output"})

        results_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            '..', 'results', self.job.job_id
        )
        os.makedirs(results_dir, exist_ok=True)
        self.job.results_dir = results_dir

        from output_writer import OutputWriter
        writer = OutputWriter(results_dir, self.job.reference_build)
        files = writer.write_all(self.job, clusters, self._chrom_lengths)

        # Write debug evidence TSV if collected
        if getattr(self, '_debug_evidence', None):
            debug_path = os.path.join(results_dir, 'debug_evidence.tsv')
            _COLS = [
                "read_hash", "chrom_a", "pos_a", "chrom_b", "pos_b",
                "mapq", "sa_mapq", "evidence_type", "cigar", "orientation",
                "accepted", "rejection_reason",
            ]
            with open(debug_path, 'w') as f:
                f.write('\t'.join(_COLS) + '\n')
                for row in self._debug_evidence:
                    vals = [str(row.get(c, '')) for c in _COLS]
                    f.write('\t'.join(vals) + '\n')
            files['debug_evidence'] = debug_path
            self._log(
                "Debug evidence: %d rows to %s",
                len(self._debug_evidence), debug_path,
            )
