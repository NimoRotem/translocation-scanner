"""Pipeline v2 orchestrator — dual-track architecture.

Track 1: Blind discovery (custom) — primary
Track 2: External caller validation (Manta/DELLY/GRIDSS) — corroboration

Stages:
  1. SV-read extraction (parallel, single-pass)
  2. Sparse 2D bin map per chrom pair (1kb bins)
  3. Cluster promotion + breakpoint refinement
  4. Per-cluster annotation (all metrics from spec)
  5. External caller validation (Track 2)
  6. Local NB background model
  7. Filtering (hard excludes + reject)
  8. Tier assignment
  9. Output
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import threading
import time
import traceback
from collections import defaultdict
from typing import Optional

import numpy as np

from models import (
    ScanJob, JobStatus, ScanStage, EvidenceCluster, Tier
)

logger = logging.getLogger(__name__)


class CancelledException(Exception):
    pass


class PipelineV2:
    """Dual-track translocation detection pipeline."""

    def __init__(self, job: ScanJob, event_callback=None, cancel_event=None):
        self.job = job
        self.emit = event_callback or (lambda e: None)
        self._cancel_event = cancel_event
        self.settings = job.settings or {}
        self._chrom_lengths: dict[str, int] = {}
        self._report_data: dict = {"timings": {}, "filter_breakdown": {}}
        self._debug_evidence: list[dict] = []
        self._job_log: Optional[logging.Logger] = None
        self._job_log_handler: Optional[logging.Handler] = None
        self._masks = None

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
        """Execute the full dual-track pipeline."""
        self._setup_job_log()
        try:
            self.job.status = JobStatus.RUNNING
            self.job.started_at = time.time()
            self._log("Pipeline v2 started for job %s", self.job.job_id)
            self._log("Input: %s", self.job.file_path)
            self._log("Settings: %s", self.settings)
            self.emit({"type": "scan.started", "job_id": self.job.job_id})

            # Load masks (from cache, never fetch at runtime)
            self._load_masks()

            # ==========================================
            # TRACK 1: BLIND DISCOVERY
            # ==========================================

            # Stage 1: Extract SV reads
            t0 = time.time()
            self._stage_extraction()
            dt = time.time() - t0
            self._report_data["timings"]["extraction"] = round(dt, 1)
            self._log("Extraction: %.1fs — %d disc, %d split",
                       dt, self.job.discordant_count, self.job.split_count)
            self._check_cancel()

            # Stage 2: Clustering (bin map + promotion + refinement)
            t0 = time.time()
            clusters = self._stage_clustering()
            dt = time.time() - t0
            self._report_data["timings"]["clustering"] = round(dt, 1)
            self._report_data.setdefault("cluster_counts", {})["raw_clusters_formed"] = len(clusters)
            self._log("Clustering: %.1fs — %d clusters", dt, len(clusters))
            self._check_cancel()

            # Stage 3: Clip realignment
            t0 = time.time()
            clusters = self._stage_clip_realignment(clusters)
            self._report_data["timings"]["clip_realignment"] = round(time.time() - t0, 1)
            self._check_cancel()

            # Stage 4: Per-cluster annotation
            t0 = time.time()
            self._annotate_clusters(clusters)
            self._report_data["timings"]["annotation"] = round(time.time() - t0, 1)
            self._check_cancel()

            # ==========================================
            # TRACK 2: EXTERNAL CALLER VALIDATION
            # ==========================================
            t0 = time.time()
            external_bnds = self._stage_external_callers()
            self._report_data["timings"]["external_callers"] = round(time.time() - t0, 1)
            self._match_external_calls(clusters, external_bnds)
            self._check_cancel()

            # ==========================================
            # SCORING & FILTERING
            # ==========================================

            # Stage 5: Local NB background model
            t0 = time.time()
            clusters = self._stage_background_model(clusters)
            dt = time.time() - t0
            self._report_data["timings"]["background_model"] = round(dt, 1)
            self._check_cancel()

            # Stage 6: Hard excludes + reject + tier assignment
            t0 = time.time()
            clusters = self._stage_filtering_and_tiers(clusters)
            self._report_data["timings"]["filtering"] = round(time.time() - t0, 1)
            self._report_data.setdefault("cluster_counts", {})["candidates_retained"] = sum(1 for c in clusters if c.tier != Tier.FILTERED)
            # Build filter breakdown from filter_flags
            fb: dict[str, int] = defaultdict(int)
            for c in clusters:
                for flag in c.filter_flags:
                    fb[flag] += 1
            self._report_data["filter_breakdown"] = dict(fb)
            self._check_cancel()

            # Stage 7: Breakpoint window aggregation
            t0 = time.time()
            merge_window = self.settings.get("breakpoint_merge_window", 5000)
            clusters = self._breakpoint_window_aggregation(clusters, merge_window)
            self._report_data["timings"]["window_aggregation"] = round(time.time() - t0, 1)

            # Reciprocal consolidation
            clusters = self._dedup_reciprocal(clusters)

            # Re-run tier assignment after merges
            from filters_v2 import FilterEngineV2
            fe = FilterEngineV2(masks=self._masks)
            for c in clusters:
                if c.tier != Tier.FILTERED:
                    fe._assign_tier(c)

            # Score clusters for within-tier ranking
            self._compute_scores(clusters)

            # Top-N selection
            top_n = self.settings.get("top_n_candidates", 100)
            all_clusters = clusters
            clusters, near_misses = self._select_top_n(clusters, top_n)

            # Tier breakdown
            tier_counts = {}
            for t_name in ["confirmed", "validated", "likely", "strong_candidate", "candidate", "filtered"]:
                n = sum(1 for c in all_clusters if c.tier.value == t_name)
                tier_counts[t_name] = n
                if n:
                    self._log("  %s: %d", t_name, n)
            self._report_data["tier_counts"] = tier_counts

            # Emit results
            self.emit({"type": "scan.completed", "job_id": self.job.job_id})
            self.emit({"type": "validation.started", "job_id": self.job.job_id})
            for c in clusters:
                self.emit({"type": "validation.call_emitted", "call": c.to_dict()})

            # Stage 8: Output
            t0 = time.time()
            self._stage_output(all_clusters)
            self._report_data["timings"]["output"] = round(time.time() - t0, 1)

            self.emit({"type": "validation.completed", "job_id": self.job.job_id,
                        "num_calls": len(clusters)})

            self.job.validated_calls = [c.to_dict() for c in clusters]
            self.job.settings["_near_misses"] = [c.to_dict() for c in near_misses[:50]]
            self.job.status = JobStatus.COMPLETED
            self.job.stage = ScanStage.COMPLETED
            self.job.completed_at = time.time()

            # Store mask manifest version in report
            if self._masks:
                self._report_data["mask_manifest_version"] = self._masks.manifest_version
            self.job.settings["_report"] = self._report_data

            elapsed = self.job.completed_at - self.job.started_at
            self._report_data["timings"]["total"] = round(elapsed, 1)
            self._log("Pipeline completed in %.1fs — %s", elapsed, tier_counts)

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
    # Mask loading
    # ------------------------------------------------------------------

    def _load_masks(self):
        """Load mask tracks from cache (never fetch at runtime)."""
        mask_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'masks', 'data'
        )
        if os.path.isdir(mask_dir):
            try:
                from masks import MaskSet
                self._masks = MaskSet.load(mask_dir)
                self._log("Loaded masks from %s (version: %s)",
                           mask_dir, self._masks.manifest_version)
            except Exception as e:
                self._log("Failed to load masks: %s", e, level=logging.WARNING)
        else:
            self._log("No mask directory at %s — run download_masks.py first",
                       mask_dir, level=logging.WARNING)

    # ------------------------------------------------------------------
    # Stage 1: Extraction
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
            min_clip_length=self.settings.get("min_clip_length", 15),
            min_split_aligned=self.settings.get("min_split_aligned", 30),
            min_pileup_depth=self.settings.get("min_pileup_depth", 4),
            pileup_window=self.settings.get("pileup_window", 5),
            exclude_chrM=True,
        )
        debug_a = self.settings.get("debug_region_a")
        debug_b = self.settings.get("debug_region_b")
        if debug_a or debug_b:
            extractor.set_debug_regions(debug_a, debug_b,
                margin=self.settings.get("debug_margin", 2_000_000))

        if self.settings.get("parallel_extraction", True):
            result = extractor.extract_parallel(
                num_workers=self.settings.get("num_workers") or None,
                cancel_event=self._cancel_event,
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

    # ------------------------------------------------------------------
    # Stage 2: Clustering
    # ------------------------------------------------------------------

    def _stage_clustering(self) -> list[EvidenceCluster]:
        self.job.stage = ScanStage.CLUSTERING
        self.emit({"type": "scan.stage_changed", "stage": "clustering"})

        from clustering import ClusterEngine
        engine = ClusterEngine(
            merge_distance=self.settings.get("merge_distance", 500),
            min_cluster_support=self.settings.get("min_cluster_support", 3),
        )
        return engine.cluster(
            self._discordant_reads,
            self._split_reads,
            self._clip_pileups,
            callback=self.emit,
            cancel_event=self._cancel_event,
        )

    # ------------------------------------------------------------------
    # Stage 3: Clip realignment (same as v1)
    # ------------------------------------------------------------------

    def _stage_clip_realignment(self, clusters):
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

            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False, prefix='clip_') as fa:
                fa_path = fa.name
                for i, pileup in enumerate(clips_with_seqs):
                    for j, seq in enumerate(pileup.clip_seqs[:10]):
                        fa.write(f">{pileup.chrom}_{pileup.pos}_{i}_{j}\n{seq}\n")

            sam_path = fa_path.replace('.fa', '.sam')
            cmd = ['minimap2', '-a', '-x', 'sr', '--secondary=no', self.job.reference_path, fa_path]
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
                    parts = fields[0].split('_')
                    if len(parts) >= 3:
                        src_chrom, src_pos = parts[0], int(parts[1])
                        target_chrom, target_pos = fields[2], int(fields[3])
                        mapq = int(fields[4])
                        if mapq >= 10 and target_chrom != src_chrom:
                            partner_map[(src_chrom, src_pos)] = (target_chrom, target_pos)

            for pileup in self._clip_pileups:
                key = (pileup.chrom, pileup.pos)
                if key in partner_map:
                    pileup.partner_chrom, pileup.partner_pos = partner_map[key]

            for cluster in clusters:
                for pileup in self._clip_pileups:
                    if (pileup.partner_chrom
                        and pileup.chrom == cluster.chrom_a
                        and abs(pileup.pos - cluster.pos_a) < 1000
                        and pileup.partner_chrom == cluster.chrom_b
                        and abs(pileup.partner_pos - cluster.pos_b) < 1000):
                        cluster.clipped_count += pileup.depth
                        cluster.ci_a = (-5, 5)
                        cluster.ci_b = (-5, 5)

            for p in [fa_path, sam_path]:
                try:
                    os.unlink(p)
                except OSError:
                    pass

        except FileNotFoundError:
            self._log("minimap2 not found, skipping clip realignment", level=logging.WARNING)
        except subprocess.TimeoutExpired:
            self._log("minimap2 timed out", level=logging.WARNING)
        except Exception:
            self._log("Clip realignment failed (non-fatal)", level=logging.WARNING)

        return clusters

    # ------------------------------------------------------------------
    # Stage 4: Per-cluster annotation
    # ------------------------------------------------------------------

    def _annotate_clusters(self, clusters: list[EvidenceCluster]) -> None:
        """Compute all per-cluster annotations from spec §1.5."""
        self.emit({"type": "scan.stage_changed", "stage": "annotation"})

        # Build spatial index for O(1) local discordant count lookups
        # Key: (chrom, 10kb_bin) -> count of discordant reads in that bin
        disc_bin_index: dict[tuple[str, int], int] = defaultdict(int)
        pair_counts: dict[tuple[str, str], int] = defaultdict(int)
        chrom_bin_counts: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        for r in self._discordant_reads:
            # Index read position
            b = r.pos // 10000
            disc_bin_index[(r.chrom, b)] += 1
            chrom_bin_counts[r.chrom][b] += 1
            # Index mate position
            if r.mate_chrom and r.mate_pos is not None:
                mb = r.mate_pos // 10000
                disc_bin_index[(r.mate_chrom, mb)] += 1
                if r.mate_chrom != r.chrom:
                    chrom_bin_counts[r.mate_chrom][mb] += 1
            # Pair counts for enrichment
            if r.mate_chrom and r.mate_chrom != r.chrom:
                a, b_chr = (r.chrom, r.mate_chrom) if r.chrom <= r.mate_chrom else (r.mate_chrom, r.chrom)
                pair_counts[(a, b_chr)] += 1

        total_interchrom = sum(pair_counts.values())

        # Precompute per-chrom median coverage from bin counts
        chrom_median_cov: dict[str, float] = {}
        for chrom, bin_dict in chrom_bin_counts.items():
            if bin_dict:
                counts = list(bin_dict.values())
                chrom_median_cov[chrom] = float(np.median(counts)) if counts else 1.0
            else:
                chrom_median_cov[chrom] = 1.0

        self._log("Built annotation index: %d bins, %d chrom pair types", len(disc_bin_index), len(pair_counts))

        # Promiscuous loci: positions appearing in >N clusters
        side_a_counts: dict[tuple[str, int], int] = defaultdict(int)
        side_b_counts: dict[tuple[str, int], int] = defaultdict(int)
        for c in clusters:
            bin_a = c.pos_a // 10000
            bin_b = c.pos_b // 10000
            side_a_counts[(c.chrom_a, bin_a)] += 1
            side_b_counts[(c.chrom_b, bin_b)] += 1

        promiscuous_threshold = self.settings.get("promiscuous_threshold", 5)

        for c in clusters:
            # Duplicate fraction
            if c.reads:
                dup_count = sum(1 for r in c.reads if r.flag & 0x400)
                c.duplicate_fraction = dup_count / len(c.reads) if c.reads else 0.0
            else:
                c.duplicate_fraction = 0.0

            # Per-side median MAPQ
            if c.reads:
                mapqs_a = [r.mapq for r in c.reads]
                c.median_mapq_a = float(np.median(mapqs_a)) if mapqs_a else 0.0
                c.median_mapq_b = c.median_mapq_a  # approximation from available data

            # Local coverage ratio — use indexed bin lookup
            med_a = chrom_median_cov.get(c.chrom_a, 1.0)
            med_b = chrom_median_cov.get(c.chrom_b, 1.0)
            rate_a = float(disc_bin_index.get((c.chrom_a, c.pos_a // 10000), 0))
            rate_b = float(disc_bin_index.get((c.chrom_b, c.pos_b // 10000), 0))
            c.local_coverage_ratio_a = rate_a / max(med_a, 1.0)
            c.local_coverage_ratio_b = rate_b / max(med_b, 1.0)

            # Chromosome pair enrichment
            pair_key = (c.chrom_a, c.chrom_b) if c.chrom_a <= c.chrom_b else (c.chrom_b, c.chrom_a)
            observed_pair = pair_counts.get(pair_key, 0)
            if total_interchrom > 0:
                n_pairs = len(self._chrom_lengths) * (len(self._chrom_lengths) - 1) / 2
                expected = total_interchrom / max(n_pairs, 1)
                c.chrom_pair_enrichment = observed_pair / max(expected, 1.0)

            # Orientation distribution
            if c.reads:
                orient_dist: dict[str, int] = defaultdict(int)
                for r in c.reads:
                    sa = "-" if r.is_reverse else "+"
                    sb = "-" if r.mate_is_reverse else "+"
                    orient_dist[sa + sb] += 1
                c.orientation_distribution = dict(orient_dist)

            # Promiscuous hotspot
            bin_a = c.pos_a // 10000
            bin_b = c.pos_b // 10000
            if (side_a_counts.get((c.chrom_a, bin_a), 0) > promiscuous_threshold
                    or side_b_counts.get((c.chrom_b, bin_b), 0) > promiscuous_threshold):
                c.promiscuous_hotspot = True

        self._log("Annotated %d clusters", len(clusters))

    # ------------------------------------------------------------------
    # Track 2: External callers
    # ------------------------------------------------------------------

    def _stage_external_callers(self) -> list[dict]:
        """Run Manta/DELLY/GRIDSS and collect interchromosomal BNDs."""
        self.emit({"type": "scan.stage_changed", "stage": "external_callers"})

        if self.settings.get("skip_external_callers", False):
            self._log("External callers skipped by settings")
            return []

        bnds = []

        # Try Manta
        manta_bnds = self._run_manta()
        bnds.extend(manta_bnds)

        # Try DELLY
        delly_bnds = self._run_delly()
        bnds.extend(delly_bnds)

        self._log("External callers: %d total BNDs (%d Manta, %d DELLY)",
                   len(bnds), len(manta_bnds), len(delly_bnds))
        return bnds

    def _run_manta(self) -> list[dict]:
        """Run Manta and parse interchromosomal BNDs."""
        if not shutil.which("configManta.py"):
            self._log("Manta not installed, skipping", level=logging.WARNING)
            return []

        try:
            workdir = tempfile.mkdtemp(prefix="manta_")
            cmd_config = [
                "configManta.py",
                "--bam", self.job.file_path,
                "--referenceFasta", self.job.reference_path,
                "--runDir", workdir,
            ]
            subprocess.run(cmd_config, check=True, capture_output=True, timeout=120)

            cmd_run = [os.path.join(workdir, "runWorkflow.py"), "-j", "4"]
            subprocess.run(cmd_run, check=True, capture_output=True, timeout=7200)

            # Parse diploidSV VCF
            vcf_path = os.path.join(workdir, "results", "variants", "diploidSV.vcf.gz")
            if not os.path.exists(vcf_path):
                vcf_path = os.path.join(workdir, "results", "variants", "diploidSV.vcf")

            bnds = self._parse_bnd_vcf(vcf_path, "manta")
            shutil.rmtree(workdir, ignore_errors=True)
            self._log("Manta: %d interchromosomal BNDs", len(bnds))
            return bnds

        except FileNotFoundError:
            return []
        except subprocess.TimeoutExpired:
            self._log("Manta timed out", level=logging.WARNING)
            return []
        except Exception as e:
            self._log("Manta failed: %s", e, level=logging.WARNING)
            return []

    def _run_delly(self) -> list[dict]:
        """Run DELLY and parse interchromosomal BNDs."""
        if not shutil.which("delly"):
            self._log("DELLY not installed, skipping", level=logging.WARNING)
            return []

        try:
            with tempfile.NamedTemporaryFile(suffix=".bcf", delete=False) as tmp:
                out_path = tmp.name

            # Build exclude file for non-reference chromosomes
            exclude_path = None
            ref_chroms = set(self._chrom_lengths.keys())
            if ref_chroms:
                import pysam
                with pysam.AlignmentFile(self.job.file_path, "rb") as af:
                    bam_chroms = set(af.references)
                non_ref = bam_chroms - ref_chroms
                if non_ref:
                    exclude_path = out_path + ".exclude.tsv"
                    with open(exclude_path, "w") as ef:
                        for chrom in sorted(non_ref):
                            ef.write(f"{chrom}\n")
                    self._log("DELLY exclude: %d non-reference chroms", len(non_ref))

            cmd = [
                "delly", "call",
                "-g", self.job.reference_path,
                "-o", out_path,
            ]
            if exclude_path:
                cmd.extend(["-x", exclude_path])
            cmd.append(self.job.file_path)
            subprocess.run(cmd, check=True, capture_output=True, timeout=7200)

            # Convert BCF to VCF for parsing
            vcf_path = out_path.replace(".bcf", ".vcf")
            subprocess.run(
                ["bcftools", "view", "-o", vcf_path, out_path],
                check=True, capture_output=True, timeout=120,
            )

            bnds = self._parse_bnd_vcf(vcf_path, "delly")
            for p in [out_path, vcf_path, exclude_path]:
                if p:
                    try:
                        os.unlink(p)
                    except OSError:
                        pass

            self._log("DELLY: %d interchromosomal BNDs", len(bnds))
            return bnds

        except FileNotFoundError:
            return []
        except subprocess.TimeoutExpired:
            self._log("DELLY timed out", level=logging.WARNING)
            return []
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.decode("utf-8", errors="replace")[:500] if e.stderr else ""
            self._log("DELLY failed (exit %s): %s", e.returncode, stderr, level=logging.WARNING)
            return []
        except Exception as e:
            self._log("DELLY failed: %s", e, level=logging.WARNING)
            return []

    @staticmethod
    def _parse_bnd_vcf(vcf_path: str, caller: str) -> list[dict]:
        """Parse BND records from a VCF, keep only interchromosomal."""
        bnds = []
        if not os.path.exists(vcf_path):
            return bnds

        import gzip
        opener = gzip.open if vcf_path.endswith(".gz") else open
        with opener(vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue
                info = fields[7]
                if "SVTYPE=BND" not in info:
                    continue
                chrom = fields[0]
                pos = int(fields[1])
                alt = fields[4]
                # Parse BND alt to get mate chrom/pos
                mate_chrom, mate_pos = _parse_bnd_alt(alt)
                if mate_chrom and mate_chrom != chrom:
                    bnds.append({
                        "caller": caller,
                        "chrom_a": chrom,
                        "pos_a": pos,
                        "chrom_b": mate_chrom,
                        "pos_b": mate_pos,
                        "orientation": _bnd_orientation(alt),
                        "info": info,
                    })
        return bnds

    def _match_external_calls(
        self, clusters: list[EvidenceCluster], external_bnds: list[dict],
        window: int = 1000,
    ) -> None:
        """Match external caller BNDs to Track 1 clusters."""
        if not external_bnds:
            return

        for cluster in clusters:
            for bnd in external_bnds:
                # Check both orientations of the match
                match_ab = (
                    cluster.chrom_a == bnd["chrom_a"]
                    and cluster.chrom_b == bnd["chrom_b"]
                    and abs(cluster.pos_a - bnd["pos_a"]) <= window
                    and abs(cluster.pos_b - bnd["pos_b"]) <= window
                )
                match_ba = (
                    cluster.chrom_a == bnd["chrom_b"]
                    and cluster.chrom_b == bnd["chrom_a"]
                    and abs(cluster.pos_a - bnd["pos_b"]) <= window
                    and abs(cluster.pos_b - bnd["pos_a"]) <= window
                )
                if match_ab or match_ba:
                    caller = bnd["caller"]
                    if caller not in cluster.external_callers:
                        cluster.external_callers.append(caller)

        matched = sum(1 for c in clusters if c.external_callers)
        self._log("External caller matching: %d/%d clusters have external support",
                   matched, len(clusters))

    # ------------------------------------------------------------------
    # Stage 5: Background model
    # ------------------------------------------------------------------

    def _stage_background_model(self, clusters):
        self.job.stage = ScanStage.BACKGROUND_MODEL
        self.emit({"type": "scan.stage_changed", "stage": "background_model"})

        from background_model_v2 import BackgroundModelV2
        model = BackgroundModelV2(
            window_size=self.settings.get("bg_window_size", 10_000),
            bin_size=self.settings.get("bg_bin_size", 1_000),
        )
        model.compute_local_rates(
            self._discordant_reads, self._chrom_lengths,
            total_reads=self.job.reads_processed,
        )
        self.job.chimeric_rate = model.chimeric_rate
        clusters = model.score_clusters(clusters, self._chrom_lengths, callback=self.emit)
        return clusters

    # ------------------------------------------------------------------
    # Stage 6: Filtering + tier assignment
    # ------------------------------------------------------------------

    def _stage_filtering_and_tiers(self, clusters):
        self.job.stage = ScanStage.FILTERING
        self.emit({"type": "scan.stage_changed", "stage": "filtering"})

        from filters_v2 import FilterEngineV2
        engine = FilterEngineV2(masks=self._masks)
        return engine.apply_all(clusters, self._chrom_lengths, callback=self.emit)

    # ------------------------------------------------------------------
    # Window aggregation + reciprocal dedup (same as v1)
    # ------------------------------------------------------------------

    def _breakpoint_window_aggregation(self, clusters, merge_window=5000):
        passing = [c for c in clusters if c.tier != Tier.FILTERED]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        groups: dict[tuple, list] = defaultdict(list)
        for c in passing:
            groups[(c.chrom_a, c.chrom_b, c.orientation)].append(c)

        merged_all = []
        for key, group in groups.items():
            group.sort(key=lambda c: c.pos_a)
            merged = [group[0]]
            for c in group[1:]:
                last = merged[-1]
                if abs(c.pos_a - last.pos_a) <= merge_window and abs(c.pos_b - last.pos_b) <= merge_window:
                    last.discordant_count += c.discordant_count
                    last.split_count += c.split_count
                    last.clipped_count += c.clipped_count
                    last.reads.extend(c.reads)
                    last.merged_subclusters.append(c.cluster_id or c.event_id)
                    if last.reads:
                        last.unique_starts_a = len(set(r.pos_a for r in last.reads))
                        last.unique_starts_b = len(set(r.pos_b for r in last.reads))
                    last.median_mapq = max(last.median_mapq, c.median_mapq)
                    last.background_p = min(last.background_p, c.background_p)
                    # Merge external callers
                    for ec in c.external_callers:
                        if ec not in last.external_callers:
                            last.external_callers.append(ec)
                else:
                    merged.append(c)
            merged_all.extend(merged)

        return merged_all + filtered

    @staticmethod
    def _dedup_reciprocal(clusters):
        _recip = {"++": "--", "--": "++", "+-": "-+", "-+": "+-"}
        passing = [c for c in clusters if c.tier != Tier.FILTERED]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        index: dict[tuple, list] = defaultdict(list)
        for c in passing:
            index[(c.chrom_a, c.chrom_b)].append(c)

        absorbed = set()
        for c in passing:
            if c.event_id in absorbed:
                continue
            exp = _recip.get(c.orientation)
            if not exp:
                continue
            for other in index.get((c.chrom_a, c.chrom_b), []):
                if other is c or other.event_id in absorbed:
                    continue
                if other.orientation != exp:
                    continue
                if abs(other.pos_a - c.pos_a) <= 1000 and abs(other.pos_b - c.pos_b) <= 1000:
                    keeper = c if c.score >= other.score else other
                    donor = other if keeper is c else c
                    keeper.reciprocal_support = donor.total_support
                    keeper.discordant_count += donor.discordant_count
                    keeper.split_count += donor.split_count
                    keeper.clipped_count += donor.clipped_count
                    keeper.reads.extend(donor.reads)
                    keeper.merged_subclusters.append(donor.cluster_id or donor.event_id)
                    for ec in donor.external_callers:
                        if ec not in keeper.external_callers:
                            keeper.external_callers.append(ec)
                    absorbed.add(donor.event_id)
                    break

        deduped = [c for c in passing if c.event_id not in absorbed]
        deduped.sort(key=lambda c: c.score, reverse=True)
        for idx, c in enumerate(deduped, 1):
            c.cluster_id = f"CLU_{idx:03d}"
        return deduped + filtered

    @staticmethod
    def _select_top_n(clusters, top_n=100):
        high_tier = {Tier.CONFIRMED, Tier.VALIDATED, Tier.LIKELY, Tier.STRONG_CANDIDATE}
        always_show = [c for c in clusters if c.tier in high_tier]
        candidates = [c for c in clusters if c.tier == Tier.CANDIDATE]
        filtered = [c for c in clusters if c.tier == Tier.FILTERED]

        candidates.sort(key=lambda c: c.score, reverse=True)
        remaining = max(0, top_n - len(always_show))
        selected = always_show + candidates[:remaining]
        near_misses = candidates[remaining:]

        interesting_filtered = [c for c in filtered if c.total_support >= 3]
        interesting_filtered.sort(key=lambda c: c.total_support, reverse=True)

        selected.sort(key=lambda c: c.score, reverse=True)
        for idx, c in enumerate(selected, 1):
            c.cluster_id = f"CLU_{idx:03d}"

        return selected, near_misses + interesting_filtered[:20]

    # ------------------------------------------------------------------
    # Scoring: compute composite score for within-tier ranking
    # ------------------------------------------------------------------

    def _compute_scores(self, clusters: list[EvidenceCluster]) -> None:
        """Compute a composite score for each cluster.

        Score components (all non-negative):
          - support_pr: log2(PR+1) * 5              (discordant pair count)
          - support_sr: log2(SR+1) * 10             (split read count, 2x weight)
          - support_clip: log2(clip+1) * 2
          - pvalue_a: -log10(pval_a) capped at 20   (statistical significance side A)
          - pvalue_b: -log10(pval_b) capped at 20   (statistical significance side B)
          - unique_starts: min(unique_starts_a, unique_starts_b) * 0.5
          - reciprocal: 10 if reciprocal_support > 0 else 0
          - external: 15 * len(external_callers)
          - mapq: median_mapq * 0.2
        """
        import math
        for c in clusters:
            if c.tier == Tier.FILTERED:
                c.score = 0.0
                continue

            comps = {}
            comps["support_pr"] = round(math.log2(c.discordant_count + 1) * 5, 1)
            comps["support_sr"] = round(math.log2(c.split_count + 1) * 10, 1)
            comps["support_clip"] = round(math.log2(c.clipped_count + 1) * 2, 1)

            pa = max(c.local_nb_pvalue_a, 1e-300)
            pb = max(c.local_nb_pvalue_b, 1e-300)
            comps["pvalue_a"] = round(min(-math.log10(pa), 20), 1)
            comps["pvalue_b"] = round(min(-math.log10(pb), 20), 1)

            comps["unique_starts"] = round(min(c.unique_starts_a, c.unique_starts_b) * 0.5, 1)
            comps["reciprocal"] = 10.0 if c.reciprocal_support > 0 else 0.0
            comps["external"] = round(15.0 * len(c.external_callers), 1)
            comps["mapq"] = round(c.median_mapq * 0.2, 1)

            c.score_components = comps
            c.score = round(sum(comps.values()), 1)

    # ------------------------------------------------------------------
    # Stage 8: Output
    # ------------------------------------------------------------------

    def _stage_output(self, clusters):
        self.job.stage = ScanStage.OUTPUT
        self.emit({"type": "scan.stage_changed", "stage": "output"})

        results_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..', 'results', self.job.job_id
        )
        os.makedirs(results_dir, exist_ok=True)
        self.job.results_dir = results_dir

        from output_writer import OutputWriter
        writer = OutputWriter(results_dir, self.job.reference_build)
        writer.write_all(self.job, clusters, self._chrom_lengths)

        # Write debug evidence
        if self._debug_evidence:
            debug_path = os.path.join(results_dir, 'debug_evidence.tsv')
            cols = ["read_hash", "chrom_a", "pos_a", "chrom_b", "pos_b",
                    "mapq", "sa_mapq", "evidence_type", "cigar", "orientation",
                    "accepted", "rejection_reason"]
            with open(debug_path, 'w') as f:
                f.write('\t'.join(cols) + '\n')
                for row in self._debug_evidence:
                    f.write('\t'.join(str(row.get(c, '')) for c in cols) + '\n')


# ---------------------------------------------------------------------------
# BND VCF parsing helpers
# ---------------------------------------------------------------------------

def _parse_bnd_alt(alt: str) -> tuple[Optional[str], int]:
    """Parse BND ALT field to extract mate chrom and pos."""
    import re
    # Patterns: N[chr:pos[  N]chr:pos]  ]chr:pos]N  [chr:pos[N
    m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, 0


def _bnd_orientation(alt: str) -> str:
    """Derive orientation from BND ALT field."""
    if alt.startswith(('N[', 'n[')):
        return "++"
    elif alt.startswith(('N]', 'n]')):
        return "+-"
    elif alt.endswith(('N', 'n')) and alt.startswith(']'):
        return "-+"
    elif alt.endswith(('N', 'n')) and alt.startswith('['):
        return "--"
    return "++"
