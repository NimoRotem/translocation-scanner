"""VCF BND, BEDPE, and JSON output generation for the translocation scanner.

Produces standards-compliant output files:
- VCF 4.2 with BND notation for interchromosomal breakpoints
- BEDPE with 0-based half-open coordinates
- JSON with full cluster detail
- Summary JSON with scan-level metrics
"""
from __future__ import annotations

import json
import os
import time
from typing import Optional

from models import EvidenceCluster, ScanJob, Tier


# Canonical chromosome sort order: chr1-chr22, chrX, chrY, chrM, then anything else
_CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25})
# Also support non-prefixed names
_CHROM_ORDER.update({str(i): i for i in range(1, 23)})
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a sort key that orders chromosomes canonically.

    Chromosomes in _CHROM_ORDER sort first (by their numeric rank),
    everything else sorts after, alphabetically.
    """
    rank = _CHROM_ORDER.get(chrom)
    if rank is not None:
        return (rank, "")
    return (999, chrom)


def _sort_key_for_vcf_row(row: dict) -> tuple[int, str, int]:
    """Sort key for a VCF row dict: (chrom_rank, chrom_fallback, pos)."""
    chrom_key = _chrom_sort_key(row["chrom"])
    return (chrom_key[0], chrom_key[1], row["pos"])


def _orientation_to_strands(orientation: str) -> tuple[str, str]:
    """Convert orientation string ('++', '+-', '-+', '--') to strand pair."""
    if len(orientation) == 2 and all(c in "+-" for c in orientation):
        return orientation[0], orientation[1]
    return "+", "+"


def _bnd_alt_field(
    ref_base: str,
    mate_chrom: str,
    mate_pos: int,
    strand_local: str,
    strand_remote: str,
) -> str:
    """Build the VCF BND ALT field for a given breakend.

    BND notation encodes the orientation of each breakend:
        strand_local  strand_remote   ALT
        +             +               N[chr:pos[
        +             -               N]chr:pos]
        -             +               ]chr:pos]N
        -             -               [chr:pos[N

    Args:
        ref_base: Reference allele character (typically 'N').
        mate_chrom: Chromosome of the mate breakend.
        mate_pos: 1-based position of the mate breakend.
        strand_local: Strand of this breakend ('+' or '-').
        strand_remote: Strand of the remote breakend ('+' or '-').

    Returns:
        The formatted BND ALT string.
    """
    target = f"{mate_chrom}:{mate_pos}"
    if strand_local == "+" and strand_remote == "+":
        return f"{ref_base}[{target}["
    elif strand_local == "+" and strand_remote == "-":
        return f"{ref_base}]{target}]"
    elif strand_local == "-" and strand_remote == "+":
        return f"]{target}]{ref_base}"
    else:  # -- case
        return f"[{target}[{ref_base}"


class OutputWriter:
    """Generates VCF BND, BEDPE, JSON, and summary output files.

    All output files are written to `output_dir` which is created if it
    does not already exist.

    Coordinate conventions:
        - VCF: 1-based, inclusive positions
        - BEDPE / JSON: 0-based, half-open intervals
    """

    def __init__(self, output_dir: str, reference_build: str = "GRCh38") -> None:
        self.output_dir = output_dir
        self.reference_build = reference_build
        os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def write_all(
        self,
        job: ScanJob,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict,
    ) -> dict:
        """Write every output format and return a dict of generated file paths.

        Args:
            job: The completed ScanJob with runtime statistics.
            clusters: List of evidence clusters (calls) to write.
            chrom_lengths: Mapping of chromosome name to length in bases,
                used for the VCF contig headers.

        Returns:
            Dict with keys 'vcf', 'bedpe', 'json', 'summary' mapping to
            absolute file paths.
        """
        # Only write non-filtered clusters to main output files
        passing = [c for c in clusters if c.tier != Tier.FILTERED]
        paths = {}
        paths["vcf"] = self.write_vcf(passing, chrom_lengths)
        paths["bedpe"] = self.write_bedpe(passing)
        paths["json"] = self.write_json(job, passing)
        # Summary still gets ALL clusters for complete stats
        paths["summary"] = self.write_summary(job, clusters, paths)
        return paths

    def write_vcf(
        self,
        clusters: list[EvidenceCluster],
        chrom_lengths: dict,
    ) -> str:
        """Write a VCF 4.2 file with BND records for each translocation.

        Each translocation cluster produces two VCF rows (mate breakends).
        Rows are sorted by chromosome order then position.

        Args:
            clusters: Evidence clusters to emit.
            chrom_lengths: Chromosome name -> length for contig headers.

        Returns:
            Absolute path to the written VCF file.
        """
        path = os.path.join(self.output_dir, "translocations.vcf")
        rows: list[dict] = []

        for idx, cluster in enumerate(clusters, start=1):
            label = f"TRA_{idx:03d}"
            id_a = f"{label}_a"
            id_b = f"{label}_b"
            strand_a, strand_b = _orientation_to_strands(cluster.orientation)

            # Positions in VCF are 1-based.  The model stores 0-based
            # positions, so we add 1 for VCF output.
            pos_a_vcf = cluster.pos_a + 1
            pos_b_vcf = cluster.pos_b + 1

            # Determine FILTER value
            filt = self._vcf_filter(cluster)

            # Shared INFO fields
            info_common = self._vcf_info_fields(cluster)

            # Row A: breakend on chrom_a pointing to chrom_b
            alt_a = _bnd_alt_field("N", cluster.chrom_b, pos_b_vcf, strand_a, strand_b)
            info_a = f"SVTYPE=BND;MATEID={id_b};{info_common}"
            rows.append({
                "chrom": cluster.chrom_a,
                "pos": pos_a_vcf,
                "id": id_a,
                "ref": "N",
                "alt": alt_a,
                "qual": ".",
                "filter": filt,
                "info": info_a,
            })

            # Row B: breakend on chrom_b pointing back to chrom_a
            alt_b = _bnd_alt_field("N", cluster.chrom_a, pos_a_vcf, strand_b, strand_a)
            info_b = f"SVTYPE=BND;MATEID={id_a};{info_common}"
            rows.append({
                "chrom": cluster.chrom_b,
                "pos": pos_b_vcf,
                "id": id_b,
                "ref": "N",
                "alt": alt_b,
                "qual": ".",
                "filter": filt,
                "info": info_b,
            })

        # Sort rows by chromosome order then position
        rows.sort(key=_sort_key_for_vcf_row)

        with open(path, "w") as fh:
            self._write_vcf_header(fh, chrom_lengths)
            for row in rows:
                line = "\t".join([
                    row["chrom"],
                    str(row["pos"]),
                    row["id"],
                    row["ref"],
                    row["alt"],
                    row["qual"],
                    row["filter"],
                    row["info"],
                ])
                fh.write(line + "\n")

        return os.path.abspath(path)

    def write_bedpe(self, clusters: list[EvidenceCluster]) -> str:
        """Write a BEDPE file with 0-based half-open coordinates.

        Each cluster produces one BEDPE row.  Rows are sorted by
        chromosome order of chrom1, then start1.

        Args:
            clusters: Evidence clusters to emit.

        Returns:
            Absolute path to the written BEDPE file.
        """
        path = os.path.join(self.output_dir, "translocations.bedpe")

        bedpe_rows: list[dict] = []
        for idx, cluster in enumerate(clusters, start=1):
            strand_a, strand_b = _orientation_to_strands(cluster.orientation)
            name = f"TRA_{idx:03d}"
            # 0-based half-open: start = pos, end = pos + 1
            bedpe_rows.append({
                "chrom1": cluster.chrom_a,
                "start1": cluster.pos_a,
                "end1": cluster.pos_a + 1,
                "chrom2": cluster.chrom_b,
                "start2": cluster.pos_b,
                "end2": cluster.pos_b + 1,
                "name": name,
                "score": f"{cluster.score:.4f}",
                "strand1": strand_a,
                "strand2": strand_b,
                "tier": cluster.tier.value,
                "support": str(cluster.total_support),
                "bg_pval": f"{cluster.background_p:.6e}",
            })

        # Sort by chromosome order of chrom1 then start1
        bedpe_rows.sort(key=lambda r: (
            _chrom_sort_key(r["chrom1"]),
            r["start1"],
            _chrom_sort_key(r["chrom2"]),
            r["start2"],
        ))

        with open(path, "w") as fh:
            fh.write(
                "#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\t"
                "name\tscore\tstrand1\tstrand2\ttier\tsupport\tbg_pval\n"
            )
            for row in bedpe_rows:
                fields = [
                    row["chrom1"],
                    str(row["start1"]),
                    str(row["end1"]),
                    row["chrom2"],
                    str(row["start2"]),
                    str(row["end2"]),
                    row["name"],
                    row["score"],
                    row["strand1"],
                    row["strand2"],
                    row["tier"],
                    row["support"],
                    row["bg_pval"],
                ]
                fh.write("\t".join(fields) + "\n")

        return os.path.abspath(path)

    def write_json(
        self,
        job: ScanJob,
        clusters: list[EvidenceCluster],
    ) -> str:
        """Write a JSON file with full cluster detail.

        All coordinates are 0-based. Clusters are sorted by chromosome
        order of chrom_a then pos_a.

        Args:
            job: The ScanJob for sample metadata.
            clusters: Evidence clusters to emit.

        Returns:
            Absolute path to the written JSON file.
        """
        path = os.path.join(self.output_dir, "translocations.json")

        sorted_clusters = sorted(
            clusters,
            key=lambda c: (
                _chrom_sort_key(c.chrom_a),
                c.pos_a,
                _chrom_sort_key(c.chrom_b),
                c.pos_b,
            ),
        )

        calls = [c.to_dict() for c in sorted_clusters]

        # Build scan summary from job
        elapsed = 0.0
        if job.started_at:
            end = job.completed_at or time.time()
            elapsed = round(end - job.started_at, 1)

        data = {
            "coordinate_system": "0-based",
            "reference_build": self.reference_build,
            "sample": os.path.basename(job.file_path),
            "scan_summary": {
                "job_id": job.job_id,
                "total_reads": job.total_reads,
                "chimeric_rate": job.chimeric_rate,
                "insert_size_median": job.insert_size_median,
                "insert_size_std": job.insert_size_std,
                "discordant_count": job.discordant_count,
                "split_count": job.split_count,
                "clip_count": job.clip_count,
                "elapsed_seconds": elapsed,
                "num_calls": len(calls),
            },
            "calls": calls,
        }

        with open(path, "w") as fh:
            json.dump(data, fh, indent=2)

        return os.path.abspath(path)

    def write_summary(
        self,
        job: ScanJob,
        clusters: list[EvidenceCluster],
        file_paths: Optional[dict] = None,
    ) -> str:
        """Write a compact summary JSON with scan-level metrics.

        Args:
            job: The completed ScanJob.
            clusters: Evidence clusters produced by the pipeline.
            file_paths: Optional dict of already-generated file paths
                (keys: 'vcf', 'bedpe', 'json') to include in the summary.

        Returns:
            Absolute path to the written summary JSON file.
        """
        path = os.path.join(self.output_dir, "summary.json")

        # Count calls by tier
        calls_by_tier: dict[str, int] = {}
        for tier_val in Tier:
            calls_by_tier[tier_val.value] = 0
        for cluster in clusters:
            calls_by_tier[cluster.tier.value] = (
                calls_by_tier.get(cluster.tier.value, 0) + 1
            )

        elapsed = 0.0
        if job.started_at:
            end = job.completed_at or time.time()
            elapsed = round(end - job.started_at, 1)

        summary = {
            "job_id": job.job_id,
            "file_path": job.file_path,
            "reference_build": self.reference_build,
            "total_reads": job.total_reads,
            "chimeric_rate": job.chimeric_rate,
            "insert_size_median": job.insert_size_median,
            "num_calls": len(clusters),
            "calls_by_tier": calls_by_tier,
            "elapsed_seconds": elapsed,
            "files": file_paths or {},
        }

        with open(path, "w") as fh:
            json.dump(summary, fh, indent=2)

        return os.path.abspath(path)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _write_vcf_header(self, fh, chrom_lengths: dict) -> None:
        """Write the VCF 4.2 header block including meta-information lines.

        Args:
            fh: Open file handle to write to.
            chrom_lengths: Mapping of chromosome name to length for
                ##contig lines.
        """
        lines = [
            "##fileformat=VCFv4.2",
            f"##reference={self.reference_build}",
            # INFO fields
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">',
            '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CI around POS">',
            '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Total supporting reads">',
            '##INFO=<ID=DISCORDANT,Number=1,Type=Integer,Description="Discordant pair support">',
            '##INFO=<ID=SPLIT,Number=1,Type=Integer,Description="Split read support">',
            '##INFO=<ID=CLIPPED,Number=1,Type=Integer,Description="Clipped read support">',
            '##INFO=<ID=BG_PVAL,Number=1,Type=Float,Description="Background model p-value">',
            '##INFO=<ID=TIER,Number=1,Type=String,Description="Call tier">',
            '##INFO=<ID=SCORE,Number=1,Type=Float,Description="Raw score">',
            # FILTER fields
            '##FILTER=<ID=LOW_MAPQ,Description="Median MAPQ below 20">',
            '##FILTER=<ID=BG_NOISE,Description="Background p-value above 0.001">',
            '##FILTER=<ID=CENTROMERE,Description="Both breakpoints near centromere">',
        ]

        # Add contig lines sorted by canonical chromosome order
        sorted_chroms = sorted(chrom_lengths.keys(), key=_chrom_sort_key)
        for chrom in sorted_chroms:
            length = chrom_lengths[chrom]
            lines.append(f"##contig=<ID={chrom},length={length}>")

        # Column header
        lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

        for line in lines:
            fh.write(line + "\n")

    def _vcf_filter(self, cluster: EvidenceCluster) -> str:
        """Determine the VCF FILTER column value for a cluster.

        Clusters with tier CONFIRMED, VALIDATED, or LIKELY get PASS.
        FILTERED or CANDIDATE clusters list the specific filter flags
        that were applied.  If a cluster has no filter_flags but is not
        passing, it gets a '.' (missing).

        Args:
            cluster: The evidence cluster to evaluate.

        Returns:
            FILTER column string ('PASS', semicolon-joined flags, or '.').
        """
        if cluster.tier in (Tier.CONFIRMED, Tier.VALIDATED, Tier.LIKELY):
            return "PASS"

        if cluster.filter_flags:
            return ";".join(cluster.filter_flags)

        return "."

    def _vcf_info_fields(self, cluster: EvidenceCluster) -> str:
        """Build the shared INFO fields string (everything except SVTYPE and MATEID).

        Args:
            cluster: The evidence cluster.

        Returns:
            Semicolon-delimited INFO key=value pairs.
        """
        ci_a = f"{cluster.ci_a[0]},{cluster.ci_a[1]}"
        ci_b = f"{cluster.ci_b[0]},{cluster.ci_b[1]}"

        fields = [
            f"CIPOS={ci_a}",
            f"SUPPORT={cluster.total_support}",
            f"DISCORDANT={cluster.discordant_count}",
            f"SPLIT={cluster.split_count}",
            f"CLIPPED={cluster.clipped_count}",
            f"BG_PVAL={cluster.background_p:.6e}",
            f"TIER={cluster.tier.value}",
            f"SCORE={cluster.score:.4f}",
        ]

        return ";".join(fields)
