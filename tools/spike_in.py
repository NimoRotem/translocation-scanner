#!/usr/bin/env python3
"""Create a BAM with injected translocation-supporting reads.

Two modes:
  --mode standalone  Create a small BAM with background + spike-in (fast, for unit tests)
  --mode merge       Create spike-in reads only, then merge into the full source BAM
                     (slow due to IO, but gives realistic background for validation)

Usage:
    # Fast standalone mode (for quick iteration):
    python spike_in.py --source /data/aligned_bams/Nimo.bam \
        --reference /data/genom-nimo/reference.fasta \
        --output /tmp/spike_standalone.bam --mode standalone

    # Full merge mode (for realistic validation):
    python spike_in.py --source /data/aligned_bams/Nimo.bam \
        --reference /data/genom-nimo/reference.fasta \
        --output /tmp/spike_merged.bam --mode merge
"""

import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile

import pysam


def random_seq(length: int) -> str:
    return ''.join(random.choice('ACGT') for _ in range(length))


def random_qual(length: int, min_q: int = 30, max_q: int = 40) -> list[int]:
    return [random.randint(min_q, max_q) for _ in range(length)]


def make_read_name(prefix: str, idx: int) -> str:
    h = hashlib.md5(f"{prefix}_{idx}".encode()).hexdigest()[:16]
    return f"SPIKE_{prefix}_{h}"


def fetch_ref(fasta: pysam.FastaFile, chrom: str, start: int, length: int) -> str:
    try:
        seq = fasta.fetch(chrom, max(0, start), start + length).upper()
        if len(seq) < length:
            seq += random_seq(length - len(seq))
        return seq
    except (KeyError, ValueError):
        return random_seq(length)


def inject_discordant_pair(
    header, fasta, chrom_a, pos_a, chrom_b, pos_b,
    read_name, read_length=150, mapq=60,
):
    """Discordant read pair: R1 at breakpoint A, R2 at breakpoint B."""
    ja = random.randint(-200, 200)
    jb = random.randint(-200, 200)
    pa = max(0, pos_a + ja)
    pb = max(0, pos_b + jb)

    tid_a = header.get_tid(chrom_a)
    tid_b = header.get_tid(chrom_b)
    if tid_a < 0 or tid_b < 0:
        raise ValueError(f"Chrom not in header: {chrom_a}={tid_a}, {chrom_b}={tid_b}")

    seq_a = fetch_ref(fasta, chrom_a, pa, read_length)
    seq_b = fetch_ref(fasta, chrom_b, pb, read_length)

    r1 = pysam.AlignedSegment(header)
    r1.query_name = read_name
    r1.query_sequence = seq_a
    r1.flag = 0x1 | 0x40 | 0x20  # paired, first-in-pair, mate-reverse
    r1.reference_id = tid_a
    r1.reference_start = pa
    r1.mapping_quality = mapq
    r1.cigar = [(0, read_length)]
    r1.query_qualities = pysam.qualitystring_to_array('I' * read_length)
    r1.next_reference_id = tid_b
    r1.next_reference_start = pb
    r1.template_length = 0

    r2 = pysam.AlignedSegment(header)
    r2.query_name = read_name
    r2.query_sequence = seq_b
    r2.flag = 0x1 | 0x80 | 0x10  # paired, second-in-pair, self-reverse
    r2.reference_id = tid_b
    r2.reference_start = pb
    r2.mapping_quality = mapq
    r2.cigar = [(0, read_length)]
    r2.query_qualities = pysam.qualitystring_to_array('I' * read_length)
    r2.next_reference_id = tid_a
    r2.next_reference_start = pa
    r2.template_length = 0

    return r1, r2


def inject_split_read(
    header, fasta, chrom_a, pos_a, chrom_b, pos_b,
    read_name, read_length=150, split_pos=75, mapq=60,
):
    """Split read: primary maps partly to A with SA tag pointing to B.

    The SA tag CIGAR must match the supplementary alignment exactly.
    Primary: first `split_pos` bases map to A, rest clipped (75M75S)
    Supplementary: first `split_pos` bases clipped, rest maps to B (75S75M)
    """
    ja = random.randint(-50, 50)
    jb = random.randint(-50, 50)
    pa = max(0, pos_a + ja)
    pb = max(0, pos_b + jb)

    tid_a = header.get_tid(chrom_a)
    tid_b = header.get_tid(chrom_b)

    seq_part_a = fetch_ref(fasta, chrom_a, pa, split_pos)
    seq_part_b = fetch_ref(fasta, chrom_b, pb, read_length - split_pos)
    full_seq = seq_part_a + seq_part_b

    qual = pysam.qualitystring_to_array('I' * read_length)

    # Primary alignment: first half maps to A, second half soft-clipped
    primary = pysam.AlignedSegment(header)
    primary.query_name = read_name
    primary.query_sequence = full_seq
    primary.flag = 0  # unpaired primary (simpler — extractor checks SA tag)
    primary.reference_id = tid_a
    primary.reference_start = pa
    primary.mapping_quality = mapq
    primary.cigar = [(0, split_pos), (4, read_length - split_pos)]
    primary.query_qualities = qual
    # SA tag must describe the supplementary alignment exactly
    supp_cigar = f"{split_pos}S{read_length - split_pos}M"
    primary.set_tag('SA', f"{chrom_b},{pb + 1},+,{supp_cigar},{mapq},0;")

    # Supplementary alignment: first half clipped, second half maps to B
    supp = pysam.AlignedSegment(header)
    supp.query_name = read_name
    supp.query_sequence = full_seq
    supp.flag = 0x800  # supplementary
    supp.reference_id = tid_b
    supp.reference_start = pb
    supp.mapping_quality = mapq
    supp.cigar = [(4, split_pos), (0, read_length - split_pos)]
    supp.query_qualities = qual
    pri_cigar = f"{split_pos}M{read_length - split_pos}S"
    supp.set_tag('SA', f"{chrom_a},{pa + 1},+,{pri_cigar},{mapq},0;")

    return primary, supp


def inject_clipped_read(
    header, fasta, chrom_a, pos_a, chrom_b, pos_b,
    read_name, read_length=150, clip_len=40, mapq=60,
):
    """Soft-clipped read at breakpoint (no SA tag — clip-only evidence)."""
    actual_pos = max(0, pos_a + random.randint(-30, 30))
    tid_a = header.get_tid(chrom_a)

    aligned_len = read_length - clip_len
    seq_aligned = fetch_ref(fasta, chrom_a, actual_pos, aligned_len)
    seq_clipped = fetch_ref(fasta, chrom_b, max(0, pos_b + random.randint(-30, 30)), clip_len)
    full_seq = seq_aligned + seq_clipped

    read = pysam.AlignedSegment(header)
    read.query_name = read_name
    read.query_sequence = full_seq
    read.flag = 0
    read.reference_id = tid_a
    read.reference_start = actual_pos
    read.mapping_quality = mapq
    read.cigar = [(0, aligned_len), (4, clip_len)]
    read.query_qualities = pysam.qualitystring_to_array('I' * read_length)

    return read


def generate_spike_reads(header, fasta, events):
    """Generate all spike-in reads for a list of translocation events."""
    reads = []
    for event in events:
        ca, pa = event["chrom_a"], event["pos_a"]
        cb, pb = event["chrom_b"], event["pos_b"]
        label = event.get("label", f"{ca}:{pa}-{cb}:{pb}")
        n_disc = event.get("discordant", 10)
        n_split = event.get("split", 5)
        n_clip = event.get("clipped", 3)

        print(f"  Generating {label}: {n_disc}d + {n_split}s + {n_clip}c")

        for i in range(n_disc):
            name = make_read_name(f"disc_{label}", i)
            r1, r2 = inject_discordant_pair(header, fasta, ca, pa, cb, pb, name)
            reads.extend([r1, r2])

        for i in range(n_split):
            name = make_read_name(f"split_{label}", i)
            pri, sup = inject_split_read(header, fasta, ca, pa, cb, pb, name)
            reads.extend([pri, sup])

        for i in range(n_clip):
            name = make_read_name(f"clip_{label}", i)
            r = inject_clipped_read(header, fasta, ca, pa, cb, pb, name)
            reads.append(r)

    return reads


def create_spike_only_bam(header, spike_reads, output_path):
    """Write spike-in reads to a sorted, indexed BAM."""
    tmp = output_path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for r in spike_reads:
            out.write(r)
    pysam.sort("-@", "4", "-o", output_path, tmp)
    pysam.index("-@", "4", output_path)
    os.unlink(tmp)
    return output_path


def mode_standalone(source_path, ref_path, output_path, events, seed):
    """Create standalone BAM with LOCAL background around breakpoints + spike-in.

    For each event, extracts 5 Mb of real reads around both breakpoints
    from the source BAM.  Also extracts background from 3 random regions
    per chromosome to give the statistical model realistic baselines.
    """
    random.seed(seed)
    fasta = pysam.FastaFile(ref_path)
    source = pysam.AlignmentFile(source_path, "rb")
    header = source.header.copy()

    # 1. Extract LOCAL background around each breakpoint (5 Mb window)
    bg_reads = []
    seen_regions = set()
    margin = 2_500_000  # 2.5 Mb on each side

    for event in events:
        for side in [("chrom_a", "pos_a"), ("chrom_b", "pos_b")]:
            chrom = event[side[0]]
            pos = event[side[1]]
            start = max(0, pos - margin)
            end = pos + margin
            key = f"{chrom}:{start}-{end}"
            if key in seen_regions:
                continue
            seen_regions.add(key)
            try:
                count = 0
                for read in source.fetch(chrom, start, end):
                    bg_reads.append(read)
                    count += 1
                print(f"  Local background {chrom}:{start}-{end}: {count} reads")
            except Exception as e:
                print(f"  Warning: could not fetch {key}: {e}")

    # 2. Add general background from each chromosome (1 Mb mid-chromosome)
    for chrom in [str(c) for c in range(1, 23)] + ['X']:
        try:
            clen = header.get_reference_length(chrom)
            mid = clen // 2
            count = 0
            for read in source.fetch(chrom, mid - 500_000, mid + 500_000):
                if count >= 10000:
                    break
                bg_reads.append(read)
                count += 1
        except Exception:
            pass

    print(f"Total background reads: {len(bg_reads)}")

    spike_reads = generate_spike_reads(header, fasta, events)
    print(f"Generated {len(spike_reads)} spike-in reads")

    # Write combined
    tmp = output_path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for r in bg_reads:
            out.write(r)
        for r in spike_reads:
            out.write(r)

    print("Sorting and indexing...")
    pysam.sort("-@", "4", "-o", output_path, tmp)
    pysam.index("-@", "4", output_path)
    os.unlink(tmp)

    total = len(bg_reads) + len(spike_reads)
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Done: {output_path} ({total} reads, {size_mb:.0f} MB)")
    source.close()
    fasta.close()


def mode_merge(source_path, ref_path, output_path, events, seed):
    """Merge spike-in reads into the full source BAM."""
    random.seed(seed)
    fasta = pysam.FastaFile(ref_path)
    source = pysam.AlignmentFile(source_path, "rb")
    header = source.header.copy()
    source.close()

    spike_reads = generate_spike_reads(header, fasta, events)
    print(f"Generated {len(spike_reads)} spike-in reads")
    fasta.close()

    # Write spike-only BAM
    spike_bam = tempfile.mktemp(suffix="_spike.bam")
    create_spike_only_bam(header, spike_reads, spike_bam)
    print(f"Spike-only BAM: {spike_bam}")

    # Merge using samtools (preserves all original reads + adds spike-in)
    print(f"Merging with {source_path} (this will take a while for large BAMs)...")
    cmd = [
        "samtools", "merge", "-@", "4", "-f",
        output_path, source_path, spike_bam,
    ]
    subprocess.run(cmd, check=True)
    print("Indexing merged BAM...")
    subprocess.run(["samtools", "index", "-@", "4", output_path], check=True)

    # Cleanup
    os.unlink(spike_bam)
    os.unlink(spike_bam + ".bai")
    print(f"Done: {output_path}")


DEFAULT_EVENTS = [
    {
        "chrom_a": "9", "pos_a": 130854064,
        "chrom_b": "22", "pos_b": 23632600,
        "label": "BCR-ABL_t(9;22)",
        "discordant": 20, "split": 10, "clipped": 5,
    },
    {
        "chrom_a": "15", "pos_a": 73027451,
        "chrom_b": "17", "pos_b": 40343345,
        "label": "PML-RARA_t(15;17)",
        "discordant": 15, "split": 8, "clipped": 3,
    },
    {
        "chrom_a": "11", "pos_a": 118353210,
        "chrom_b": "14", "pos_b": 105862702,
        "label": "IGH-CCND1_t(11;14)",
        "discordant": 8, "split": 5, "clipped": 2,
    },
    {
        "chrom_a": "8", "pos_a": 128750686,
        "chrom_b": "14", "pos_b": 105862702,
        "label": "IGH-MYC_t(8;14)",
        "discordant": 8, "split": 5, "clipped": 2,
    },
    {
        "chrom_a": "2", "pos_a": 29415640,
        "chrom_b": "5", "pos_b": 170837543,
        "label": "ALK-NPM1_t(2;5)",
        "discordant": 5, "split": 3, "clipped": 1,
    },
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translocation scanner spike-in BAM generator")
    parser.add_argument("--source", required=True, help="Source BAM for background/header")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--output", required=True, help="Output BAM path")
    parser.add_argument("--events", help="JSON file with events (defaults to built-in set)")
    parser.add_argument("--mode", choices=["standalone", "merge"], default="standalone",
                        help="standalone=small BAM, merge=full BAM with spike-in")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    if args.events:
        with open(args.events) as f:
            events = json.load(f)
    else:
        events = DEFAULT_EVENTS

    if args.mode == "standalone":
        mode_standalone(args.source, args.reference, args.output, events, args.seed)
    else:
        mode_merge(args.source, args.reference, args.output, events, args.seed)
