# Hardware & Reference Audit

## Server: genom-beast-gpu
- **Instance**: GCE us-central1-c, IP 34.135.47.236
- **CPU**: 32 cores, Intel Xeon @ 2.30GHz
- **RAM**: 117 GB total (~16 GB used typical)
- **Disk**: 2.9 TB SSD, 2.0 TB free (mounted at /dev/sda1)
- **GPU**: Tesla T4, 15360 MiB VRAM
- **OS**: Debian (genom-beast-gpu)

## Resource Budget
- Pipeline stages: 28 cores, 100 GB RAM
- OS/service headroom: 4 cores, 17 GB RAM
- DELLY: OMP_NUM_THREADS=14
- Manta: REMOVED (requires Python 2, not available)
- samtools threading: -@ 8 for sort/index, -@ 4 for view
- Extraction workers: min(nproc - 4, 24) = 24

## Reference FASTA
- **Path**: /data/refs/hs38DH.fa
- **SHA256**: 3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50
- **Contig count**: 3366 (chr-prefixed, includes HLA/alt/decoy)
- **Primary contigs**: chr1-22, chrX, chrY (24 total)
- **Symlink**: /data/refs/GRCh38.fa -> hs38DH.fa
- **Numeric reference**: /data/genom-nimo/reference.fasta (1-22, X, Y naming)

## Tools (conda genomics env)
- **Python**: /home/nimo/miniconda3/envs/genomics/bin/python 3.13.12
- **samtools**: present (lib warning but functional)
- **bcftools**: 1.22
- **DELLY**: v1.7.3
- **Manta**: configManta.py present but BROKEN (requires python2, not installed)
- **minimap2**: check availability before clip_realignment stage

## BAM Files
- Nimo.bam: 93 GB, indexed (.bai present)
- Also: B2XH (58G), B3XH (58G), Chichi (57G), Efi (54G), Mina (52G)

## Exclude BED
- **Path**: /data/masks/exclude_grch38.bed
- **Contents**: All non-primary contigs from hs38DH.fa (3342 entries)
- **Used by**: DELLY (-x flag)
- **Generated**: At deploy time from reference .fai
