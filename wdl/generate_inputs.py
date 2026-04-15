#!/usr/bin/env python3
"""
Generate MEGA-xTEA WDL inputs.json from a BAM directory or sample list.

Usage:
  # From a directory of BAMs:
  python generate_inputs.py --bam-dir /path/to/bams/ --ref-dir /path/to/reference/ \
      --kmer-dir /path/to/kmer/ --ref-name GRCh38.fa > inputs.json

  # From a sample list file (one sample name per line, BAMs named {sample}.bam):
  python generate_inputs.py --sample-list samples.txt --bam-dir /path/to/bams/ \
      --ref-dir /path/to/reference/ --kmer-dir /path/to/kmer/ --ref-name GRCh38.fa > inputs.json
"""

import argparse
import json
import os
import sys
from pathlib import Path


def find_bam_index(bam_path: str) -> str:
    """Find the BAM index file (.bai or .bam.bai)."""
    for suffix in [".bai", ".bam.bai"]:
        candidate = bam_path.replace(".bam", suffix) if suffix == ".bai" else bam_path + ".bai"
        if os.path.exists(candidate):
            return candidate
    # Default: assume .bam.bai
    return bam_path + ".bai"


def main():
    parser = argparse.ArgumentParser(description="Generate MEGA-xTEA WDL inputs.json")
    parser.add_argument("--bam-dir", required=True, help="Directory containing BAM files")
    parser.add_argument("--ref-dir", required=True, help="Directory containing reference + BLAST DB + RM output")
    parser.add_argument("--kmer-dir", required=True, help="Directory containing .mk and .mi files")
    parser.add_argument("--ref-name", default="GRCh38.fa", help="Reference FASTA filename (default: GRCh38.fa)")
    parser.add_argument("--kmer-prefix", default=None, help="K-mer file prefix (default: ref name without .fa)")
    parser.add_argument("--repeat-lib", default=None, help="Repeat library path (default: ref-dir/Dfam3.2_human.rep)")
    parser.add_argument("--sample-list", default=None, help="File with sample names (one per line)")
    parser.add_argument("--threads", type=int, default=8, help="Threads per sample")
    parser.add_argument("--memory-gb", type=int, default=16, help="Memory (GB) per sample")
    parser.add_argument("--disk-gb", type=int, default=100, help="Disk (GB) per sample")
    parser.add_argument("--docker", default="mega-xtea:latest", help="Docker image")
    parser.add_argument("--no-sva-filter", action="store_true", help="Disable SVA filtering")
    parser.add_argument("--gaussian-genotype", action="store_true", help="Use Gaussian instead of ML")
    parser.add_argument("--no-deletion", action="store_true", help="Disable deletion detection")
    parser.add_argument("--model-tar", default=None, help="Deep Forest model .tar.gz path")
    parser.add_argument("--model-pkl", default=None, help="Sklearn .pkl model path")
    args = parser.parse_args()

    bam_dir = os.path.abspath(args.bam_dir)
    ref_dir = os.path.abspath(args.ref_dir)
    kmer_dir = os.path.abspath(args.kmer_dir)
    ref_name = args.ref_name
    kmer_prefix = args.kmer_prefix or ref_name.replace(".fa", "").replace(".fasta", "")

    # Discover samples
    if args.sample_list:
        with open(args.sample_list) as f:
            samples = [line.strip() for line in f if line.strip()]
    else:
        samples = sorted(
            Path(p).stem
            for p in os.listdir(bam_dir)
            if p.endswith(".bam") and not p.endswith(".unmapped.bam")
        )

    if not samples:
        print("ERROR: No samples found.", file=sys.stderr)
        sys.exit(1)

    print(f"# Found {len(samples)} samples", file=sys.stderr)

    bams = [os.path.join(bam_dir, f"{s}.bam") for s in samples]
    bais = [find_bam_index(b) for b in bams]

    repeat_lib = args.repeat_lib or os.path.join(ref_dir, "Dfam3.2_human.rep")

    inputs = {
        "mega_xtea_batch.input_bams": bams,
        "mega_xtea_batch.input_bam_indices": bais,
        "mega_xtea_batch.sample_names": samples,
        "mega_xtea_batch.reference_fasta": os.path.join(ref_dir, ref_name),
        "mega_xtea_batch.reference_fasta_index": os.path.join(ref_dir, f"{ref_name}.fai"),
        "mega_xtea_batch.reference_blast_nhr": os.path.join(ref_dir, f"{ref_name}.nhr"),
        "mega_xtea_batch.reference_blast_nin": os.path.join(ref_dir, f"{ref_name}.nin"),
        "mega_xtea_batch.reference_blast_nsq": os.path.join(ref_dir, f"{ref_name}.nsq"),
        "mega_xtea_batch.repeat_masker_out": os.path.join(ref_dir, f"{ref_name}.out"),
        "mega_xtea_batch.repeat_library": repeat_lib,
        "mega_xtea_batch.kmer_mk": os.path.join(kmer_dir, f"{kmer_prefix}.mk"),
        "mega_xtea_batch.kmer_mi": os.path.join(kmer_dir, f"{kmer_prefix}.mi"),
        "mega_xtea_batch.sva_filter": not args.no_sva_filter,
        "mega_xtea_batch.ml_genotype": not args.gaussian_genotype,
        "mega_xtea_batch.detect_deletion": not args.no_deletion,
        "mega_xtea_batch.threads_per_sample": args.threads,
        "mega_xtea_batch.memory_gb_per_sample": args.memory_gb,
        "mega_xtea_batch.disk_gb_per_sample": args.disk_gb,
        "mega_xtea_batch.docker_image": args.docker,
    }

    if args.model_tar:
        inputs["mega_xtea_batch.ml_model_tar"] = os.path.abspath(args.model_tar)
    if args.model_pkl:
        inputs["mega_xtea_batch.ml_model_pkl"] = os.path.abspath(args.model_pkl)

    json.dump(inputs, sys.stdout, indent=4)
    print()  # trailing newline


if __name__ == "__main__":
    main()
