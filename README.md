# MEGA-xTEA

## Overview

MEGA-xTEA merges the best of two leading transposable element detection tools:

- **MEGAnE** (speed): C++ multi-threaded BAM processing, k-mer filtering, deletion detection
- **xTea** (accuracy): SVA VNTR-aware filtering, ML-based genotyping, transduction detection

Designed for **high-speed, large-scale germline MEI detection** in population cohorts (>10,000 samples).

### Key Features

- 8-10x faster than xTea with comparable SVA detection accuracy
- TE insertion detection: LINE1, Alu, SVA, HERV-K
- TE deletion detection (Absent ME) -- unique capability from MEGAnE
- SVA-specific VNTR-aware filtering from xTea
- ML-based genotyping (Random Forest) with Gaussian fallback
- Scalable joint calling for massive cohorts (>10,000 samples)
- VCF output compatible with phasing/imputation tools

## Architecture

MEGA-xTEA combines MEGAnE's C++ read-extraction engine with xTea's Python-based
post-filtering and genotyping modules into a unified pipeline:

```
BAM/CRAM ──> [C++ Core: MEGAnE] ──> Candidate Sites
                   |                        |
                   |-- extract_discordant    |-- SVA Post-Filter [xTea]
                   |-- extract_unmapped      |-- ML Genotyping [xTea]
                   +-- k-mer filtering       +-- Transduction [xTea]
                                                    |
                                             ──> Filtered VCF
                                                    |
              [Deletion Detection: MEGAnE] ────> Absent ME calls
                                                    |
              [Joint Calling: MEGAnE] ─────> Population VCF
```

The five compiled C++ shared-object modules handle the performance-critical work:

| Module | Purpose |
|--------|---------|
| `extract_discordant.so` | Extract discordant read pairs from BAM/CRAM |
| `extract_unmapped.so` | Extract unmapped reads with soft-clipped anchors |
| `convert_rep_to_2bit_k11.so` | Build 2-bit encoded k=11 repeat k-mer set |
| `remove_multimapping_reads_from_fa.so` | Remove ambiguous multi-mapping reads |
| `save_redundant_kmers.so` | Identify and store redundant k-mers for filtering |

Python modules in `megaxtea/` provide the accuracy-critical post-processing:

| Module | Origin | Purpose |
|--------|--------|---------|
| `sva_filter` | xTea | SVA VNTR-aware consistency relaxation and filtering |
| `ml_genotype` | xTea | Random Forest genotyping (15-dimensional features) |
| `polyA_detector` | xTea | Structure-aware polyA signal detection |
| `transduction` | xTea | 3' transduction detection |
| `te_classifier` | MEGAnE + xTea | Enhanced TE classification with binary encoding |
| `config` | Both | Central parameter configuration |

## Requirements

- Python >= 3.8
- GCC >= 7.0 with C++11 support
- htslib >= 1.10
- BLAST+ (`ncbi-blast+`)
- bedtools >= 2.29
- samtools >= 1.10
- Python packages: numpy, scipy, scikit-learn, pysam

## Installation

### Docker (Recommended)

```bash
git clone https://github.com/Mirror-fish/MEGA-xTEA.git
cd MEGA-xTEA
docker build -t mega-xtea:latest .

# Verify
docker run --rm mega-xtea:latest --version
```

See [docs/docker-guide.md](docs/docker-guide.md) for complete Docker usage including batch processing, Singularity conversion, Snakemake integration, and HPC deployment.

### Native Install

```bash
git clone https://github.com/Mirror-fish/MEGA-xTEA.git
cd MEGA-xTEA
make                        # Compile C++ modules
pip install -r requirements.txt
```

### Building C++ modules

The C++ modules link against htslib. Make sure htslib is installed and its
headers/libraries are available under `external/htslib` (or adjust the Makefile
`LHTS` variable):

```bash
make clean && make
# Verify that all five shared objects were built
ls cpp/*.so
```

Expected output:

```
cpp/convert_rep_to_2bit_k11.so
cpp/extract_discordant.so
cpp/extract_unmapped.so
cpp/remove_multimapping_reads_from_fa.so
cpp/save_redundant_kmers.so
```

### Verifying installation

```bash
python mega-xtea.py --help
```

## Quick Start

### Step 0: Build k-mer set (one-time per reference genome)

```bash
python mega-xtea.py build-kmer \
  -r /path/to/GRCh38.fa \
  -o /path/to/kmer_set/ \
  -t 4
```

### Step 1: Individual calling

```bash
python mega-xtea.py call \
  -i sample.bam \
  -r GRCh38.fa \
  -k /path/to/kmer_set/kmer.mk \
  -R /path/to/repeat_lib/ \
  -o /path/to/output/ \
  -t 8 \
  --sva-filter \
  --ml-genotype \
  --detect-deletion
```

### Step 2: Joint calling (multi-sample)

```bash
python mega-xtea.py joint-call \
  --vcf-list sample_vcfs.txt \
  -o /path/to/joint_output/ \
  -t 8
```

### Step 3: Reshape VCF for imputation (optional)

```bash
python mega-xtea.py reshape-vcf \
  -i joint_output.vcf.gz \
  -o reshaped.vcf.gz
```

## Detailed Usage

### Controlling SVA filtering

The xTea-derived SVA filter is enabled by default. It applies VNTR-aware
consistency relaxation (cluster-diff cutoff relaxed to 200 bp for SVA) and
structure-aware polyA detection. To disable:

```bash
python mega-xtea.py call ... --no-sva-filter
```

Key SVA filter parameters (see `megaxtea/config.py` `SVAParams`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `REP_SVA_CNS_HEAD` | 400 | SVA consensus 5' head length |
| `REP_SVA_POLYA_START` | 1900 | PolyA region start in SVA consensus |
| `REP_SVA_MIN_LEN` | 300 | Minimum SVA insertion length |
| `sva_clip_cluster_diff_cutoff` | 200 | Relaxed cluster-diff for SVA VNTR regions |

### Genotyping methods

Two methods are available:

- `--ml-genotype` (default): xTea's Random Forest classifier using 15-dimensional
  feature vectors (clip counts, discordant reads, coverage, polyA signals).
  Trained with 20 estimators and 70/30 train-test split.
- `--gaussian-genotype`: MEGAnE's three-component Gaussian mixture fitting over
  allele-evidence reads (faster, less accurate). Requires a minimum of 3
  evidence reads per site.

### Deletion detection

Detect absence of reference TEs (polymorphic deletions):

```bash
python mega-xtea.py call ... --detect-deletion
```

Output: `absent_MEs.bed` and `absent_MEs_transduction.bed`

### Transduction detection

3' transduction detection uses xTea's simplified algorithm with these key
parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FLANK_WINDOW_SVA` | 5000 | SVA transduction search window (bp) |
| `FLANK_WINDOW_LINE1` | 5000 | LINE1 transduction search window (bp) |
| `MIN_DISC_CUTOFF` | 2 | Minimum discordant reads for transduction call |
| `F_MIN_TRSDCT_DISC_MAP_RATIO` | 0.65 | Minimum mapped-read ratio |

### Large cohort analysis (>1000 samples)

For population-scale studies, use chunked joint calling:

```bash
python mega-xtea.py joint-call \
  --vcf-list all_samples.txt \
  --chunk-size 500 \
  -o population_output/ \
  -t 16
```

## Output Files

| File | Description |
|------|-------------|
| `*.mei.vcf.gz` | MEI insertion calls with genotypes (LINE1, Alu, SVA, HERV-K) |
| `*.absent_MEs.bed` | TE deletion calls (Absent ME) |
| `*.absent_MEs_transduction.bed` | Deletions with 3' transduction evidence |
| `joint_calls.vcf.gz` | Joint-called population VCF across all samples |

VCF ALT allele representations follow standard ME notation:

- `INS:ME:LINE1`
- `INS:ME:ALU`
- `INS:ME:SVA`
- `INS:ME:HERV-K`

## Parameter Reference

Below is a summary of key configurable parameters. All defaults are defined in
`megaxtea/config.py` and can be adjusted programmatically via `MegaXTeaConfig`.

| Category | Parameter | Default | Description |
|----------|-----------|---------|-------------|
| SVA | `REP_SVA_CNS_HEAD` | 400 | SVA consensus 5' head length |
| SVA | `REP_SVA_MIN_LEN` | 300 | Minimum SVA insertion length |
| SVA | `sva_clip_cluster_diff_cutoff` | 200 | Relaxed cluster-diff for SVA VNTR |
| PolyA | `POLYA_RATIO` | 0.4 | Minimum polyA base ratio |
| PolyA | `N_MIN_A_T` | 5 | Minimum consecutive A/T bases |
| ML Genotyping | `n_features` | 15 | Feature vector dimensionality |
| ML Genotyping | `n_estimators` | 20 | Random Forest estimators |
| MEGAnE | `overhang_evalue_threshold` | 1e-5 | BLAST e-value for classification |
| MEGAnE | `gaussian_min_evidence` | 3 | Min reads for Gaussian genotyping |
| Clip/Disc | `MINIMUM_CLIP_MAPQ` | 12 | Minimum MAPQ for clipped reads |
| Clip/Disc | `MINIMUM_DISC_MAPQ` | 20 | Minimum MAPQ for discordant reads |
| Transduction | `FLANK_WINDOW_SVA` | 5000 | SVA transduction search window |
| Transduction | `MIN_DISC_CUTOFF` | 2 | Min discordant reads for transduction |

## Citation

If you use MEGA-xTEA in your research, please cite the original tools:

- **MEGAnE**: Kojima S, Koyama S, Ka M, et al. "Mobile element genomic analysis of
  cohorts (MEGAnE): a robust and efficient method for detection and genotyping of
  mobile element insertions from large cohort short-read resequencing data."
  *Nature Communications*, 2023.
- **xTea**: Chu C, Borges-Monroy R, Viswanadham VV, et al. "Comprehensive
  identification of transposable element insertions using multiple sequencing
  technologies." *Nature Communications*, 2021.

## License

MIT License -- see [LICENSE](LICENSE) for details.

This project incorporates code from MEGAnE (MIT License) and xTea (MIT License).

## Acknowledgments

MEGA-xTEA builds upon the foundational work of the MEGAnE and xTea development
teams. We gratefully acknowledge their contributions to the transposable element
detection community and their release of open-source implementations that made
this integration possible.
