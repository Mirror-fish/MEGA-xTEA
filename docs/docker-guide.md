# MEGA-xTEA Docker Guide

## Quick Reference

```bash
# Build
docker build -t mega-xtea:latest .

# Run (basic)
docker run --rm mega-xtea:latest --help

# Run (full analysis)
docker run --rm \
  -v /path/to/bams:/data/input \
  -v /path/to/reference:/data/reference \
  -v /path/to/kmer:/data/kmer \
  -v /path/to/output:/data/output \
  mega-xtea:latest call \
  -i /data/input/sample.bam \
  -r /data/reference/GRCh38.fa \
  -k /data/kmer/kmer.mk \
  -R /data/reference/Dfam_custom.ref \
  -o /data/output \
  -t 8
```

---

## 1. Build the Docker Image

### From the repository

```bash
git clone https://github.com/Mirror-fish/MEGA-xTEA.git
cd MEGA-xTEA
docker build -t mega-xtea:latest .
```

The multi-stage build will:
1. Compile htslib 1.19 from source
2. Compile all 5 C++ shared libraries (`.so`)
3. Create a slim runtime image (~800 MB) with all dependencies

### Build options

```bash
# Specify platform (for M1/M2 Mac building x86 images)
docker build --platform linux/amd64 -t mega-xtea:latest .

# With build cache for faster rebuilds
docker build --build-arg BUILDKIT_INLINE_CACHE=1 -t mega-xtea:latest .

# Tag with version
docker build -t mega-xtea:0.1.0 -t mega-xtea:latest .
```

### Verify the build

```bash
docker run --rm mega-xtea:latest --version
# Expected: mega-xtea 0.1.0

docker run --rm mega-xtea:latest call --help
```

---

## 2. Directory Mount Convention

The container expects data through volume mounts at these paths:

| Host path | Container path | Purpose |
|-----------|---------------|---------|
| BAM/CRAM files | `/data/input` | Input sequencing data |
| Reference genome | `/data/reference` | FASTA + index files |
| K-mer set | `/data/kmer` | Pre-built `.mk` file |
| Results | `/data/output` | All output files |

**Important**: Mount the reference genome directory (not just the file), because samtools/BLAST need access to `.fai`, `.dict`, and index files alongside the FASTA.

---

## 3. Complete Workflow

### Step 0: Build k-mer set (one-time per reference genome)

```bash
docker run --rm \
  -v /path/to/reference:/data/reference \
  -v /path/to/kmer:/data/kmer \
  mega-xtea:latest build-kmer \
  -r /data/reference/GRCh38.fa \
  -o /data/kmer \
  -t 4
```

This produces the `.mk` k-mer set file. Only needs to run once per reference genome. Takes ~10 min and ~48 GB RAM for GRCh38.

### Step 1: Individual calling

```bash
docker run --rm \
  -v /home/user/bams:/data/input:ro \
  -v /home/user/ref:/data/reference:ro \
  -v /home/user/kmer:/data/kmer:ro \
  -v /home/user/results/sample1:/data/output \
  mega-xtea:latest call \
  -i /data/input/sample1.bam \
  -r /data/reference/GRCh38.fa \
  -k /data/kmer/kmer.mk \
  -R /data/reference/Dfam_custom.ref \
  -o /data/output \
  -t 8 \
  --sva-filter \
  --ml-genotype \
  --detect-deletion
```

Flags:
- `--sva-filter` — Enable xTea's VNTR-aware SVA filtering (recommended, default ON)
- `--no-sva-filter` — Disable SVA filtering for faster runs
- `--ml-genotype` — Use Random Forest genotyping (default)
- `--gaussian-genotype` — Use Gaussian fitting (faster, less accurate)
- `--detect-deletion` — Enable absent ME deletion detection (default ON)
- `--no-deletion` — Skip deletion detection
- `-t 8` — Thread count (match your available cores)

### Step 2: Joint calling (multi-sample)

First, create a text file listing all individual VCF paths (container paths):

```bash
# On host: create vcf_list.txt
echo "/data/input/sample1/sample1.mei.vcf.gz" > vcf_list.txt
echo "/data/input/sample2/sample2.mei.vcf.gz" >> vcf_list.txt
# ... add all samples
```

Then run:

```bash
docker run --rm \
  -v /home/user/results:/data/input:ro \
  -v /home/user/joint_output:/data/output \
  -v /home/user/vcf_list.txt:/data/vcf_list.txt:ro \
  mega-xtea:latest joint-call \
  --vcf-list /data/vcf_list.txt \
  -o /data/output \
  -t 16
```

### Step 3: Reshape VCF for imputation (optional)

```bash
docker run --rm \
  -v /home/user/joint_output:/data/input:ro \
  -v /home/user/imputation:/data/output \
  mega-xtea:latest reshape-vcf \
  -i /data/input/joint_calls.vcf.gz \
  -o /data/output/reshaped.vcf.gz
```

---

## 4. Batch Processing

### Using a shell loop

```bash
REF_DIR=/home/user/ref
KMER_DIR=/home/user/kmer
BAM_DIR=/home/user/bams
OUT_DIR=/home/user/results

for bam in ${BAM_DIR}/*.bam; do
  sample=$(basename "$bam" .bam)
  mkdir -p "${OUT_DIR}/${sample}"

  docker run --rm \
    -v ${BAM_DIR}:/data/input:ro \
    -v ${REF_DIR}:/data/reference:ro \
    -v ${KMER_DIR}:/data/kmer:ro \
    -v ${OUT_DIR}/${sample}:/data/output \
    mega-xtea:latest call \
    -i "/data/input/${sample}.bam" \
    -r /data/reference/GRCh38.fa \
    -k /data/kmer/kmer.mk \
    -R /data/reference/Dfam_custom.ref \
    -o /data/output \
    -t 4 \
    --sample-name "${sample}"
done
```

### Using GNU parallel (recommended for HPC)

```bash
# Create sample list
ls /home/user/bams/*.bam | sed 's|.*/||;s|\.bam$||' > samples.txt

# Run with 4 concurrent containers, 4 threads each
cat samples.txt | parallel -j 4 \
  docker run --rm \
    -v /home/user/bams:/data/input:ro \
    -v /home/user/ref:/data/reference:ro \
    -v /home/user/kmer:/data/kmer:ro \
    -v /home/user/results/{}:/data/output \
    mega-xtea:latest call \
    -i /data/input/{}.bam \
    -r /data/reference/GRCh38.fa \
    -k /data/kmer/kmer.mk \
    -R /data/reference/Dfam_custom.ref \
    -o /data/output \
    -t 4 \
    --sample-name {}
```

### Singularity (for HPC clusters without Docker)

```bash
# Convert Docker image to Singularity
singularity build mega-xtea.sif docker://mirror-fish/mega-xtea:latest

# Or from local Docker image
singularity build mega-xtea.sif docker-daemon://mega-xtea:latest

# Run with Singularity
singularity exec \
  --bind /path/to/bams:/data/input \
  --bind /path/to/ref:/data/reference \
  --bind /path/to/kmer:/data/kmer \
  --bind /path/to/output:/data/output \
  mega-xtea.sif \
  python3 /opt/mega-xtea/mega-xtea.py call \
  -i /data/input/sample.bam \
  -r /data/reference/GRCh38.fa \
  -k /data/kmer/kmer.mk \
  -R /data/reference/Dfam_custom.ref \
  -o /data/output \
  -t 8
```

---

## 5. Resource Recommendations

| Cohort size | Threads/sample | RAM/sample | Disk/sample | Strategy |
|-------------|---------------|------------|-------------|----------|
| 1-10 | 8 | 16 GB | 20 GB | Sequential |
| 10-100 | 4 | 8 GB | 20 GB | 4 parallel containers |
| 100-1000 | 4 | 8 GB | 20 GB | GNU parallel or Snakemake |
| 1000-10000+ | 4 | 8 GB | 20 GB | HPC scheduler (SLURM/LSF) |

### Memory notes
- K-mer set building (Step 0): requires ~48 GB RAM for GRCh38
- Individual calling (Step 1): ~8 GB for 30x WGS
- Joint calling (Step 2): scales with sample count, use `--chunk-size` for large cohorts

---

## 6. Snakemake Integration

Example `Snakefile` for cluster deployment:

```python
SAMPLES = glob_wildcards("bams/{sample}.bam").sample

rule all:
    input: "joint_output/joint_calls.vcf.gz"

rule individual_call:
    input: bam="bams/{sample}.bam"
    output: vcf="results/{sample}/{sample}.mei.vcf.gz"
    threads: 4
    resources: mem_mb=8000
    singularity: "mega-xtea.sif"
    shell:
        """
        python3 /opt/mega-xtea/mega-xtea.py call \
          -i {input.bam} \
          -r reference/GRCh38.fa \
          -k kmer/kmer.mk \
          -R reference/Dfam_custom.ref \
          -o results/{wildcards.sample} \
          -t {threads} \
          --sample-name {wildcards.sample}
        """

rule joint_call:
    input: expand("results/{sample}/{sample}.mei.vcf.gz", sample=SAMPLES)
    output: "joint_output/joint_calls.vcf.gz"
    threads: 16
    resources: mem_mb=32000
    singularity: "mega-xtea.sif"
    shell:
        """
        ls results/*/*.mei.vcf.gz > vcf_list.txt
        python3 /opt/mega-xtea/mega-xtea.py joint-call \
          --vcf-list vcf_list.txt \
          -o joint_output \
          -t {threads}
        """
```

Run with SLURM:

```bash
snakemake --use-singularity --jobs 100 \
  --cluster "sbatch -p normal -t 2:00:00 -c {threads} --mem={resources.mem_mb}"
```

---

## 7. Troubleshooting

### "permission denied" on output directory
```bash
# Ensure the output directory is writable
chmod 777 /path/to/output
# Or run container as current user
docker run --rm --user $(id -u):$(id -g) ...
```

### "libhts.so not found"
This shouldn't happen with the multi-stage build. If it does:
```bash
docker run --rm mega-xtea:latest bash -c "ldconfig -p | grep hts"
```

### Out of memory during k-mer building
The `build-kmer` step needs ~48 GB for GRCh38. Ensure Docker has enough memory:
```bash
# Docker Desktop: Settings → Resources → Memory → 50 GB+
# Or use --memory flag
docker run --rm --memory=50g ...
```

### Container exits without output
Check logs:
```bash
docker run mega-xtea:latest call ... 2>&1 | tee mega-xtea.log
```

### Using CRAM instead of BAM
Mount the reference genome with read access — CRAM decoding requires it:
```bash
-v /path/to/reference:/data/reference:ro
```
