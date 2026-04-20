version 1.0

## ============================================================================
## MEGA-xTEA Individual Calling Workflow
##
## Scatter over multiple samples, running Step 1 (individual variant calling
## + genotyping) for each BAM independently. Designed for local Cromwell or
## miniwdl execution with a pre-built Docker image.
## ============================================================================

# ---------------------------------------------------------------------------
# Task: Run mega-xtea call on a single sample
# ---------------------------------------------------------------------------
task mega_xtea_call {
    input {
        # ---- Required inputs ----
        File input_bam                       # .bam or .cram
        File input_bam_index              # .bai/.csi or .crai
        String sample_name

        # ---- Reference bundle (shared across all samples) ----
        File reference_fasta
        File reference_fasta_index        # .fai
        File reference_blast_nhr          # makeblastdb outputs
        File reference_blast_nin
        File reference_blast_nsq
        File? repeat_masker_out           # RepeatMasker .out file (optional if repeat_masker_bed provided)
        File? repeat_masker_bed           # Pre-converted BED (alternative to .out)
        File repeat_library               # Dfam FASTA (e.g. Dfam3.2_human.rep)

        # ---- K-mer set (pre-built via build-kmer) ----
        File kmer_mk                      # .mk binary
        File kmer_mi                      # .mi index

        # ---- Optional companion files (bundled in image under /opt/mega-xtea/docs/) ----
        File? non_me_rep_headers          # human_non_ME_rep_headers.txt
        File? me_with_polya              # human_ME_with_polyA_tail.txt
        File? main_chr_list               # e.g. hg38_human_main_chrs_ucsc_style.txt

        # ---- Analysis options ----
        Boolean sva_filter = true
        Boolean ml_genotype = true
        Boolean detect_deletion = true
        Int sva_breakpoint_gap = 150      # SVA breakpoint pairing gap (bp). Try 100-200.
        File? ml_model_pkl              # Pre-trained sklearn .pkl model (optional)
        File? ml_model_tar              # Deep Forest model directory as .tar.gz (optional, preferred)

        # ---- Runtime ----
        Int threads = 8
        Int memory_gb = 16
        Int disk_gb = 100
        String docker_image = "mega-xtea:latest"
    }

    # Derive the reference prefix for localization
    String ref_basename = basename(reference_fasta)

    command <<<
        set -euo pipefail

        # ----------------------------------------------------------------
        # Stage reference files into a single directory so BLAST can find
        # the database files alongside the FASTA.
        # ----------------------------------------------------------------
        REF_DIR="$(pwd)/reference"
        mkdir -p "${REF_DIR}"
        ln -s ~{reference_fasta}       "${REF_DIR}/~{ref_basename}"
        ln -s ~{reference_fasta_index} "${REF_DIR}/~{ref_basename}.fai"
        ln -s ~{reference_blast_nhr}   "${REF_DIR}/~{ref_basename}.nhr"
        ln -s ~{reference_blast_nin}   "${REF_DIR}/~{ref_basename}.nin"
        ln -s ~{reference_blast_nsq}   "${REF_DIR}/~{ref_basename}.nsq"
        
        if [ "~{defined(repeat_masker_out)}" == "true" ]; then
            ln -s ~{repeat_masker_out} "${REF_DIR}/~{ref_basename}.out"
        fi

        # Stage BAM/CRAM + index together (auto-detect format)
        BAM_DIR="$(pwd)/bam"
        mkdir -p "${BAM_DIR}"
        INPUT_FILE="~{input_bam}"
        if [[ "${INPUT_FILE}" == *.cram ]]; then
            ln -s ~{input_bam}       "${BAM_DIR}/~{sample_name}.cram"
            ln -s ~{input_bam_index} "${BAM_DIR}/~{sample_name}.cram.crai"
            INPUT_PATH="${BAM_DIR}/~{sample_name}.cram"
        else
            ln -s ~{input_bam}       "${BAM_DIR}/~{sample_name}.bam"
            ln -s ~{input_bam_index} "${BAM_DIR}/~{sample_name}.bam.bai"
            INPUT_PATH="${BAM_DIR}/~{sample_name}.bam"
        fi

        # Stage k-mer files together
        KMER_DIR="$(pwd)/kmer"
        mkdir -p "${KMER_DIR}"
        KMER_BASE=$(basename ~{kmer_mk} .mk)
        ln -s ~{kmer_mk} "${KMER_DIR}/${KMER_BASE}.mk"
        ln -s ~{kmer_mi} "${KMER_DIR}/${KMER_BASE}.mi"

        # Output directory
        OUT_DIR="$(pwd)/output/~{sample_name}"
        mkdir -p "${OUT_DIR}"

        # ----------------------------------------------------------------
        # Build command
        # ----------------------------------------------------------------
        CMD="python3 /opt/mega-xtea/mega-xtea.py call"
        CMD="${CMD} -i ${INPUT_PATH}"
        CMD="${CMD} -r ${REF_DIR}/~{ref_basename}"
        CMD="${CMD} -k ${KMER_DIR}/${KMER_BASE}.mk"
        CMD="${CMD} -R ~{repeat_library}"
        CMD="${CMD} -o ${OUT_DIR}"
        CMD="${CMD} -t ~{threads}"
        CMD="${CMD} --sample-name ~{sample_name}"

        # Pass repeat masker data (BED preferred over .out)
        if [ "~{defined(repeat_masker_bed)}" == "true" ]; then
            CMD="${CMD} --repeat-masker-bed ~{repeat_masker_bed}"
        else
            CMD="${CMD} --repeat-masker-out ${REF_DIR}/~{ref_basename}.out"
        fi

        # SVA filter
        if [ "~{sva_filter}" = "true" ]; then
            CMD="${CMD} --sva-filter"
        else
            CMD="${CMD} --no-sva-filter"
        fi
        CMD="${CMD} --sva-breakpoint-gap ~{sva_breakpoint_gap}"

        # Genotyping method
        if [ "~{ml_genotype}" = "true" ]; then
            CMD="${CMD} --ml-genotype"
        else
            CMD="${CMD} --gaussian-genotype"
        fi

        # Deletion detection
        if [ "~{detect_deletion}" = "true" ]; then
            CMD="${CMD} --detect-deletion"
        else
            CMD="${CMD} --no-deletion"
        fi

        # Optional companion files
        if ~{if defined(non_me_rep_headers) then "true" else "false"} = "true"; then
            CMD="${CMD} --non-me-rep ~{select_first([non_me_rep_headers, 'NONE'])}"
        fi
        if ~{if defined(me_with_polya) then "true" else "false"} = "true"; then
            CMD="${CMD} --me-with-pa ~{select_first([me_with_polya, 'NONE'])}"
        fi
        if ~{if defined(main_chr_list) then "true" else "false"} = "true"; then
            CMD="${CMD} --main-chr ~{select_first([main_chr_list, 'NONE'])}"
        fi
        if ~{if defined(ml_model_tar) then "true" else "false"} = "true"; then
            # Extract Deep Forest model directory from tar.gz
            MODEL_DIR="$(pwd)/ml_model"
            mkdir -p "${MODEL_DIR}"
            tar -xzf ~{select_first([ml_model_tar, 'NONE'])} -C "${MODEL_DIR}"
            # Find the directory containing param.pkl (may be nested)
            MODEL_PATH=$(dirname $(find "${MODEL_DIR}" -name "param.pkl" -type f | head -1))
            CMD="${CMD} --model-path ${MODEL_PATH}"
        elif ~{if defined(ml_model_pkl) then "true" else "false"} = "true"; then
            CMD="${CMD} --model-path ~{select_first([ml_model_pkl, 'NONE'])}"
        fi

        CMD="${CMD} -v"

        echo "============================================"
        echo "MEGA-xTEA call: ~{sample_name}"
        echo "Command: ${CMD}"
        echo "============================================"

        eval ${CMD}

        # ----------------------------------------------------------------
        # Collect outputs — tar the full result directory for portability
        # ----------------------------------------------------------------
        tar -czf "~{sample_name}.results.tar.gz" -C "$(pwd)/output" "~{sample_name}"
    >>>

    output {
        File results_tar = "~{sample_name}.results.tar.gz"
        Array[File] vcf_files = glob("output/~{sample_name}/*.vcf*")
        Array[File] bed_files = glob("output/~{sample_name}/*.bed")
        Array[File] tsv_files = glob("output/~{sample_name}/*.tsv")
        File stdout_log = stdout()
        File stderr_log = stderr()
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        maxRetries: 1
    }

    meta {
        description: "Run MEGA-xTEA individual calling on a single BAM sample"
        author: "MEGA-xTEA"
    }
}

# ---------------------------------------------------------------------------
# Workflow: Scatter over sample array
# ---------------------------------------------------------------------------
workflow mega_xtea_batch {
    input {
        # ---- Sample list ----
        Array[File] input_bams
        Array[File] input_bam_indices
        Array[String] sample_names

        # ---- Shared reference bundle ----
        File reference_fasta
        File reference_fasta_index
        File reference_blast_nhr
        File reference_blast_nin
        File reference_blast_nsq
        File? repeat_masker_out
        File? repeat_masker_bed
        File repeat_library

        # ---- K-mer set ----
        File kmer_mk
        File kmer_mi

        # ---- Optional companion files ----
        File? non_me_rep_headers
        File? me_with_polya
        File? main_chr_list

        # ---- Options ----
        Boolean sva_filter = true
        Boolean ml_genotype = true
        Boolean detect_deletion = true
        Int sva_breakpoint_gap = 150
        File? ml_model_pkl
        File? ml_model_tar

        # ---- Runtime per sample ----
        Int threads_per_sample = 8
        Int memory_gb_per_sample = 16
        Int disk_gb_per_sample = 100
        String docker_image = "mega-xtea:latest"
    }

    # Scatter over all samples in parallel
    scatter (idx in range(length(input_bams))) {
        call mega_xtea_call {
            input:
                input_bam           = input_bams[idx],
                input_bam_index     = input_bam_indices[idx],
                sample_name         = sample_names[idx],
                reference_fasta     = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_blast_nhr = reference_blast_nhr,
                reference_blast_nin = reference_blast_nin,
                reference_blast_nsq = reference_blast_nsq,
                repeat_masker_out   = repeat_masker_out,
                repeat_masker_bed   = repeat_masker_bed,
                repeat_library      = repeat_library,
                kmer_mk             = kmer_mk,
                kmer_mi             = kmer_mi,
                non_me_rep_headers  = non_me_rep_headers,
                me_with_polya       = me_with_polya,
                main_chr_list       = main_chr_list,
                sva_filter          = sva_filter,
                ml_genotype         = ml_genotype,
                detect_deletion     = detect_deletion,
                sva_breakpoint_gap  = sva_breakpoint_gap,
                ml_model_pkl        = ml_model_pkl,
                ml_model_tar        = ml_model_tar,
                threads             = threads_per_sample,
                memory_gb           = memory_gb_per_sample,
                disk_gb             = disk_gb_per_sample,
                docker_image        = docker_image
        }
    }

    output {
        Array[File] all_results_tar = mega_xtea_call.results_tar
        Array[Array[File]] all_vcf_files = mega_xtea_call.vcf_files
        Array[Array[File]] all_bed_files = mega_xtea_call.bed_files
        Array[Array[File]] all_tsv_files = mega_xtea_call.tsv_files
        Array[File] all_stdout_logs = mega_xtea_call.stdout_log
        Array[File] all_stderr_logs = mega_xtea_call.stderr_log
    }

    meta {
        description: "MEGA-xTEA batch individual calling — scatter over multiple samples"
        author: "MEGA-xTEA"
        version: "0.1.0"
    }
}
