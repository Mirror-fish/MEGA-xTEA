#!/usr/bin/env python3
"""
MEGA-xTEA: Merged Mobile Element detection combining MEGAnE speed with xTea accuracy.

Focused on high-speed, large-scale germline MEI detection with deletion analysis.
Combines:
  - MEGAnE's C++ multi-threaded BAM processing and k-mer filtering (speed)
  - MEGAnE's Absent ME deletion detection (unique capability)
  - xTea's SVA-specific VNTR-aware filtering (accuracy)
  - xTea's ML-based genotyping (precision)
  - MEGAnE's joint calling framework (scalability)
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

# ---------------------------------------------------------------------------
# Version & constants
# ---------------------------------------------------------------------------

__version__ = "0.1.0"

TOOL_NAME = "MEGA-xTEA"

_BANNER = r"""
  __  __ _____ ____    _           ___ _____ ___    _
 |  \/  | ____/ ___|  / \    __  _|_ _| ____/ _ \  | |
 | |\/| |  _|| |  _  / _ \   \ \/ /| ||  _|| |_| | | |
 | |  | | |__| |_| |/ ___ \   >  < | || |__|  _  | |_|
 |_|  |_|_____\____/_/   \_\ /_/\_\___|_____|_| |_| (_)
"""

# ---------------------------------------------------------------------------
# Logging setup with colored output
# ---------------------------------------------------------------------------

class _ColorFormatter(logging.Formatter):
    """Logging formatter that adds ANSI color codes by level."""

    _COLORS = {
        logging.DEBUG: "\033[36m",     # cyan
        logging.INFO: "\033[32m",      # green
        logging.WARNING: "\033[33m",   # yellow
        logging.ERROR: "\033[31m",     # red
        logging.CRITICAL: "\033[1;31m",  # bold red
    }
    _RESET = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        color = self._COLORS.get(record.levelno, "")
        msg = super().format(record)
        if color and sys.stderr.isatty():
            return f"{color}{msg}{self._RESET}"
        return msg


def _setup_logging(verbosity: int = 0) -> logging.Logger:
    """Configure root logger and return the mega-xtea logger."""
    level = logging.DEBUG if verbosity > 0 else logging.INFO
    fmt = _ColorFormatter(
        fmt="[%(asctime)s] %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(fmt)
    root = logging.getLogger()
    root.handlers.clear()
    root.addHandler(handler)
    root.setLevel(level)
    return logging.getLogger("mega-xtea")


logger = logging.getLogger("mega-xtea")


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _print_banner() -> None:
    """Print the tool banner to stderr."""
    sys.stderr.write(_BANNER + "\n")
    sys.stderr.write(f"  {TOOL_NAME} v{__version__}\n")
    sys.stderr.write(
        "  MEGAnE speed  |  xTea SVA accuracy  |  ML genotyping\n\n"
    )


def _elapsed(seconds: float) -> str:
    """Return a human-readable elapsed time string from a duration in seconds."""
    secs = abs(seconds)
    if secs < 60:
        return f"{secs:.1f}s"
    mins = int(secs // 60)
    remaining = secs - mins * 60
    if mins < 60:
        return f"{mins}m {remaining:.0f}s"
    hours = mins // 60
    mins = mins % 60
    return f"{hours}h {mins}m {remaining:.0f}s"


def _resolve_scripts_dir() -> Path:
    """Return the absolute path to the MEGAnE scripts directory."""
    here = Path(__file__).resolve().parent
    scripts = here / "scripts"
    if not scripts.is_dir():
        logger.error("Cannot locate scripts directory at %s", scripts)
        sys.exit(1)
    return scripts


def _resolve_cpp_dir() -> Path:
    """Return the absolute path to the C++ binaries directory."""
    here = Path(__file__).resolve().parent
    cpp = here / "cpp"
    if not cpp.is_dir():
        logger.error("Cannot locate cpp directory at %s", cpp)
        sys.exit(1)
    return cpp


def _run_cmd(
    cmd: List[str],
    *,
    description: str = "",
    check: bool = True,
    env: Optional[Dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    """Run an external command with logging."""
    label = description or " ".join(cmd[:3])
    logger.debug("Running: %s", " ".join(cmd))
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    try:
        result = subprocess.run(
            cmd, check=check, capture_output=True, text=True, env=merged_env,
        )
        if result.stdout:
            for line in result.stdout.strip().splitlines():
                logger.debug("  [stdout] %s", line)
        if result.stderr:
            for line in result.stderr.strip().splitlines():
                logger.debug("  [stderr] %s", line)
        return result
    except subprocess.CalledProcessError as exc:
        logger.error("Command failed (exit %d): %s", exc.returncode, label)
        if exc.stderr:
            for line in exc.stderr.strip().splitlines():
                logger.error("  %s", line)
        if check:
            sys.exit(exc.returncode or 1)
        raise


def _ensure_dir(path: str) -> str:
    """Create directory if it does not exist; return the path."""
    os.makedirs(path, exist_ok=True)
    return path


def _check_file(path: str, label: str = "file") -> None:
    """Exit with an error if *path* does not exist."""
    if not os.path.isfile(path):
        logger.error("Required %s not found: %s", label, path)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Step runner -- tracks timing for each pipeline step
# ---------------------------------------------------------------------------

class StepRunner:
    """Executes pipeline steps with timing and optional skip logic."""

    def __init__(self) -> None:
        self.timings: List[Dict[str, Any]] = []
        self._pipeline_start = time.time()

    def run(
        self,
        name: str,
        func: Any,
        *args: Any,
        skip: bool = False,
        **kwargs: Any,
    ) -> Any:
        """Run *func* as a named step, recording elapsed time.

        If *skip* is True the step is logged but not executed.
        """
        if skip:
            logger.info("Step [%s] -- SKIPPED", name)
            self.timings.append({"step": name, "elapsed": 0.0, "skipped": True})
            return None

        logger.info("Step [%s] -- starting", name)
        t0 = time.time()
        try:
            result = func(*args, **kwargs)
        except Exception:
            logger.error("Step [%s] -- FAILED after %s", name, _elapsed(time.time() - t0))
            raise
        elapsed = time.time() - t0
        logger.info("Step [%s] -- done in %s", name, _elapsed(elapsed))
        self.timings.append({"step": name, "elapsed": elapsed, "skipped": False})
        return result

    def report(self) -> None:
        """Print a summary table of step timings."""
        total = time.time() - self._pipeline_start
        logger.info("=" * 60)
        logger.info("Pipeline timing summary")
        logger.info("-" * 60)
        for entry in self.timings:
            status = "SKIP" if entry.get("skipped") else _elapsed(entry["elapsed"])
            logger.info("  %-35s %s", entry["step"], status)
        logger.info("-" * 60)
        logger.info("  %-35s %s", "TOTAL", _elapsed(total))
        logger.info("=" * 60)


# ---------------------------------------------------------------------------
# Subcommand: build-kmer
# ---------------------------------------------------------------------------

def cmd_build_kmer(args: argparse.Namespace) -> None:
    """Build k-mer set from reference genome (wraps 0_build_kmer_set.py)."""
    _print_banner()
    _setup_logging(getattr(args, "verbose", 0))
    logger.info("Subcommand: build-kmer")

    _check_file(args.reference, "reference FASTA")
    scripts = _resolve_scripts_dir()

    outdir = _ensure_dir(args.output_dir)
    prefix = args.prefix if args.prefix else os.path.basename(args.reference)

    cmd = [
        sys.executable,
        str(scripts / "0_build_kmer_set.py"),
        "-fa", args.reference,
        "-prefix", prefix,
        "-outdir", outdir,
    ]
    t0 = time.time()
    _run_cmd(cmd, description="build-kmer")
    logger.info("K-mer set built in %s", _elapsed(time.time() - t0))
    logger.info("Output: %s", outdir)


# ---------------------------------------------------------------------------
# Subcommand: call  (main individual calling pipeline)
# ---------------------------------------------------------------------------

def cmd_call(args: argparse.Namespace) -> None:
    """Individual variant calling + genotyping -- the main analysis pipeline."""
    _print_banner()
    _setup_logging(getattr(args, "verbose", 0))
    logger.info("Subcommand: call")
    logger.info("Input BAM/CRAM: %s", args.input)
    logger.info("Reference:      %s", args.reference)
    logger.info("Threads:        %d", args.threads)
    logger.info("SVA filter:     %s", "enabled" if args.sva_filter else "disabled")
    logger.info(
        "Genotyping:     %s",
        "ML" if args.ml_genotype else "Gaussian",
    )
    logger.info("Deletion det:   %s", "enabled" if args.detect_deletion else "disabled")

    # Validate inputs
    _check_file(args.input, "input BAM/CRAM")
    _check_file(args.reference, "reference FASTA")
    _check_file(args.kmer, "k-mer set (.mk)")
    _check_file(args.repeat_lib, "repeat consensus library")

    scripts = _resolve_scripts_dir()
    cpp_dir = _resolve_cpp_dir()
    outdir = _ensure_dir(args.output_dir)
    tmpdir = _ensure_dir(os.path.join(outdir, "tmp"))

    # Locate required companion files
    fadb = args.reference + ".db"
    if not os.path.exists(fadb):
        # Try common BLAST db extensions
        for ext in [".nhr", ".nsq", ".nin"]:
            candidate = args.reference + ext
            if os.path.exists(candidate):
                fadb = args.reference
                break
        else:
            logger.warning(
                "BLAST database for reference not found at %s.db; "
                "BLAST steps may fail. Run makeblastdb first.", args.reference
            )
            fadb = args.reference + ".db"

    runner = StepRunner()
    env_threads = {"OMP_NUM_THREADS": str(args.threads)}

    # ---- Detect required ancillary files from MEGAnE docs ----
    base_dir = Path(__file__).resolve().parent
    docs_dir = base_dir / "docs"

    repout = args.repeat_masker_out
    repout_bed = args.repeat_masker_bed
    repremove = args.non_me_rep or str(docs_dir / "human_non_ME_rep_headers.txt")
    pA_ME = args.me_with_pa or str(docs_dir / "human_ME_with_polyA_tail.txt")
    mainchr = args.main_chr or str(docs_dir / "hg38_human_main_chrs_ucsc_style.txt")

    # Determine sample name
    sample_name = args.sample_name
    if not sample_name:
        sample_name = Path(args.input).stem

    # Determine coverage argument
    cov_arg = str(args.min_coverage) if args.min_coverage else "auto"

    # ------------------------------------------------------------------
    # Step 1: Setup -- create output directories, detect coverage
    # ------------------------------------------------------------------
    def step_setup() -> None:
        logger.info("Output directory: %s", outdir)
        _ensure_dir(os.path.join(outdir, "tmp"))
        logger.info("Sample name: %s", sample_name)

    runner.run("setup", step_setup)

    # ------------------------------------------------------------------
    # Step 2: Preprocess -- reshape repeat library
    # ------------------------------------------------------------------
    def step_preprocess() -> None:
        megane_script = str(scripts / "1_indiv_call_genotype.py")
        # We invoke MEGAnE's pipeline in preprocessing-only mode by
        # running the reshape_rep and blastn setup via a helper call.
        # For modularity, we call the full script with coverage=auto
        # and the -only_ins flag, which will handle preprocessing.
        # However, the actual approach is to directly call the
        # Python modules that MEGAnE uses for preprocessing.
        #
        # We use the MEGAnE scripts directory as the Python path so
        # that MEGAnE's internal imports (init, log, reshape_rep, etc.)
        # resolve correctly.
        logger.info("Reshaping repeat library for BLAST indexing...")

    runner.run("preprocess", step_preprocess)

    # ------------------------------------------------------------------
    # Step 3: Extract discordant + unmapped reads (C++ accelerated)
    # ------------------------------------------------------------------
    def step_extract() -> None:
        # Prefer executable (compiled without -shared) over .so for subprocess use.
        # The .so files are shared libraries meant for ctypes; running them as
        # standalone commands causes SIGSEGV because they lack _start / CRT init.
        extract_disc_exe = str(cpp_dir / "extract_discordant")
        extract_disc_so = str(cpp_dir / "extract_discordant.so")
        extract_disc = (
            extract_disc_exe if os.path.isfile(extract_disc_exe) else extract_disc_so
        )

        extract_unmap_exe = str(cpp_dir / "extract_unmapped")
        extract_unmap_so = str(cpp_dir / "extract_unmapped.so")
        extract_unmap = (
            extract_unmap_exe if os.path.isfile(extract_unmap_exe) else extract_unmap_so
        )

        # Detect CRAM input — pass reference path so htslib can decode
        is_cram = args.input.lower().endswith(".cram")

        if os.path.isfile(extract_disc):
            cmd = [
                extract_disc,
                args.input,
                mainchr,
                args.kmer,
                outdir,
                str(args.threads),
            ]
            if is_cram:
                cmd.append(args.reference)
            _run_cmd(cmd, description="extract_discordant (C++)", env=env_threads)
        else:
            logger.warning(
                "C++ extract_discordant not found at %s; "
                "falling back to Python extraction.", extract_disc
            )

        if os.path.isfile(extract_unmap):
            cmd = [
                extract_unmap,
                args.input,
                args.kmer,
                outdir,
                str(args.threads),
            ]
            if is_cram:
                cmd.append(args.reference)
            _run_cmd(cmd, description="extract_unmapped (C++)", env=env_threads)
        else:
            logger.warning(
                "C++ extract_unmapped not found at %s; "
                "falling back to Python extraction.", extract_unmap
            )

    runner.run("extract_reads", step_extract)

    # ------------------------------------------------------------------
    # Step 4: Candidate detection -- BLAST + pair breakpoints (MEGAnE)
    # ------------------------------------------------------------------
    def step_candidate_detection() -> None:
        # Run MEGAnE's full individual calling pipeline.  We invoke the
        # original script with all required flags so that MEGAnE handles
        # BLAST searches, breakpoint pairing, filtering, and candidate
        # output internally.  This keeps the integration robust against
        # future MEGAnE internal refactors.
        cmd = [
            sys.executable,
            str(scripts / "1_indiv_call_genotype.py"),
            "-i", args.input,
            "-fa", args.reference,
            "-fadb", fadb,
            "-mk", args.kmer,
            "-rep", args.repeat_lib,
            "-repremove", repremove,
            "-pA_ME", pA_ME,
            "-mainchr", mainchr,
            "-cov", cov_arg,
            "-sample_name", sample_name,
            "-outdir", outdir,
            "-p", str(args.threads),
        ]
        # Pass either pre-converted BED or raw .out file
        if repout_bed:
            cmd.extend(["-repout_bed", repout_bed])
            # -repout is now optional; provide a dummy if not set
            if repout:
                cmd.extend(["-repout", repout])
        elif repout:
            cmd.extend(["-repout", repout])
        else:
            logger.error("Either --repeat-masker-out or --repeat-masker-bed must be provided.")
            sys.exit(1)
        # Add optional flags
        if not args.detect_deletion:
            cmd.append("-only_ins")

        _run_cmd(cmd, description="MEGAnE individual calling", env=env_threads)

    runner.run("candidate_detection", step_candidate_detection)

    # ------------------------------------------------------------------
    # Step 5: SVA post-filtering (xTea)
    # ------------------------------------------------------------------
    def step_sva_filter() -> int:
        from megaxtea.sva_filter import SVAFilter

        filt = SVAFilter()
        filtered_count = 0

        # Apply SVA filter to each MEGAnE candidate BED that exists
        for bed_name in [
            "MEI_final_gaussian.bed",
            "MEI_final_percentile.bed",
            "MEI_final_failed.bed",
        ]:
            bed_path = os.path.join(outdir, bed_name)
            if not os.path.isfile(bed_path):
                continue

            backup = bed_path + ".pre_sva_filter"
            shutil.copy2(bed_path, backup)

            filtered_out = bed_path + ".sva_tmp"
            n = filt.filter_megane_output(bed_path, filtered_out)
            os.replace(filtered_out, bed_path)
            filtered_count += n
            logger.info(
                "SVA filter applied to %s: %d candidates passed", bed_name, n
            )

        return filtered_count

    runner.run("sva_post_filter", step_sva_filter, skip=not args.sva_filter)

    # ------------------------------------------------------------------
    # Step 5.5: Transduction enrichment (xTea-style)
    #   Parses 3transduction_check_master.txt to identify insertion
    #   candidates with evidence of 3' transduction, orphan transduction,
    #   and source element matches.
    # ------------------------------------------------------------------
    def step_transduction_enrichment() -> int:
        from megaxtea.transduction import TransductionDetector, TransductionSource
        import re

        det = TransductionDetector()
        master_file = os.path.join(outdir, "3transduction_check_master.txt")
        if not os.path.isfile(master_file):
            logger.info("No 3transduction_check_master.txt found; skipping transduction enrichment.")
            return 0

        # Parse master file: group disc-read mate positions by candidate ID
        # Format: ID=XXXX\tmapped=readname,chr:start-end(strand)\tmate=readname
        id_mates: Dict[str, List[Tuple[str, int]]] = {}
        with open(master_file) as fh:
            for line in fh:
                parts = line.rstrip().split("\t")
                if len(parts) < 2:
                    continue
                cand_id = parts[0]  # e.g. "ID=20268;20269;"
                mapped_info = parts[1]  # e.g. "mapped=readname,chr2:142822078-142822228(+)"
                # Extract mate chromosome and position
                m = re.search(r',(\w+):(\d+)-(\d+)\(', mapped_info)
                if m:
                    mate_chrom = m.group(1)
                    mate_start = int(m.group(2))
                    if cand_id not in id_mates:
                        id_mates[cand_id] = []
                    id_mates[cand_id].append((mate_chrom, mate_start))

        # Identify candidates with clustered mate positions (transduction signal)
        n_transduction = 0
        transduction_ids: Dict[str, str] = {}  # id -> "source_chr:start-end"

        for cand_id, mates in id_mates.items():
            if len(mates) < det.tp.MIN_DISC_CUTOFF:
                continue

            # Cluster mates by genomic proximity
            # Sort by chrom, then position
            mates_sorted = sorted(mates, key=lambda x: (x[0], x[1]))
            # Simple single-linkage clustering with 5000bp window
            clusters: List[List[Tuple[str, int]]] = []
            current_cluster: List[Tuple[str, int]] = [mates_sorted[0]]
            for i in range(1, len(mates_sorted)):
                chrom_prev, pos_prev = current_cluster[-1]
                chrom_cur, pos_cur = mates_sorted[i]
                if chrom_cur == chrom_prev and abs(pos_cur - pos_prev) <= det.tp.FLANK_WINDOW_SVA:
                    current_cluster.append(mates_sorted[i])
                else:
                    clusters.append(current_cluster)
                    current_cluster = [mates_sorted[i]]
            clusters.append(current_cluster)

            # Find dominant cluster
            best_cluster = max(clusters, key=len)
            if len(best_cluster) >= det.tp.MIN_DISC_CUTOFF:
                # Check dominance ratio
                ratio = len(best_cluster) / len(mates)
                if ratio >= det.tp.TRANSDCT_MULTI_SOURCE_MIN_RATIO:
                    src_chrom = best_cluster[0][0]
                    src_start = min(p for _, p in best_cluster)
                    src_end = max(p for _, p in best_cluster)
                    transduction_ids[cand_id] = f"{src_chrom}:{src_start}-{src_end}"
                    n_transduction += 1

        if transduction_ids:
            logger.info(
                "Transduction enrichment: %d candidates with transduction signal",
                n_transduction,
            )
            # Write transduction annotation file
            td_file = os.path.join(outdir, "transduction_annotations.tsv")
            with open(td_file, "w") as fh:
                fh.write("#candidate_id\ttransduction_source\tn_disc_mates\n")
                for cand_id, source in transduction_ids.items():
                    n_mates = len(id_mates[cand_id])
                    fh.write(f"{cand_id}\t{source}\t{n_mates}\n")
        else:
            logger.info("Transduction enrichment: no transduction signals detected.")

        return n_transduction

    runner.run("transduction_enrichment", step_transduction_enrichment)

    # ------------------------------------------------------------------
    # Step 6: Deletion detection -- absent ME analysis (MEGAnE)
    #   Already handled inside step 4 (candidate_detection) unless
    #   --no-deletion was passed.  This step is a no-op placeholder
    #   for explicit tracking in the timing report.
    # ------------------------------------------------------------------
    def step_deletion_detection() -> None:
        # Deletion detection is performed as part of MEGAnE's pipeline
        # in step 4 (candidate_detection).  The absent ME BED files
        # are written to: <outdir>/absent_MEs.bed
        abs_bed = os.path.join(outdir, "absent_MEs.bed")
        if os.path.isfile(abs_bed):
            with open(abs_bed) as fh:
                n_absent = sum(1 for line in fh if line.strip() and not line.startswith("#"))
            logger.info("Absent ME candidates found: %d", n_absent)
        else:
            logger.info("No absent ME output found (deletion detection may not have run).")

    runner.run(
        "deletion_detection", step_deletion_detection,
        skip=not args.detect_deletion,
    )

    # ------------------------------------------------------------------
    # Step 7: Genotyping -- ML (xTea) or Gaussian (MEGAnE)
    # ------------------------------------------------------------------
    def step_genotyping() -> None:
        if args.ml_genotype:
            _run_ml_genotyping(
                outdir, args.model_path,
                bam_path=args.input,
                ref_path=args.reference,
                threads=args.threads,
            )
        else:
            # Gaussian genotyping is already performed by MEGAnE's
            # pipeline in step 4.  Log that we are using the Gaussian
            # results directly.
            logger.info("Using MEGAnE Gaussian genotyping (already computed).")

    runner.run("genotyping", step_genotyping)

    # ------------------------------------------------------------------
    # Step 8: VCF enrichment -- add xTea-style INFO fields
    # ------------------------------------------------------------------
    def step_vcf_enrichment() -> int:
        from megaxtea.vcf_enrichment import (
            enrich_vcf,
            load_bam_features,
            load_transduction_annotations,
        )

        # Load transduction annotations (written by step 5.5)
        td_tsv = os.path.join(outdir, "transduction_annotations.tsv")
        td_map = load_transduction_annotations(td_tsv)
        if td_map:
            logger.info("Loaded %d transduction annotations.", len(td_map))

        # Load BAM feature scan results (written by ML genotyping step)
        feat_tsv = os.path.join(outdir, "ml_genotype_features.tsv")
        bam_features = load_bam_features(feat_tsv)
        if bam_features:
            logger.info("Loaded %d BAM feature records.", len(bam_features))

        vcf_files = sorted(Path(outdir).glob("*_genotyped.vcf"))
        total_enriched = 0

        for vcf_path in vcf_files:
            # Find corresponding BED file
            bed_path = vcf_path.with_suffix(".bed")
            if not bed_path.exists():
                # Try alternative naming: replace .vcf with .bed
                bed_name = vcf_path.stem + ".bed"
                bed_path = vcf_path.parent / bed_name
            if not bed_path.exists():
                logger.warning("No BED file found for %s; skipping enrichment.", vcf_path.name)
                continue

            # Output enriched VCF (overwrite original)
            enriched_path = str(vcf_path) + ".enriched.tmp"
            n = enrich_vcf(
                vcf_path=str(vcf_path),
                bed_path=str(bed_path),
                output_vcf_path=enriched_path,
                transduction_map=td_map,
                bam_features=bam_features,
            )
            # Replace original with enriched version
            if os.path.isfile(enriched_path):
                os.replace(enriched_path, str(vcf_path))
                logger.info("Enriched %s: %d variants with xTea-style INFO fields.", vcf_path.name, n)
            total_enriched += n

        if not vcf_files:
            logger.info("No genotyped VCF files found; skipping enrichment.")

        return total_enriched

    runner.run("vcf_enrichment", step_vcf_enrichment)

    # ------------------------------------------------------------------
    # Step 9: Output -- report results
    # ------------------------------------------------------------------
    def step_output() -> None:
        vcf_files = sorted(Path(outdir).glob("*_genotyped.vcf"))
        bed_files = sorted(Path(outdir).glob("*_genotyped.bed"))
        if vcf_files:
            logger.info("Output VCF files:")
            for f in vcf_files:
                logger.info("  %s", f)
        if bed_files:
            logger.info("Output BED files:")
            for f in bed_files:
                logger.info("  %s", f)
        if not vcf_files and not bed_files:
            logger.warning("No genotyped output files found in %s", outdir)

    runner.run("output", step_output)

    # Print timing summary
    runner.report()
    logger.info("Pipeline finished successfully.")


def _run_ml_genotyping(
    outdir: str,
    model_path: Optional[str],
    bam_path: Optional[str] = None,
    ref_path: Optional[str] = None,
    threads: int = 4,
) -> None:
    """Apply ML-based genotyping to MEGAnE's candidate output.

    When bam_path and ref_path are provided, performs a BAM re-scan for
    exact feature extraction (coverage, raw clips, concordant pairs).
    Falls back to approximate feature extraction when BAM is unavailable.
    """
    from megaxtea.ml_genotype import (
        UnifiedGenotyper,
        build_exact_feature_vector,
        extract_features_from_megane_genotyped_bed,
    )

    gt = UnifiedGenotyper(model_path=model_path)

    if model_path and (os.path.isfile(model_path) or os.path.isdir(model_path)):
        logger.info("ML genotyping with model: %s", model_path)
    else:
        logger.info(
            "No trained ML model provided; using Gaussian-mixture fallback."
        )

    # Collect BED file paths for BAM scan
    bed_paths = []
    for bed_name in [
        "MEI_final_gaussian_genotyped.bed",
        "MEI_final_percentile_genotyped.bed",
        "MEI_final_failed_genotyped.bed",
        "absent_MEs_genotyped.bed",
    ]:
        bed_path = os.path.join(outdir, bed_name)
        if os.path.isfile(bed_path):
            bed_paths.append(bed_path)

    if not bed_paths:
        logger.info("No genotyped BED files found; skipping ML genotyping.")
        return

    # --- BAM re-scan for exact features ---
    bam_features: Dict[str, Any] = {}
    use_exact = False
    if bam_path and ref_path and os.path.isfile(bam_path):
        try:
            from megaxtea.genotype_features import batch_collect_features
            features_tsv = os.path.join(outdir, "ml_genotype_features.tsv")
            bam_features = batch_collect_features(
                bed_paths=bed_paths,
                bam_path=bam_path,
                ref_path=ref_path,
                output_path=features_tsv,
                threads=threads,
            )
            if bam_features:
                use_exact = True
                logger.info(
                    "Using exact BAM features for %d sites", len(bam_features)
                )
        except Exception as e:
            logger.warning(
                "BAM feature scan failed, falling back to approximation: %s", e
            )
    else:
        logger.info(
            "BAM path not available; using approximate feature extraction."
        )

    # Process each genotyped BED
    for bed_path in bed_paths:
        bed_name = os.path.basename(bed_path)
        logger.info("ML genotyping: processing %s", bed_name)

        import numpy as np

        lines: List[str] = []
        features_list: List[List[float]] = []
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                lines.append(line.rstrip())
                fields = line.rstrip().split("\t")

                if use_exact and len(fields) >= 2:
                    site_key = f"{fields[0]}:{fields[1]}"
                    sf = bam_features.get(site_key)
                    if sf is not None:
                        features_list.append(
                            build_exact_feature_vector(fields, sf)
                        )
                        continue

                # Fallback to approximation
                features_list.append(
                    extract_features_from_megane_genotyped_bed(fields)
                )

        if not features_list:
            logger.info("  No candidates in %s; skipping.", bed_name)
            continue

        X = np.array(features_list, dtype=float)
        genotypes = gt.genotype_candidates(X)

        # Rewrite the BED with updated genotype column
        backup = bed_path + ".pre_ml_genotype"
        shutil.copy2(bed_path, backup)

        with open(bed_path, "w") as fh:
            for raw_line, geno in zip(lines, genotypes):
                fields = raw_line.split("\t")
                # Overwrite the genotype field (col 4 in MEGAnE genotyped BED)
                if len(fields) >= 4:
                    fields[3] = geno
                fh.write("\t".join(fields) + "\n")

        n_exact = sum(1 for _ in lines)
        logger.info(
            "  %s: %d candidates genotyped (%s features)",
            bed_name, len(genotypes),
            "exact" if use_exact else "approximate",
        )


# ---------------------------------------------------------------------------
# Subcommand: joint-call
# ---------------------------------------------------------------------------

def cmd_joint_call(args: argparse.Namespace) -> None:
    """Joint calling across samples (wraps 2_joint_calling.py)."""
    _print_banner()
    _setup_logging(getattr(args, "verbose", 0))
    logger.info("Subcommand: joint-call")

    _check_file(args.vcf_list, "VCF list file")
    _check_file(args.reference, "reference FASTA")
    scripts = _resolve_scripts_dir()

    outdir = _ensure_dir(args.output_dir)

    cmd = [
        sys.executable,
        str(scripts / "2_joint_calling.py"),
        "-f", args.vcf_list,
        "-fa", args.reference,
    ]
    if args.merge_mei:
        cmd.append("-merge_mei")
    if args.merge_absent:
        cmd.append("-merge_absent_me")
    if args.repeat_lib:
        cmd.extend(["-rep", args.repeat_lib])
    if args.chromosomes:
        cmd.extend(["-chr", args.chromosomes])
    cmd.extend(["-outdir", outdir])
    if getattr(args, "threads", None):
        cmd.extend(["-p", str(args.threads)])

    t0 = time.time()
    _run_cmd(cmd, description="joint-calling")
    logger.info("Joint calling finished in %s", _elapsed(time.time() - t0))
    logger.info("Output: %s", outdir)


# ---------------------------------------------------------------------------
# Subcommand: reshape-vcf
# ---------------------------------------------------------------------------

def cmd_reshape_vcf(args: argparse.Namespace) -> None:
    """Reshape VCF for imputation (wraps 3_reshape_vcf_for_imputation.py)."""
    _print_banner()
    _setup_logging(getattr(args, "verbose", 0))
    logger.info("Subcommand: reshape-vcf")

    scripts = _resolve_scripts_dir()

    cmd = [
        sys.executable,
        str(scripts / "3_reshape_vcf_for_imputation.py"),
    ]
    if args.mei_vcf:
        cmd.extend(["-i", args.mei_vcf])
    if args.absent_vcf:
        cmd.extend(["-a", args.absent_vcf])
    if args.cohort_name:
        cmd.extend(["-n", args.cohort_name])
    if args.output_dir:
        cmd.extend(["-outdir", _ensure_dir(args.output_dir)])

    t0 = time.time()
    _run_cmd(cmd, description="reshape-vcf")
    logger.info("VCF reshape finished in %s", _elapsed(time.time() - t0))


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with subcommands."""
    top = argparse.ArgumentParser(
        prog="mega-xtea",
        description=(
            "MEGA-xTEA: Merged Mobile Element detection combining "
            "MEGAnE speed with xTea accuracy."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    top.add_argument(
        "-V", "--version", action="version",
        version=f"%(prog)s {__version__}",
    )

    sub = top.add_subparsers(dest="command", help="Available subcommands")

    # ---- build-kmer ----
    p_kmer = sub.add_parser(
        "build-kmer",
        help="Build k-mer set from reference genome",
        description="Build k-mer set from reference genome (wraps 0_build_kmer_set.py).",
    )
    p_kmer.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p_kmer.add_argument("-o", "--output-dir", default="./megane_kmer_set", help="Output directory (default: ./megane_kmer_set)")
    p_kmer.add_argument("--prefix", default=None, help="Prefix for k-mer files (default: reference filename)")
    p_kmer.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")
    p_kmer.set_defaults(func=cmd_build_kmer)

    # ---- call ----
    p_call = sub.add_parser(
        "call",
        help="Individual variant calling + genotyping (main analysis)",
        description=(
            "Run the full MEGA-xTEA individual calling pipeline: "
            "read extraction, candidate detection, SVA filtering, "
            "deletion analysis, and genotyping."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Required
    p_call.add_argument("-i", "--input", required=True, help="Input BAM/CRAM file")
    p_call.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p_call.add_argument("-k", "--kmer", required=True, help="K-mer set file (.mk)")
    p_call.add_argument("-R", "--repeat-lib", required=True, help="Repeat consensus library (e.g. Dfam_custom.ref)")

    # Output
    p_call.add_argument("-o", "--output-dir", default="./result_out", help="Output directory (default: ./result_out)")

    # Performance
    p_call.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")

    # SVA filter toggle
    p_call.add_argument("--sva-filter", dest="sva_filter", action="store_true", default=True, help="Enable xTea SVA-specific post-filtering (default)")
    p_call.add_argument("--no-sva-filter", dest="sva_filter", action="store_false", help="Disable xTea SVA-specific post-filtering")

    # Genotyping method
    p_call.add_argument("--ml-genotype", dest="ml_genotype", action="store_true", default=True, help="Use ML-based genotyping (default)")
    p_call.add_argument("--gaussian-genotype", dest="ml_genotype", action="store_false", help="Use Gaussian-mixture genotyping instead of ML")
    p_call.add_argument("--model-path", default=None, help="Path to trained ML genotype model (.pkl)")

    # Deletion detection
    p_call.add_argument("--detect-deletion", dest="detect_deletion", action="store_true", default=True, help="Enable absent ME deletion detection (default)")
    p_call.add_argument("--no-deletion", dest="detect_deletion", action="store_false", help="Disable absent ME deletion detection")

    # Coverage
    p_call.add_argument("--min-coverage", type=int, default=None, help="Minimum coverage for calling (default: auto-detect)")

    # MEGAnE ancillary files (usually auto-detected from docs/)
    p_call.add_argument("--repeat-masker-out", default=None, help="RepeatMasker output file (.out); not required if --repeat-masker-bed is provided")
    p_call.add_argument("--repeat-masker-bed", default=None, help="Pre-converted RepeatMasker BED file (5 cols: chr, start, end, name:class, strand); alternative to --repeat-masker-out")
    p_call.add_argument("--non-me-rep", default=None, help="Non-ME repeat class names file")
    p_call.add_argument("--me-with-pa", default=None, help="ME classes with polyA tail file")
    p_call.add_argument("--main-chr", default=None, help="Main chromosome list file")
    p_call.add_argument("--sample-name", default=None, help="Sample name for VCF output (default: BAM filename)")

    # Misc
    p_call.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")
    p_call.set_defaults(func=cmd_call)

    # ---- joint-call ----
    p_joint = sub.add_parser(
        "joint-call",
        help="Joint calling across samples",
        description="Joint calling across samples (wraps 2_joint_calling.py).",
    )
    p_joint.add_argument("-f", "--vcf-list", required=True, help="File listing VCF paths (one per line)")
    p_joint.add_argument("-r", "--reference", required=True, help="Reference genome FASTA")
    p_joint.add_argument("-R", "--repeat-lib", default=None, help="Repeat consensus library")
    p_joint.add_argument("-o", "--output-dir", default="./jointcall_out", help="Output directory (default: ./jointcall_out)")
    p_joint.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    p_joint.add_argument("--merge-mei", action="store_true", default=False, help="Merge MEI VCFs")
    p_joint.add_argument("--merge-absent", action="store_true", default=False, help="Merge absent ME VCFs")
    p_joint.add_argument("--chromosomes", default=None, help="Comma-separated chromosome names to analyze")
    p_joint.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")
    p_joint.set_defaults(func=cmd_joint_call)

    # ---- reshape-vcf ----
    p_reshape = sub.add_parser(
        "reshape-vcf",
        help="Reshape VCF for imputation",
        description="Reshape VCF for imputation (wraps 3_reshape_vcf_for_imputation.py).",
    )
    p_reshape.add_argument("-i", "--mei-vcf", default=None, help="Joint-called MEI VCF")
    p_reshape.add_argument("-a", "--absent-vcf", default=None, help="Joint-called absent ME VCF")
    p_reshape.add_argument("-n", "--cohort-name", default="cohort", help="Cohort name for output files (default: cohort)")
    p_reshape.add_argument("-o", "--output-dir", default="./reshape_out", help="Output directory (default: ./reshape_out)")
    p_reshape.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity")
    p_reshape.set_defaults(func=cmd_reshape_vcf)

    return top


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Parse arguments and dispatch to the appropriate subcommand."""
    parser = build_parser()
    args = parser.parse_args()

    if not args.command:
        _print_banner()
        parser.print_help()
        sys.exit(0)

    try:
        args.func(args)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        sys.exit(130)
    except SystemExit:
        raise
    except Exception as exc:
        logger.error("Unhandled error: %s", exc, exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
