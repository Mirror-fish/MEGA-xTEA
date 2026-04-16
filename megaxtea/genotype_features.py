"""
genotype_features.py -- BAM re-scan for exact ML genotyping features.

Ported from xTea's x_genotype_feature.py (Simon Chu, Harvard DBMI).
For each candidate insertion site, performs a targeted pysam fetch to
compute the exact raw values needed for the 15-dim feature vector.

This replaces the approximation approach in ml_genotype.py by directly
counting reads from the BAM file around each breakpoint.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from multiprocessing import Pool
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants (from xTea global_values.py)
# ---------------------------------------------------------------------------

TSD_CUTOFF = 100             # window for raw clip counting
CLIP_EXACT_CLIP_SLACK = 3    # window for effective (AF) clip counting
DFT_IS = 550                 # default insert size for disc/concordant
DISC_THRESHOLD = 2000        # distance threshold for discordant pairs
BWA_HALF_READ_MIN_SCORE = 45 # for distinguishing left/right dominant full maps
CK_POLYA_CLIP_WIN = 25       # polyA detection window
CK_POLYA_SEQ_MAX = 20        # max polyA sequence length to check
N_MIN_A_T = 5                # min consecutive A/T for polyA


# ---------------------------------------------------------------------------
# Per-site feature result
# ---------------------------------------------------------------------------

@dataclass
class SiteFeatures:
    """Raw feature counts for a single candidate site."""
    chrom: str = ""
    ins_pos: int = 0
    n_af_clip: int = 0        # effective clip reads (max of left/right)
    n_full_map: int = 0       # fully mapped reads spanning breakpoint
    n_raw_lclip: int = 0      # raw left clipped reads within TSD_CUTOFF
    n_raw_rclip: int = 0      # raw right clipped reads within TSD_CUTOFF
    n_disc_pairs: int = 0     # discordant read pairs
    n_concd_pairs: int = 0    # concordant read pairs
    left_coverage: float = 0  # reads covering left breakpoint region
    right_coverage: float = 0 # reads covering right breakpoint region
    n_polyA: int = 0          # polyA reads detected in clipped sequences


# ---------------------------------------------------------------------------
# PolyA helper (simplified from xTea x_polyA)
# ---------------------------------------------------------------------------

def _is_consecutive_polyA_T(seq: str, min_len: int = N_MIN_A_T) -> bool:
    """Check if sequence contains consecutive polyA or polyT."""
    if not seq:
        return False
    count_a = 0
    count_t = 0
    max_a = 0
    max_t = 0
    for c in seq.upper():
        if c == 'A':
            count_a += 1
            max_a = max(max_a, count_a)
            count_t = 0
        elif c == 'T':
            count_t += 1
            max_t = max(max_t, count_t)
            count_a = 0
        else:
            count_a = 0
            count_t = 0
    return max_a >= min_len or max_t >= min_len


# ---------------------------------------------------------------------------
# Core: collect features for one site
# ---------------------------------------------------------------------------

def collect_features_one_site(
    bam_path: str,
    ref_path: str,
    chrom: str,
    ins_pos: int,
    extnd: int = DFT_IS,
) -> SiteFeatures:
    """Collect exact genotyping features for one insertion site.

    Performs a pysam fetch around the insertion position and counts:
    - Raw left/right clipped reads (within TSD_CUTOFF)
    - Effective clip reads (within CLIP_EXACT_CLIP_SLACK)
    - Fully mapped reads spanning the breakpoint
    - Discordant and concordant read pairs
    - Local coverage on left and right sides
    - PolyA reads

    Ported from xTea x_genotype_feature.py::collect_features_one_site().
    """
    import pysam

    sf = SiteFeatures(chrom=chrom, ins_pos=ins_pos)

    start_pos = max(1, ins_pos - extnd)
    end_pos = ins_pos + extnd

    try:
        samfile = pysam.AlignmentFile(bam_path, "rb",
                                      reference_filename=ref_path)
    except Exception as e:
        logger.warning("Cannot open BAM %s: %s", bam_path, e)
        return sf

    # Determine if BAM uses "chr" prefix
    bam_refs = set(samfile.references)
    if chrom in bam_refs:
        chrm_in_bam = chrom
    elif chrom.startswith("chr") and chrom[3:] in bam_refs:
        chrm_in_bam = chrom[3:]
    elif ("chr" + chrom) in bam_refs:
        chrm_in_bam = "chr" + chrom
    else:
        chrm_in_bam = chrom

    n_l_af_clip = 0
    n_r_af_clip = 0
    n_l_full_map = 0
    n_r_full_map = 0
    n_full_map = 0
    n_l_raw_clip = 0
    n_r_raw_clip = 0
    n_disc_pairs = 0
    n_polyA = 0
    left_reads = 0   # reads covering left side of breakpoint
    right_reads = 0  # reads covering right side of breakpoint

    clip_qnames = set()
    concordant_candidates = []

    try:
        for aln in samfile.fetch(chrm_in_bam, start_pos, end_pos):
            if aln.is_duplicate or aln.is_unmapped:
                continue
            cigar = aln.cigartuples
            if not cigar:
                continue

            map_pos = aln.reference_start
            query_seq = aln.query_sequence or ""

            # Mate info
            mate_chrm = "*"
            mate_pos = 0
            if (not aln.mate_is_unmapped and aln.next_reference_id >= 0
                    and aln.next_reference_name):
                mate_chrm = aln.next_reference_name
                mate_pos = aln.next_reference_start

            # --- Coverage counting ---
            # Simple: if read covers left side of ins_pos → left_reads++
            # If covers right side → right_reads++
            read_end = aln.reference_end or (map_pos + 1)
            if map_pos < ins_pos:
                left_reads += 1
            if read_end > ins_pos:
                right_reads += 1

            # --- Fully mapped reads ---
            b_fully_mapped = (len(cigar) == 1 and cigar[0][0] == 0)
            if b_fully_mapped:
                i_map_end = map_pos + cigar[0][1]
                if ins_pos >= map_pos and ins_pos <= i_map_end:
                    n_full_map += 1
                    if abs(map_pos - ins_pos) < BWA_HALF_READ_MIN_SCORE:
                        n_r_full_map += 1
                    elif abs(i_map_end - ins_pos) < BWA_HALF_READ_MIN_SCORE:
                        n_l_full_map += 1

            # --- Left clipped reads ---
            s_clip_seq_ck = ""
            if cigar[0][0] == 4:  # soft-clip at left
                if aln.is_supplementary or aln.is_secondary:
                    continue
                clip_len = cigar[0][1]
                clipped_seq = query_seq[:clip_len] if query_seq else ""

                if abs(map_pos - ins_pos) < TSD_CUTOFF:
                    n_l_raw_clip += 1
                    clip_qnames.add(aln.query_name)

                if abs(map_pos - ins_pos) < CK_POLYA_CLIP_WIN and clipped_seq:
                    if not aln.is_reverse:
                        s_clip_seq_ck = clipped_seq[-CK_POLYA_SEQ_MAX:]
                    else:
                        s_clip_seq_ck = clipped_seq[:CK_POLYA_SEQ_MAX]

                if abs(map_pos - ins_pos) < CLIP_EXACT_CLIP_SLACK:
                    n_l_af_clip += 1

            # --- Right clipped reads ---
            if cigar[-1][0] == 4:  # soft-clip at right
                # Calculate exact clip position
                right_clip_pos = map_pos
                for (op, length) in cigar[:-1]:
                    if op in (4, 5, 1):  # soft-clip, hard-clip, insertion
                        continue
                    right_clip_pos += length

                if aln.is_supplementary or aln.is_secondary:
                    continue

                clip_len = cigar[-1][1]
                clipped_seq = query_seq[-clip_len:] if query_seq else ""

                if abs(right_clip_pos - ins_pos) < TSD_CUTOFF:
                    n_r_raw_clip += 1
                    clip_qnames.add(aln.query_name)

                if abs(right_clip_pos - ins_pos) < CK_POLYA_CLIP_WIN and clipped_seq:
                    if not aln.is_reverse:
                        s_clip_seq_ck = clipped_seq[:CK_POLYA_SEQ_MAX]
                    else:
                        s_clip_seq_ck = clipped_seq[-CK_POLYA_SEQ_MAX:]

                if abs(right_clip_pos - ins_pos) < CLIP_EXACT_CLIP_SLACK:
                    n_r_af_clip += 1

            # --- PolyA detection ---
            if s_clip_seq_ck and _is_consecutive_polyA_T(s_clip_seq_ck):
                n_polyA += 1

            # --- Discordant / concordant ---
            if mate_chrm == "*":
                continue

            if _is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos):
                if abs(map_pos - ins_pos) <= DFT_IS:
                    n_disc_pairs += 1
            else:
                concordant_candidates.append(
                    (aln.query_name, chrm_in_bam, map_pos,
                     mate_chrm, mate_pos)
                )

    except ValueError:
        # Region not found in BAM
        pass

    samfile.close()

    # Count concordant pairs (excluding clip reads, as in xTea)
    n_concd_pairs = 0
    for qname, chrm_r, pos_r, mate_c, mate_p in concordant_candidates:
        if qname in clip_qnames:
            continue
        if _is_concordant(chrm_r, pos_r, mate_c, mate_p, ins_pos):
            n_concd_pairs += 1

    # Select dominant clip side (as xTea does)
    n_af_clip = n_l_af_clip
    adj_full_map = n_full_map - n_l_full_map
    if n_l_af_clip < n_r_af_clip:
        n_af_clip = n_r_af_clip
        adj_full_map = n_full_map - n_r_full_map
    adj_full_map = max(adj_full_map, 0)

    sf.n_af_clip = n_af_clip
    sf.n_full_map = adj_full_map
    sf.n_raw_lclip = n_l_raw_clip
    sf.n_raw_rclip = n_r_raw_clip
    sf.n_disc_pairs = n_disc_pairs
    sf.n_concd_pairs = n_concd_pairs
    sf.left_coverage = max(float(left_reads), 0.0000000001)
    sf.right_coverage = max(float(right_reads), 0.0000000001)
    sf.n_polyA = n_polyA

    return sf


def _is_discordant(chrm: str, pos: int, mate_chrm: str, mate_pos: int) -> bool:
    """Check if a read pair is discordant."""
    if chrm != mate_chrm:
        return True
    if abs(mate_pos - pos) > DISC_THRESHOLD:
        return True
    return False


def _is_concordant(
    chrm: str, pos: int, mate_chrm: str, mate_pos: int,
    ins_pos: int
) -> bool:
    """Check if a read pair is concordant (spanning the insertion)."""
    if chrm != mate_chrm:
        return False
    i_start = min(pos, mate_pos)
    i_end = max(pos, mate_pos)
    if (ins_pos >= i_start and ins_pos <= i_end
            and (ins_pos - i_start) <= DFT_IS
            and (i_end - ins_pos) <= DFT_IS):
        return True
    return False


# ---------------------------------------------------------------------------
# Worker function for multiprocessing
# ---------------------------------------------------------------------------

def _worker(args_tuple):
    """Multiprocessing worker for collect_features_one_site."""
    bam_path, ref_path, chrom, ins_pos = args_tuple
    return collect_features_one_site(bam_path, ref_path, chrom, ins_pos)


# ---------------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------------

def batch_collect_features(
    bed_paths: List[str],
    bam_path: str,
    ref_path: str,
    output_path: str,
    threads: int = 4,
) -> Dict[str, SiteFeatures]:
    """Collect features for all candidates in BED files.

    Args:
        bed_paths: List of MEGAnE BED file paths.
        bam_path: Path to BAM/CRAM file.
        ref_path: Path to reference FASTA.
        output_path: Path to write features TSV.
        threads: Number of parallel workers.

    Returns:
        Dict mapping "chrom:pos" to SiteFeatures.
    """
    import re

    # Collect all candidate positions from BED files
    sites: List[Tuple[str, int, str]] = []  # (chrom, pos, site_key)
    for bed_path in bed_paths:
        if not os.path.isfile(bed_path):
            continue
        with open(bed_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                key = f"{chrom}:{pos}"
                sites.append((chrom, pos, key))

    if not sites:
        logger.info("No candidate sites to scan.")
        return {}

    logger.info("BAM feature scan: %d candidate sites, %d threads",
                len(sites), threads)

    # Build work items
    work_items = [(bam_path, ref_path, chrom, pos) for chrom, pos, _ in sites]

    # Run in parallel
    results: Dict[str, SiteFeatures] = {}
    if threads > 1:
        with Pool(processes=min(threads, len(work_items))) as pool:
            feature_list = pool.map(_worker, work_items)
        for (_, _, key), sf in zip(sites, feature_list):
            results[key] = sf
    else:
        for (chrom, pos, key) in sites:
            sf = collect_features_one_site(bam_path, ref_path, chrom, pos)
            results[key] = sf

    # Write features TSV
    with open(output_path, "w") as fh:
        fh.write("#chrom\tins_pos\tn_af_clip\tn_full_map\t"
                 "n_raw_lclip\tn_raw_rclip\tn_disc_pairs\tn_concd_pairs\t"
                 "left_coverage\tright_coverage\tn_polyA\n")
        for key, sf in results.items():
            fh.write(f"{sf.chrom}\t{sf.ins_pos}\t{sf.n_af_clip}\t"
                     f"{sf.n_full_map}\t{sf.n_raw_lclip}\t{sf.n_raw_rclip}\t"
                     f"{sf.n_disc_pairs}\t{sf.n_concd_pairs}\t"
                     f"{sf.left_coverage:.1f}\t{sf.right_coverage:.1f}\t"
                     f"{sf.n_polyA}\n")

    logger.info("BAM feature scan complete: %d sites processed", len(results))
    return results
