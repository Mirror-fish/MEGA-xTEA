"""
fp_filter.py -- Phase 1 False-Positive reduction filters (xTea-style).

Implements three key FP reduction mechanisms from xTea's post-filtering:

  1. AFConflictFilter  -- 4-ratio allele-frequency quality check
  2. Low-divergence reference TE copy filter -- RepeatMasker divergence
  3. Cluster consistency check (simplified) -- breakpoint position validation

These filters run AFTER ML genotyping and BEFORE VCF enrichment,
consuming ml_genotype_features.tsv + genotyped BED + RepeatMasker data.
"""

from __future__ import annotations

import logging
import os
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FilterDecision:
    """Per-variant filter decision."""
    chrom: str = ""
    pos: int = 0
    var_id: str = ""
    filters: List[str] = field(default_factory=list)  # filter tags to ADD

    @property
    def is_filtered(self) -> bool:
        return len(self.filters) > 0


# ---------------------------------------------------------------------------
# 1. AF Conflict Filter (xTea x_post_filter.py:1039-1185)
# ---------------------------------------------------------------------------

# [PHASE1_TUNABLE] AF conflict thresholds -- match xTea germline defaults.
# All 4 ratios must exceed their threshold (AND logic).
AF_EF_CLIP_CUTOFF = 0.075
AF_EF_DISC_CUTOFF = 0.075
AF_CLIP_FULL_CUTOFF = 0.075
AF_DISC_CONCORD_CUTOFF = 0.075


def af_conflict_check(
    left_chimeric: int,
    right_chimeric: int,
    left_hybrid: int,
    right_hybrid: int,
    n_af_clip: int,
    n_full_map: int,
    n_disc_pairs: int,
    n_concd_pairs: int,
) -> Tuple[bool, str]:
    """Check whether a candidate passes the 4-ratio AF quality filter.

    Returns (pass, reason).  pass=True means the candidate is OK.

    Ratio 1: (left+right clip-to-consensus) / total effective clip
    Ratio 2: (left+right disc-to-consensus) / total disc pairs
    Ratio 3: effective clip / (effective clip + fully mapped)
    Ratio 4: disc pairs / (disc pairs + concordant pairs)
    """
    # Ratio 1: consensus-aligned clips / effective clips
    n_cns_clip = left_chimeric + right_chimeric
    f_ef_clip = float(n_cns_clip) / float(n_af_clip) if n_af_clip > 0 else 0.0

    # Ratio 2: consensus-aligned disc / total disc
    n_cns_disc = left_hybrid + right_hybrid
    f_ef_disc = float(n_cns_disc) / float(n_disc_pairs) if n_disc_pairs > 0 else 0.0

    # Ratio 3: clip / (clip + fully-mapped)  -- AF frequency indicator
    denom3 = n_af_clip + n_full_map
    f_clip_full = float(n_af_clip) / float(denom3) if denom3 > 0 else 0.0

    # Ratio 4: disc / (disc + concordant)
    denom4 = n_disc_pairs + n_concd_pairs
    f_disc_concd = float(n_disc_pairs) / float(denom4) if denom4 > 0 else 0.0

    # All 4 must pass (AND logic, matching xTea)
    b1 = f_ef_clip > AF_EF_CLIP_CUTOFF
    b2 = f_ef_disc > AF_EF_DISC_CUTOFF
    b3 = f_clip_full > AF_CLIP_FULL_CUTOFF
    b4 = f_disc_concd > AF_DISC_CONCORD_CUTOFF

    if b1 and b2 and b3 and b4:
        return True, ""

    parts = []
    if not b1:
        parts.append("ef_clip=%.3f" % f_ef_clip)
    if not b2:
        parts.append("ef_disc=%.3f" % f_ef_disc)
    if not b3:
        parts.append("clip_full=%.3f" % f_clip_full)
    if not b4:
        parts.append("disc_concd=%.3f" % f_disc_concd)
    return False, "AF_CONFLICT(%s)" % ",".join(parts)


# ---------------------------------------------------------------------------
# 2. Low-divergence reference TE copy filter
#    (xTea x_post_filter.py:379-407 + x_annotation.py)
# ---------------------------------------------------------------------------

# [PHASE1_TUNABLE] Divergence cutoffs from xTea global_values.py
REP_DIVERGENT_CUTOFF = 15.0       # default for most TE types
REP_LOW_DIVERGENT_CUTOFF = 7.0    # stricter for LINE1


@dataclass
class RepeatAnnotation:
    """A single RepeatMasker annotation entry with divergence."""
    chrom: str
    start: int
    end: int
    name: str       # e.g. "L1HS"
    family: str     # e.g. "LINE/L1"
    strand: str
    divergence: float  # % mismatch from RepeatMasker .out column[1]


class RepeatMaskerIndex:
    """Indexed RepeatMasker annotations for fast overlap + divergence queries.

    Uses sorted-start-position arrays with bisect for O(log n) lookups,
    avoiding the need for a full IntervalTree dependency.
    """

    def __init__(self) -> None:
        # {chrom: list of RepeatAnnotation sorted by start}
        self._data: Dict[str, List[RepeatAnnotation]] = {}
        self._starts: Dict[str, List[int]] = {}  # parallel sorted start arrays

    def add(self, ann: RepeatAnnotation) -> None:
        self._data.setdefault(ann.chrom, []).append(ann)

    def build_index(self) -> None:
        """Sort by start position and build bisect arrays."""
        for chrom in self._data:
            self._data[chrom].sort(key=lambda a: a.start)
            self._starts[chrom] = [a.start for a in self._data[chrom]]

    def query_overlaps(
        self, chrom: str, pos: int, window: int = 100
    ) -> List[RepeatAnnotation]:
        """Find all annotations overlapping pos ± window."""
        if chrom not in self._starts:
            return []
        starts = self._starts[chrom]
        annotations = self._data[chrom]
        # Find entries whose start <= pos + window
        right_idx = bisect_right(starts, pos + window)
        # Scan backwards to find entries that overlap
        results = []
        for i in range(max(0, bisect_left(starts, pos - 50000)), right_idx):
            ann = annotations[i]
            if ann.end >= pos - window and ann.start <= pos + window:
                results.append(ann)
        return results

    def __len__(self) -> int:
        return sum(len(v) for v in self._data.values())


def load_repeatmasker_out(out_path: str, boundary_extnd: int = 100) -> RepeatMaskerIndex:
    """Load RepeatMasker .out file and build indexed annotation.

    Parses divergence from column[1] of the standard RepeatMasker output.
    Extends boundaries by `boundary_extnd` bp on each side (like xTea).
    """
    index = RepeatMaskerIndex()
    non_rep_headers = {"SW", "bit", "score"}

    if not os.path.isfile(out_path):
        logger.warning("RepeatMasker .out file not found: %s", out_path)
        return index

    with open(out_path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 15:
                continue
            if fields[0] in non_rep_headers:
                continue
            try:
                divergence = float(fields[1])
                chrom = fields[4]
                start = int(fields[5]) - 1 - boundary_extnd  # 0-based + extend
                end = int(fields[6]) + boundary_extnd
                name = fields[9]        # subfamily e.g. "L1HS"
                family = fields[10]     # family e.g. "LINE/L1"
                strand = "+" if fields[8] == "+" else "-"
                if start < 0:
                    start = 0
                index.add(RepeatAnnotation(chrom, start, end, name, family, strand, divergence))
            except (ValueError, IndexError):
                continue

    index.build_index()
    logger.info("Loaded %d RepeatMasker annotations with divergence from %s", len(index), out_path)
    return index


def _te_type_from_bed(bed_fields: List[str]) -> str:
    """Extract TE type label (e.g. 'LINE', 'SINE', 'SVA') from BED col3."""
    if len(bed_fields) > 3:
        return bed_fields[3].upper()
    return ""


def _te_family_matches(cand_type: str, ref_family: str) -> bool:
    """Check if candidate TE type matches reference annotation family.

    cand_type: from BED e.g. "LINE", "SINE", "SVA", "LTR"
    ref_family: from RepeatMasker e.g. "LINE/L1", "SINE/Alu", "Retroposon/SVA"
    """
    ref_upper = ref_family.upper()
    cand_upper = cand_type.upper()

    # Direct match
    if cand_upper in ref_upper:
        return True

    # Common mappings
    mapping = {
        "LINE": "LINE",
        "SINE": "SINE",
        "SVA": "SVA",
        "LTR": "LTR",
        "HERV": "LTR",
        "ALU": "SINE",
    }
    mapped = mapping.get(cand_upper, cand_upper)
    return mapped in ref_upper


def low_div_ref_te_check(
    rmsk_index: RepeatMaskerIndex,
    chrom: str,
    pos: int,
    cand_te_type: str,
) -> Tuple[bool, float, str]:
    """Check if candidate falls within a low-divergence same-type reference TE.

    Returns (is_filtered, min_divergence, reason).
    """
    if len(rmsk_index) == 0:
        return False, -1.0, ""

    # Determine divergence cutoff based on TE type
    cand_upper = cand_te_type.upper()
    if "LINE" in cand_upper or "L1" in cand_upper:
        cutoff = REP_LOW_DIVERGENT_CUTOFF  # 7% for LINE1
    else:
        cutoff = REP_DIVERGENT_CUTOFF  # 15% for others

    overlaps = rmsk_index.query_overlaps(chrom, pos, window=100)
    if not overlaps:
        return False, -1.0, ""

    # Find lowest divergence among same-type overlapping copies
    min_div = -1.0
    for ann in overlaps:
        if not _te_family_matches(cand_te_type, ann.family):
            continue
        if min_div < 0 or ann.divergence < min_div:
            min_div = ann.divergence

    if min_div < 0:
        return False, -1.0, ""

    if min_div <= cutoff:
        return True, min_div, "LOW_DIV_REF_TE(div=%.1f%%,cutoff=%.0f%%)" % (min_div, cutoff)

    return False, min_div, ""


# ---------------------------------------------------------------------------
# 3. Cluster consistency check (simplified -- Phase 1)
#    Checks breakpoint position consistency from BED col4/col5
# ---------------------------------------------------------------------------

# [PHASE1_TUNABLE] Maximum allowed gap between left and right breakpoint
# positions for a "consistent" two-side candidate.
MAX_CONSISTENT_BP_GAP = 500  # bp


def cluster_consistency_check(
    left_ref_pos: int,
    right_ref_pos: int,
    support_type: str,
) -> Tuple[bool, str]:
    """Simplified cluster consistency: check breakpoint position gap.

    For two-side candidates, the gap between left and right breakpoints
    should be within a reasonable TSD range.  Large gaps suggest
    inconsistent clustering (likely FP).

    Returns (pass, reason).
    """
    if "two_side" not in support_type.lower():
        return True, ""  # Only check two-side candidates

    if left_ref_pos <= 0 or right_ref_pos <= 0:
        return True, ""  # Missing position data, skip

    gap = abs(left_ref_pos - right_ref_pos)
    if gap > MAX_CONSISTENT_BP_GAP:
        return False, "CLUSTER_GAP(%dbp)" % gap

    return True, ""


# ---------------------------------------------------------------------------
# Main filtering orchestrator
# ---------------------------------------------------------------------------

def run_fp_filters(
    bed_path: str,
    bam_features: Dict[str, Dict[str, Any]],
    rmsk_index: Optional[RepeatMaskerIndex] = None,
) -> Dict[str, FilterDecision]:
    """Run all Phase 1 FP filters on genotyped BED candidates.

    Args:
        bed_path: Path to genotyped BED (11-column MEGAnE format).
        bam_features: {chrom:pos: {field: value}} from ml_genotype_features.tsv.
        rmsk_index: Optional RepeatMaskerIndex for divergence filtering.

    Returns:
        {variant_id: FilterDecision} for variants that should be filtered.
    """
    decisions: Dict[str, FilterDecision] = {}

    if not os.path.isfile(bed_path):
        return decisions

    n_af_filtered = 0
    n_div_filtered = 0
    n_cluster_filtered = 0
    n_total = 0

    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue

            n_total += 1
            chrom = fields[0]
            pos_0based = int(fields[1])
            cand_te_type = fields[3]
            var_id = fields[-1] if fields[-1].startswith("ID=") else ""

            decision = FilterDecision(chrom=chrom, pos=pos_0based, var_id=var_id)

            # --- Parse BED evidence for AF check ---
            left_chimeric = 0
            right_chimeric = 0
            left_hybrid = 0
            right_hybrid = 0
            left_ref_pos = 0
            right_ref_pos = 0

            if len(fields) > 4:
                m = re.search(r"chimeric=(\d+)", fields[4])
                if m:
                    left_chimeric = int(m.group(1))
                m = re.search(r"hybrid=(\d+)", fields[4])
                if m:
                    left_hybrid = int(m.group(1))
                m = re.search(r"ref_pos=(\d+)", fields[4])
                if m:
                    left_ref_pos = int(m.group(1))

            if len(fields) > 5:
                m = re.search(r"chimeric=(\d+)", fields[5])
                if m:
                    right_chimeric = int(m.group(1))
                m = re.search(r"hybrid=(\d+)", fields[5])
                if m:
                    right_hybrid = int(m.group(1))
                m = re.search(r"ref_pos=(\d+)", fields[5])
                if m:
                    right_ref_pos = int(m.group(1))

            # --- Get BAM features for this site ---
            pos_key = "%s:%d" % (chrom, pos_0based)
            bam_f = bam_features.get(pos_key, {})

            n_af_clip = int(bam_f.get("n_af_clip", 0))
            n_full_map = int(bam_f.get("n_full_map", 0))
            n_disc_pairs = int(bam_f.get("n_disc_pairs", 0))
            n_concd_pairs = int(bam_f.get("n_concd_pairs", 0))

            # --- Filter 1: AF Conflict ---
            if bam_f:  # Only if BAM features available
                af_pass, af_reason = af_conflict_check(
                    left_chimeric, right_chimeric,
                    left_hybrid, right_hybrid,
                    n_af_clip, n_full_map,
                    n_disc_pairs, n_concd_pairs,
                )
                if not af_pass:
                    decision.filters.append("AF_CONFLICT")
                    n_af_filtered += 1

            # --- Filter 2: Low-divergence reference TE copy ---
            if rmsk_index is not None and len(rmsk_index) > 0:
                div_filtered, min_div, div_reason = low_div_ref_te_check(
                    rmsk_index, chrom, pos_0based, cand_te_type,
                )
                if div_filtered:
                    decision.filters.append("LOW_DIV_REF")
                    n_div_filtered += 1

            # --- Filter 3: Cluster consistency (simplified) ---
            # Determine support_type from BED col8 (subfamily_pred)
            support_type = ""
            if len(fields) > 8:
                pred = fields[8]
                # Infer from evidence pattern
                if left_chimeric > 0 and right_chimeric > 0:
                    support_type = "two_side"

            if support_type:
                cl_pass, cl_reason = cluster_consistency_check(
                    left_ref_pos, right_ref_pos, support_type,
                )
                if not cl_pass:
                    decision.filters.append("CLUSTER_INCONSIST")
                    n_cluster_filtered += 1

            if decision.is_filtered:
                decisions[var_id or pos_key] = decision

    logger.info(
        "FP filter results: %d/%d candidates flagged "
        "(AF_CONFLICT=%d, LOW_DIV_REF=%d, CLUSTER_INCONSIST=%d)",
        len(decisions), n_total, n_af_filtered, n_div_filtered, n_cluster_filtered,
    )
    return decisions


def apply_fp_filters_to_vcf(
    vcf_path: str,
    output_vcf_path: str,
    decisions: Dict[str, FilterDecision],
) -> Tuple[int, int]:
    """Apply FP filter decisions to a VCF file.

    Adds filter tags to FILTER column. Candidates with FP filters
    are moved from PASS to the specific filter tag.

    Returns (n_filtered, n_total).
    """
    if not os.path.isfile(vcf_path):
        return 0, 0

    new_lines = []
    n_filtered = 0
    n_total = 0
    header_added = False

    # Extra FILTER header lines
    extra_headers = [
        '##FILTER=<ID=AF_CONFLICT,Description="Failed xTea-style AF quality check (4-ratio test, all must be >7.5%)">\n',
        '##FILTER=<ID=LOW_DIV_REF,Description="Falls within low-divergence same-type reference TE copy">\n',
        '##FILTER=<ID=CLUSTER_INCONSIST,Description="Breakpoint cluster positions inconsistent (gap >%dbp)">\n' % MAX_CONSISTENT_BP_GAP,
    ]

    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                new_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                if not header_added:
                    for eh in extra_headers:
                        new_lines.append(eh)
                    header_added = True
                new_lines.append(line)
                continue

            # Variant line
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                new_lines.append(line)
                continue

            n_total += 1
            var_id = parts[2]
            pos_1based = parts[1]
            chrom = parts[0]

            # Find matching decision
            decision = None
            # Try by ID
            id_nums = re.findall(r"\d+", var_id)
            for key, dec in decisions.items():
                key_nums = re.findall(r"\d+", key)
                if id_nums and key_nums and id_nums[0] == key_nums[0]:
                    decision = dec
                    break
            # Fallback by position
            if decision is None:
                pos_key = "%s:%d" % (chrom, int(pos_1based) - 1)
                decision = decisions.get(pos_key)

            if decision is not None and decision.is_filtered:
                # Add filter tags
                current_filter = parts[6]
                new_tags = ";".join(decision.filters)
                if current_filter == "PASS" or current_filter == ".":
                    parts[6] = new_tags
                else:
                    parts[6] = current_filter + ";" + new_tags
                n_filtered += 1

            new_lines.append("\t".join(parts) + "\n")

    with open(output_vcf_path, "w") as fh:
        fh.writelines(new_lines)

    return n_filtered, n_total
