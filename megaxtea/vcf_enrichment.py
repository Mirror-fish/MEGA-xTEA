"""
vcf_enrichment.py -- Post-process VCF to add xTea-style INFO fields.

After MEGAnE generates VCFs and MEGA-xTEA runs SVA filtering + transduction
enrichment, this module enriches VCF files with additional INFO fields to
match xTea's output format, including:

  - LCLIP / RCLIP / LDISC / RDISC (evidence counts)
  - LPOLYA / RPOLYA (polyA read counts)
  - TD_SRC (transduction source)
  - SUBTYPE (insertion subtype classification)
  - STRAND (insertion orientation)
  - TSDLEN (TSD length, derived from HOMLEN)
  - REF_REP (repeat annotation overlap)
  - Evidence counts from BAM feature scan (if available)
"""

from __future__ import annotations

import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Load transduction annotations
# ---------------------------------------------------------------------------

def load_transduction_annotations(tsv_path: str) -> Dict[str, str]:
    """Load transduction_annotations.tsv → {candidate_id: source_string}.

    The candidate_id in transduction_annotations.tsv looks like
    "ID=20268;20269;" which matches the ID field in BED/VCF.
    We normalize it to just the numeric part for matching.
    """
    td_map: Dict[str, str] = {}
    if not os.path.isfile(tsv_path):
        return td_map

    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            cand_id = parts[0]  # e.g. "ID=20268;20269;"
            source = parts[1]   # e.g. "chr2:142822078-142822228"
            # Normalize ID: extract numeric parts
            nums = re.findall(r'\d+', cand_id)
            for n in nums:
                td_map[n] = source
    return td_map


# ---------------------------------------------------------------------------
# Load BAM features (from ml_genotype_features.tsv)
# ---------------------------------------------------------------------------

def load_bam_features(tsv_path: str) -> Dict[str, Dict[str, Any]]:
    """Load BAM feature scan results → {chrom:pos: {field: value}}."""
    features: Dict[str, Dict[str, Any]] = {}
    if not os.path.isfile(tsv_path):
        return features

    with open(tsv_path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                continue
            if not header:
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            key = f"{parts[0]}:{parts[1]}"
            row = {}
            for i, name in enumerate(header):
                if i < len(parts):
                    try:
                        row[name] = float(parts[i])
                    except ValueError:
                        row[name] = parts[i]
            features[key] = row
    return features


# ---------------------------------------------------------------------------
# Parse evidence from BED fields
# ---------------------------------------------------------------------------

def _parse_bed_evidence(bed_fields: List[str]) -> Dict[str, str]:
    """Extract evidence counts from MEGAnE BED fields for VCF INFO.

    BED format (11 columns):
      0: chrom
      1: start
      2: end
      3: TE class
      4: MEI_left:ref_pos=X,chimeric=Y,hybrid=Z,pA=N
      5: MEI_right:ref_pos=X,chimeric=Y,hybrid=Z,pA=N
      6: confidence:high/low
      7: unique:yes/no,...
      8: subfamily_pred:status=...,MEI=SVA_E,740/764,+/+
      9: 3transduction:yes/no
     10: ID=XXX
    """
    info = {}

    # Parse left breakpoint (col4)
    if len(bed_fields) > 4:
        col4 = bed_fields[4]
        m = re.search(r'chimeric=(\d+)', col4)
        info['LCLIP'] = m.group(1) if m else '0'
        m = re.search(r'hybrid=(\d+)', col4)
        info['LDISC'] = m.group(1) if m else '0'
        m = re.search(r'pA=(\d+)', col4)
        info['LPOLYA'] = m.group(1) if m else '0'

    # Parse right breakpoint (col5)
    if len(bed_fields) > 5:
        col5 = bed_fields[5]
        m = re.search(r'chimeric=(\d+)', col5)
        info['RCLIP'] = m.group(1) if m else '0'
        m = re.search(r'hybrid=(\d+)', col5)
        info['RDISC'] = m.group(1) if m else '0'
        m = re.search(r'pA=(\d+)', col5)
        info['RPOLYA'] = m.group(1) if m else '0'

    # Parse TSD length from breakpoint positions
    if len(bed_fields) > 5:
        try:
            left_pos = int(re.search(r'ref_pos=(\d+)', bed_fields[4]).group(1))
            right_pos = int(re.search(r'ref_pos=(\d+)', bed_fields[5]).group(1))
            tsd_len = left_pos - right_pos  # positive = TSD, negative = deletion
            info['TSDLEN'] = str(abs(tsd_len)) if tsd_len > 0 else '0'
        except (AttributeError, ValueError):
            info['TSDLEN'] = '0'

    # Parse strand from col8 (subfamily_pred)
    if len(bed_fields) > 8:
        pred = bed_fields[8]
        # Pattern: +/+ or -/- or +/- etc.
        m = re.search(r'[,/]([+-])/([+-])', pred)
        if m:
            l_strand, r_strand = m.group(1), m.group(2)
            if l_strand == r_strand:
                info['STRAND'] = l_strand
            else:
                info['STRAND'] = '.'
        elif re.search(r',([+-])$', pred):
            m = re.search(r',([+-])$', pred)
            info['STRAND'] = m.group(1)
        else:
            info['STRAND'] = '.'

    # Parse transduction flag from col9
    if len(bed_fields) > 9:
        info['IS_TRANSDUCTION'] = 'yes' in bed_fields[9].lower()

    # Determine SUBTYPE from evidence pattern
    lclip = int(info.get('LCLIP', '0'))
    rclip = int(info.get('RCLIP', '0'))
    ldisc = int(info.get('LDISC', '0'))
    rdisc = int(info.get('RDISC', '0'))

    has_left_clip = lclip > 0
    has_right_clip = rclip > 0
    has_left_disc = ldisc > 0
    has_right_disc = rdisc > 0

    is_td = info.get('IS_TRANSDUCTION', False)

    if has_left_clip and has_right_clip:
        info['SUBTYPE'] = 'two_side_tprt_both'
    elif (has_left_clip and has_right_disc) or (has_right_clip and has_left_disc):
        if is_td:
            info['SUBTYPE'] = 'one_side_and_half_transduction'
        else:
            info['SUBTYPE'] = 'one_half_side_tprt_both'
    elif has_left_clip or has_right_clip:
        if is_td:
            info['SUBTYPE'] = 'one_side_and_half_transduction'
        else:
            info['SUBTYPE'] = 'one_half_side_tprt'
    else:
        info['SUBTYPE'] = 'disc_only'

    return info


# ---------------------------------------------------------------------------
# Main VCF enrichment function
# ---------------------------------------------------------------------------

def enrich_vcf(
    vcf_path: str,
    bed_path: str,
    output_vcf_path: str,
    transduction_map: Optional[Dict[str, str]] = None,
    bam_features: Optional[Dict[str, Dict[str, Any]]] = None,
) -> int:
    """Enrich a VCF file with additional INFO fields from BED and annotations.

    Args:
        vcf_path: Input VCF (from MEGAnE output_genotyped_vcf.py).
        bed_path: Corresponding genotyped BED file.
        output_vcf_path: Output enriched VCF path.
        transduction_map: {numeric_id: source_string} from transduction detection.
        bam_features: {chrom:pos: {field: value}} from BAM feature scan.

    Returns:
        Number of variants enriched.
    """
    if transduction_map is None:
        transduction_map = {}
    if bam_features is None:
        bam_features = {}

    # Load BED data keyed by ID
    bed_by_id: Dict[str, List[str]] = {}
    bed_by_pos: Dict[str, List[str]] = {}
    if os.path.isfile(bed_path):
        with open(bed_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) > 10:
                    # ID is last column (col10)
                    cand_id = fields[-1]  # e.g. "ID=1234"
                    bed_by_id[cand_id] = fields
                    # Also key by position for fallback matching
                    pos_key = f"{fields[0]}:{fields[1]}"
                    bed_by_pos[pos_key] = fields

    # Process VCF
    new_header_lines = []
    new_variant_lines = []
    n_enriched = 0

    # Additional INFO header definitions
    extra_headers = [
        '##INFO=<ID=LCLIP,Number=1,Type=Integer,Description="Number of left clipped reads (aligned to copies)">',
        '##INFO=<ID=RCLIP,Number=1,Type=Integer,Description="Number of right clipped reads (aligned to copies)">',
        '##INFO=<ID=LDISC,Number=1,Type=Integer,Description="Number of left discordant reads (aligned to copies)">',
        '##INFO=<ID=RDISC,Number=1,Type=Integer,Description="Number of right discordant reads (aligned to copies)">',
        '##INFO=<ID=LPOLYA,Number=1,Type=Integer,Description="Number of left polyA reads">',
        '##INFO=<ID=RPOLYA,Number=1,Type=Integer,Description="Number of right polyA reads">',
        '##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description="TSD length">',
        '##INFO=<ID=SUBTYPE,Number=1,Type=String,Description="Subtype based on evidence pattern">',
        '##INFO=<ID=TD_SRC,Number=1,Type=String,Description="Transduction source">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion orientation (+/- for sense/antisense)">',
        '##INFO=<ID=LCOV,Number=1,Type=Float,Description="Left focal coverage">',
        '##INFO=<ID=RCOV,Number=1,Type=Float,Description="Right focal coverage">',
        '##INFO=<ID=AF_CLIP,Number=1,Type=Integer,Description="Number of effective clip reads">',
        '##INFO=<ID=AF_FMAP,Number=1,Type=Integer,Description="Number of effective fully mapped reads">',
        '##INFO=<ID=AF_DISC,Number=1,Type=Integer,Description="Number of effective discordant pairs">',
        '##INFO=<ID=AF_CONCORDNT,Number=1,Type=Integer,Description="Number of effective concordant pairs">',
        '##INFO=<ID=LRAWCLIP,Number=1,Type=Integer,Description="Number of left raw clip reads">',
        '##INFO=<ID=RRAWCLIP,Number=1,Type=Integer,Description="Number of right raw clip reads">',
    ]

    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                new_header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                # Insert extra headers before #CHROM line
                for eh in extra_headers:
                    new_header_lines.append(eh + "\n")
                new_header_lines.append(line)
                continue

            # Variant line
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                new_variant_lines.append(line)
                continue

            chrom = parts[0]
            pos_1based = parts[1]
            var_id = parts[2]
            info_field = parts[7]

            # Find corresponding BED record
            bed_fields = None
            # Try by ID
            for bid, bf in bed_by_id.items():
                id_nums = re.findall(r'\d+', bid)
                var_id_nums = re.findall(r'\d+', var_id)
                if id_nums and var_id_nums and id_nums[0] == var_id_nums[0]:
                    bed_fields = bf
                    break
            # Fallback: by position
            if bed_fields is None:
                pos_0based = str(int(pos_1based) - 1)
                pos_key = f"{chrom}:{pos_0based}"
                bed_fields = bed_by_pos.get(pos_key)

            if bed_fields is None:
                new_variant_lines.append(line)
                continue

            # Extract evidence from BED
            evidence = _parse_bed_evidence(bed_fields)

            # Add transduction source
            cand_id_str = bed_fields[-1] if bed_fields else ""
            id_nums = re.findall(r'\d+', cand_id_str)
            td_src = "not_transduction"
            for n in id_nums:
                if n in transduction_map:
                    td_src = transduction_map[n]
                    break

            # Add BAM features if available
            pos_key = f"{chrom}:{int(pos_1based) - 1}"
            bam_f = bam_features.get(pos_key, {})

            # Build additional INFO fields
            extra_info_parts = []
            extra_info_parts.append(f"LCLIP={evidence.get('LCLIP', '0')}")
            extra_info_parts.append(f"RCLIP={evidence.get('RCLIP', '0')}")
            extra_info_parts.append(f"LDISC={evidence.get('LDISC', '0')}")
            extra_info_parts.append(f"RDISC={evidence.get('RDISC', '0')}")
            extra_info_parts.append(f"LPOLYA={evidence.get('LPOLYA', '0')}")
            extra_info_parts.append(f"RPOLYA={evidence.get('RPOLYA', '0')}")
            extra_info_parts.append(f"TSDLEN={evidence.get('TSDLEN', '0')}")
            extra_info_parts.append(f"SUBTYPE={evidence.get('SUBTYPE', 'unknown')}")
            extra_info_parts.append(f"TD_SRC={td_src}")
            extra_info_parts.append(f"STRAND={evidence.get('STRAND', '.')}")

            # BAM-derived fields
            if bam_f:
                lcov = bam_f.get('left_coverage', 0)
                rcov = bam_f.get('right_coverage', 0)
                extra_info_parts.append(f"LCOV={lcov:.1f}")
                extra_info_parts.append(f"RCOV={rcov:.1f}")
                extra_info_parts.append(f"AF_CLIP={int(bam_f.get('n_af_clip', 0))}")
                extra_info_parts.append(f"AF_FMAP={int(bam_f.get('n_full_map', 0))}")
                extra_info_parts.append(f"AF_DISC={int(bam_f.get('n_disc_pairs', 0))}")
                extra_info_parts.append(f"AF_CONCORDNT={int(bam_f.get('n_concd_pairs', 0))}")
                extra_info_parts.append(f"LRAWCLIP={int(bam_f.get('n_raw_lclip', 0))}")
                extra_info_parts.append(f"RRAWCLIP={int(bam_f.get('n_raw_rclip', 0))}")

            # Append to existing INFO
            enriched_info = info_field + ";" + ";".join(extra_info_parts)
            parts[7] = enriched_info

            new_variant_lines.append("\t".join(parts) + "\n")
            n_enriched += 1

    # Write output
    with open(output_vcf_path, "w") as fh:
        fh.writelines(new_header_lines)
        fh.writelines(new_variant_lines)

    return n_enriched
