"""
sva_filter.py -- SVA-specific post-processing filter.

This is the KEY integration module that ports xTea's SVA detection logic
(x_post_filter.py :: post_processing_SVA, lines 63-123) into a clean
Python module that can be applied as a post-processing step on MEGAnE's
candidate output.

Origin:
  - Filter hierarchy and SVA parameters: xTea (Simon Chu, Harvard)
  - Candidate input format: MEGAnE (Shohei Kojima, RIKEN)

Algorithm overview (算法概述):
  SVA 的 VNTR 区域导致 clip cluster 间距较大, 因此需要放松一致性阈值
  (sva_clip_cluster_diff_cutoff=200 vs 默认300)。过滤按四级层次执行:
    1. two_side    -- 两侧都有 clip 证据 (最可信)
    2. one_half    -- 一侧半强信号
    3. one_side    -- 仅单侧
    4. other       -- 直接丢弃
  每一级都要求 polyA 信号, 并排除 polyA 主导的假阳性。
  落在参考基因组 SVA 拷贝内的候选需额外审查 divergence rate。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from megaxtea.config import DEFAULT_CONFIG, MegaXTeaConfig, SVAParams
from megaxtea.polyA_detector import PolyADetector

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Support type labels (来自 xTea global_values 的插入子类型标签)
# ---------------------------------------------------------------------------

TWO_SIDE = "two_side"
ONE_HALF_SIDE = "one_half"
ONE_SIDE = "one_side"
OTHER = "other"

# Transduction labels
ORPHAN_TRANSDUCTION = "orphan_or_sibling_transduction"
TRANSDUCTION = "transduction"


# ---------------------------------------------------------------------------
# Candidate record helper
# 候选位点记录 -- 统一 MEGAnE 与 xTea 的输出格式
# ---------------------------------------------------------------------------

@dataclass
class SVACandidate:
    """Represents a single SVA insertion candidate.

    Fields are intentionally a superset of what MEGAnE and xTea each produce,
    so the integration function can populate whichever subset is available.
    """

    chrom: str = ""
    pos: int = 0

    # Support type (来自 xTea 的子类型分类)
    support_type: str = OTHER

    # Clip / disc evidence counts
    left_clip_cns: int = 0
    right_clip_cns: int = 0
    left_disc_cns: int = 0
    right_disc_cns: int = 0
    left_polyA: int = 0
    right_polyA: int = 0

    # Consensus hit position (用于判断 polyA 区域)
    left_cns_hit_start: int = 0
    right_cns_hit_start: int = 0

    # Annotation overlap
    in_repeat: bool = False
    repeat_divergence: float = 100.0
    repeat_subfamily: str = ""
    repeat_family: str = ""
    repeat_start: int = 0
    repeat_end: int = 0

    # Flags
    has_polyA: bool = False
    is_orphan_transduction: bool = False
    is_transduction: bool = False

    # Allele-frequency / quality filter pass (from xTea af_filter equivalent)
    af_qualified: bool = True

    # Raw line from the input file (preserved for output)
    raw_line: str = ""

    # Insertion length (from consensus alignment)
    ins_length: int = 0


# ---------------------------------------------------------------------------
# Core SVA filter class
# ---------------------------------------------------------------------------

class SVAFilter:
    """SVA post-processing filter implementing xTea's filter hierarchy.

    Usage::

        filt = SVAFilter()
        passed = filt.filter_candidates(candidates)

    或者直接对 MEGAnE 的输出文件进行过滤::

        filt = SVAFilter()
        filt.filter_megane_output("megane_candidates.tsv", "sva_filtered.tsv")
    """

    def __init__(self, config: Optional[MegaXTeaConfig] = None) -> None:
        self.cfg = config or DEFAULT_CONFIG
        self.sva = self.cfg.sva
        self.polya_det = PolyADetector(self.cfg.polya)

    # ------------------------------------------------------------------
    # polyA dominance checks (对应 xTea xtprt_filter 方法)
    # ------------------------------------------------------------------

    def is_polyA_dominant_two_side_sva(
        self, cand: SVACandidate, cns_head: int = 0
    ) -> bool:
        """Check if a two-side candidate is dominated by polyA on both sides.

        对应 xTea xtprt_filter.is_polyA_dominant_two_side_sva:
        如果两侧的 consensus 比对都落在 polyA 区域 (>= REP_SVA_POLYA_START),
        且没有足够的非 polyA clip 证据, 则判定为 polyA 假阳性。

        Args:
            cand: SVA candidate record.
            cns_head: Minimum consensus head position for non-polyA evidence.
                      Defaults to SVA CNS head (400bp).
        """
        if cns_head == 0:
            cns_head = self.sva.REP_SVA_CNS_HEAD

        polya_start = self.sva.REP_SVA_POLYA_START

        # 两侧比对都在 polyA 区域 => polyA dominant
        left_in_polya = cand.left_cns_hit_start >= polya_start
        right_in_polya = cand.right_cns_hit_start >= polya_start

        # 如果两侧都在 polyA 区域, 且 clip 证据中没有 head 区域命中
        if left_in_polya and right_in_polya:
            return True

        # 单侧 polyA + 对侧仅在 head 区域很小的 clip => 也算 polyA dominant
        if left_in_polya and cand.right_clip_cns == 0:
            return True
        if right_in_polya and cand.left_clip_cns == 0:
            return True

        return False

    def is_two_side_clip_both_polyA(self, cand: SVACandidate) -> bool:
        """Both sides clipped but both fall in polyA region.

        对应 xTea xtprt_filter.is_two_side_clip_both_polyA_sva。
        """
        polya_start = self.sva.REP_SVA_POLYA_START
        return (
            cand.left_cns_hit_start >= polya_start
            and cand.right_cns_hit_start >= polya_start
        )

    def is_two_side_clip_both_non_polyA(self, cand: SVACandidate) -> bool:
        """Both sides clipped but neither shows any polyA signal.

        对应 xTea xtprt_filter.is_two_side_clip_both_non_polyA_sva:
        两侧都有 clip 但都不在 polyA 区域 => 可能是其他类型的 SV, 不是 SVA 插入。
        """
        polya_start = self.sva.REP_SVA_POLYA_START
        return (
            cand.left_cns_hit_start < polya_start
            and cand.right_cns_hit_start < polya_start
            and not cand.has_polyA
        )

    # ------------------------------------------------------------------
    # Main filter hierarchy
    # 主过滤层次 (对应 xTea post_processing_SVA 的核心逻辑)
    # ------------------------------------------------------------------

    def evaluate_candidate(self, cand: SVACandidate) -> Tuple[bool, str]:
        """Evaluate a single SVA candidate through the filter hierarchy.

        Returns:
            (passed, reason) -- passed=True means the candidate survives filtering.

        过滤层次:
          1. two_side:     排除 polyA dominant, 需通过 AF 质量, 要求 polyA 信号
          2. one_half:     额外检查是否落在低 divergence 的参考 SVA 拷贝内
          3. one_side:     最严格 -- 还排除 orphan transduction 和 transduction
          4. other:        直接过滤
        """
        st = cand.support_type

        # ------- Level 4: other -- 直接丢弃 -------
        if st == OTHER:
            return False, "other_type_skipped"

        # ------- Level 1: two_side -------
        if st == TWO_SIDE:
            if self.is_polyA_dominant_two_side_sva(cand):
                return False, "two_side_polyA_dominant"
            if not cand.af_qualified:
                return False, "two_side_af_unqualified"
            if not cand.has_polyA:
                return False, "two_side_no_polyA"

        # ------- Level 2: one_half_side -------
        elif st == ONE_HALF_SIDE:
            # 落在同类型低 divergence 的参考重复区内 => 可能是 VNTR 扩展假阳性
            if cand.in_repeat and cand.repeat_divergence < self.sva.REP_DIVERGENT_CUTOFF:
                return False, "one_half_in_low_div_repeat"
            if not cand.has_polyA:
                return False, "one_half_no_polyA"
            if self.is_polyA_dominant_two_side_sva(cand):
                return False, "one_half_polyA_dominant"

        # ------- Level 3: one_side -------
        elif st == ONE_SIDE:
            if cand.in_repeat and cand.repeat_divergence < self.sva.REP_DIVERGENT_CUTOFF:
                return False, "one_side_in_low_div_repeat"
            if not cand.has_polyA:
                return False, "one_side_no_polyA"
            # 单侧信号不足以支持 orphan / 常规 transduction (germline 情形)
            if cand.is_orphan_transduction or cand.is_transduction:
                return False, "one_side_transduction_filtered"

        # ------- Cross-level checks (适用于所有通过上面的候选) -------
        # 落在参考 SVA 拷贝内 + 两侧 clip 都在 polyA 区域 => 假阳性
        if cand.in_repeat and self.is_two_side_clip_both_polyA(cand):
            return False, "in_rep_both_clip_polyA"

        # 落在参考拷贝内 + 两侧都没 polyA 信号 => 可能不是新插入
        if cand.in_repeat and self.is_two_side_clip_both_non_polyA(cand):
            return False, "in_rep_both_clip_non_polyA"

        return True, "PASS"

    # ------------------------------------------------------------------
    # Batch filter
    # ------------------------------------------------------------------

    def filter_candidates(
        self, candidates: Sequence[SVACandidate]
    ) -> List[SVACandidate]:
        """Filter a list of SVA candidates; return those that pass.

        Also attaches annotation info (in_SVA_copy / not_in_SVA_copy) to raw_line.
        """
        passed: List[SVACandidate] = []
        for cand in candidates:
            ok, reason = self.evaluate_candidate(cand)
            if ok:
                rep_info = "not_in_SVA_copy"
                if cand.in_repeat:
                    rep_info = f"Fall_in_SVA_copy_{cand.repeat_divergence}"
                cand.raw_line = cand.raw_line.rstrip() + "\t" + rep_info
                passed.append(cand)
            else:
                logger.debug(
                    "Filtered %s:%d reason=%s", cand.chrom, cand.pos, reason
                )
        return passed

    # ------------------------------------------------------------------
    # Integration with MEGAnE output
    # 与 MEGAnE 候选输出的集成接口
    # ------------------------------------------------------------------

    @staticmethod
    def parse_megane_candidate_line(line: str) -> SVACandidate:
        """Parse a single line from MEGAnE's candidate TSV into an SVACandidate.

        MEGAnE BED columns (pre-SVA-filter, 13 cols):
          0: chr
          1: start
          2: end
          3: TE class (e.g. "Retroposon/SVA")
          4: MEI_left:ref_pos=X,chimeric=Y,hybrid=Z,pA=N
          5: MEI_right:ref_pos=X,chimeric=Y,hybrid=Z,pA=N
          6: confidence:high / confidence:low
          7: unique:yes,...
          8: subfamily_pred:status=PASS,MEI=SVA_E,740/764,+/+
          9: 3transduction:no / 3transduction:yes
         10: ID=XXXX

        Extracts real evidence from these fields to populate SVACandidate.
        """
        import re

        fields = line.rstrip().split("\t")
        cand = SVACandidate()
        cand.raw_line = line.rstrip()

        if len(fields) < 6:
            return cand

        cand.chrom = fields[0]
        try:
            cand.pos = int(fields[1])
        except ValueError:
            pass

        # --- Parse chimeric/hybrid counts from col4 and col5 ---
        left_chimeric = 0
        left_hybrid = 0
        right_chimeric = 0
        right_hybrid = 0

        m = re.search(r'chimeric=(\d+)', fields[4])
        if m:
            left_chimeric = int(m.group(1))
        m = re.search(r'hybrid=(\d+)', fields[4])
        if m:
            left_hybrid = int(m.group(1))
        m = re.search(r'chimeric=(\d+)', fields[5])
        if m:
            right_chimeric = int(m.group(1))
        m = re.search(r'hybrid=(\d+)', fields[5])
        if m:
            right_hybrid = int(m.group(1))

        cand.left_clip_cns = left_chimeric
        cand.right_clip_cns = right_chimeric
        cand.left_disc_cns = left_hybrid
        cand.right_disc_cns = right_hybrid

        # --- Infer support_type from evidence ---
        has_left = (left_chimeric > 0)
        has_right = (right_chimeric > 0)
        has_left_disc = (left_hybrid > 0)
        has_right_disc = (right_hybrid > 0)

        if has_left and has_right:
            cand.support_type = TWO_SIDE
        elif (has_left and has_right_disc) or (has_right and has_left_disc):
            cand.support_type = ONE_HALF_SIDE
        elif has_left or has_right:
            cand.support_type = ONE_SIDE
        else:
            cand.support_type = OTHER

        # --- Parse consensus hit positions from col8 (subfamily_pred) ---
        # Patterns:
        #   MEI=SVA_E,740/764,+/+         → single hit range 740-764
        #   MEI_left_breakpoint=SVA_F,18,+ → left hit at pos 18
        #   MEI_right_breakpoint=pA        → right is polyA
        #   pT                             → polyT (reverse polyA)
        if len(fields) > 8:
            pred = fields[8]
            cand.has_polyA = False

            # Check for polyA/polyT indicators
            if 'pA' in pred or 'pT' in pred:
                cand.has_polyA = True

            # Parse single MEI= pattern: MEI=SVA_E,740/764,+/+
            m_single = re.search(r'MEI=\w+,(\d+)/(\d+),', pred)
            if m_single:
                hit_start = int(m_single.group(1))
                hit_end = int(m_single.group(2))
                # Use the range for both sides (single alignment)
                cand.left_cns_hit_start = hit_start
                cand.right_cns_hit_start = hit_end

            # Parse left breakpoint: MEI_left_breakpoint=SVA_F,18,+
            m_left = re.search(r'MEI_left_breakpoint=(?!pA|pT)(\w+),(\d+),', pred)
            if m_left:
                cand.left_cns_hit_start = int(m_left.group(2))

            # Parse right breakpoint: MEI_right_breakpoint=SVA_C,496,-
            m_right = re.search(r'MEI_right_breakpoint=(?!pA|pT)(\w+),(\d+),', pred)
            if m_right:
                cand.right_cns_hit_start = int(m_right.group(2))

            # Left is polyA/polyT → set high cns hit (in polyA region)
            if re.search(r'MEI_left_breakpoint=pA|MEI_left_breakpoint=pT', pred):
                cand.left_cns_hit_start = 2000  # beyond POLYA_START
                cand.has_polyA = True
            if re.search(r'MEI_right_breakpoint=pA|MEI_right_breakpoint=pT', pred):
                cand.right_cns_hit_start = 2000
                cand.has_polyA = True

        # --- Parse transduction flag from col9 ---
        if len(fields) > 9:
            transd = fields[9]
            if '3transduction:yes' in transd:
                cand.is_transduction = True

        return cand

    def filter_megane_output(
        self,
        input_path: str,
        output_path: str,
        *,
        blacklist_regions: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    ) -> int:
        """Read MEGAnE candidate file, apply SVA filtering, write survivors.

        Args:
            input_path:  Path to MEGAnE SVA candidate TSV.
            output_path: Path to write filtered results.
            blacklist_regions: Optional dict {chrom: [(start, end), ...]} of
                regions to exclude (like xTea's x_blklist).

        Returns:
            Number of candidates that passed filtering.

        使用方式:
            filt = SVAFilter()
            n = filt.filter_megane_output("candidates.tsv", "filtered.tsv")
        """
        sva_candidates: List[SVACandidate] = []
        non_sva_lines: List[str] = []

        with open(input_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                te_class = fields[3] if len(fields) > 3 else ""

                # Only apply SVA filtering to SVA/Retroposon entries
                if "SVA" in te_class or "Retroposon" in te_class:
                    cand = self.parse_megane_candidate_line(line)

                    # Blacklist check
                    if blacklist_regions and cand.chrom in blacklist_regions:
                        skip = False
                        for start, end in blacklist_regions[cand.chrom]:
                            if start <= cand.pos <= end:
                                skip = True
                                break
                        if skip:
                            continue

                    sva_candidates.append(cand)
                else:
                    # Non-SVA entries pass through with annotation
                    non_sva_lines.append(line.rstrip() + "\tnot_in_SVA_copy")

        passed = self.filter_candidates(sva_candidates)

        with open(output_path, "w") as fh:
            # Write non-SVA entries first (pass-through)
            for raw_line in non_sva_lines:
                fh.write(raw_line + "\n")
            # Write SVA entries that passed filtering
            for cand in passed:
                fh.write(cand.raw_line + "\n")

        logger.info(
            "SVA filter: %d / %d SVA candidates passed, %d non-SVA passed through",
            len(passed), len(sva_candidates), len(non_sva_lines),
        )
        return len(passed) + len(non_sva_lines)

    # ------------------------------------------------------------------
    # Annotation boundary extension for SVA
    # SVA 注释边界扩展 (对应 xTea SVA_ANNOTATION_EXTND)
    # ------------------------------------------------------------------

    def extend_annotation_boundary(
        self, chrom: str, pos: int, rep_start: int, rep_end: int
    ) -> Tuple[int, int]:
        """Extend repeat annotation boundaries for SVA overlap checking.

        SVA 的 VNTR 区域可能导致注释边界不够精确,
        需要向两侧各扩展 SVA_ANNOTATION_EXTND (默认200bp)。
        """
        ext = self.sva.SVA_ANNOTATION_EXTND
        return max(0, rep_start - ext), rep_end + ext

    # ------------------------------------------------------------------
    # VNTR-aware consistency relaxation
    # VNTR 感知的一致性放松
    # ------------------------------------------------------------------

    def get_relaxed_cluster_diff(self) -> int:
        """Return the relaxed clip-cluster diff cutoff for SVA.

        SVA 的 VNTR 区域会导致两侧 clip cluster 间距大于其他 TE 类型,
        因此使用更宽松的阈值 (200 vs 默认 300)。
        对应 xTea global_values.set_two_clip_cluster_diff_cutoff 的调用逻辑。
        """
        return self.sva.sva_clip_cluster_diff_cutoff


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Quick self-test with synthetic candidates."""
    filt = SVAFilter()

    # Candidate that should PASS (two_side, has polyA, not in repeat)
    c1 = SVACandidate(
        chrom="chr1", pos=100000, support_type=TWO_SIDE,
        has_polyA=True, in_repeat=False,
        left_cns_hit_start=100, right_cns_hit_start=1800,
    )
    ok1, r1 = filt.evaluate_candidate(c1)
    print(f"c1 (two_side, polyA, not_in_rep): passed={ok1}, reason={r1}")

    # Candidate that should FAIL (two_side, polyA dominant)
    c2 = SVACandidate(
        chrom="chr2", pos=200000, support_type=TWO_SIDE,
        has_polyA=True, in_repeat=False,
        left_cns_hit_start=1950, right_cns_hit_start=1960,
    )
    ok2, r2 = filt.evaluate_candidate(c2)
    print(f"c2 (two_side, polyA_dominant)   : passed={ok2}, reason={r2}")

    # Candidate that should FAIL (one_side, in low-div repeat)
    c3 = SVACandidate(
        chrom="chr3", pos=300000, support_type=ONE_SIDE,
        has_polyA=True, in_repeat=True, repeat_divergence=5.0,
    )
    ok3, r3 = filt.evaluate_candidate(c3)
    print(f"c3 (one_side, low_div_repeat)   : passed={ok3}, reason={r3}")

    # Relaxed cluster diff
    print(f"\nRelaxed cluster diff for SVA: {filt.get_relaxed_cluster_diff()}")
    print(f"Annotation extension example: {filt.extend_annotation_boundary('chr1', 5000, 4000, 6000)}")


if __name__ == "__main__":
    _cli_test()
