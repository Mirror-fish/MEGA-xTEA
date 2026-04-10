"""
polyA_detector.py -- Advanced polyA signal detection.

Ported from xTea's x_polyA.py (Simon Chu, Harvard Medical School).
Provides position-aware polyA detection for different TE types, canonical
signal search, and consecutive stretch detection with orientation.

Key concepts (关键算法说明):
  - 不同 TE 类型的 polyA 区域起始位置不同:
      SVA:   1900 bp
      Alu:   255 bp
      LINE1: 5950 bp
  - 典型 polyA 信号: AATAAA / TTTATT (正向/反向互补)
  - 变体信号: ATTAAA / TTTAAT, AGTAAA, TATAAA 等
  - 连续 polyA/T 段检测可带方向 (left-clip => polyA, right-clip => polyT)
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Set, Tuple

from megaxtea.config import DEFAULT_CONFIG, PolyAParams


# ---------------------------------------------------------------------------
# Canonical polyA signal patterns
# 典型 polyA 信号模式
# ---------------------------------------------------------------------------

# Primary signals
CANONICAL_FORWARD = ("AATAAA", "ATTAAA")   # polyA (sense)
CANONICAL_REVERSE = ("TTTATT", "TTTAAT")   # polyT (antisense / RC)

# Extended variant signals documented in xTea
VARIANT_SIGNALS_FWD = (
    "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
    "AATATA", "AATACA", "AATAGA", "AAAAAG",
    "ACTAAA", "AAGAAA", "AATGAA", "TTTAAA",
    "AAAACA", "GGGGCT",
)


class PolyADetector:
    """Advanced polyA signal detection engine.

    Combines xTea's polyA detection heuristics with configurable parameters
    suitable for post-processing MEGAnE candidates.
    """

    def __init__(self, params: Optional[PolyAParams] = None) -> None:
        self.p = params or DEFAULT_CONFIG.polya

    # ------------------------------------------------------------------
    # Basic polyA / polyT composition checks (来自 xTea PolyA class)
    # ------------------------------------------------------------------

    @staticmethod
    def is_poly_A_or_T(seq: str, ratio: float = 0.75) -> bool:
        """Return True if *seq* is dominated (>= *ratio*) by A or T.

        对应 xTea PolyA.is_poly_A_T -- 判断序列是否为高纯度 polyA/T。
        """
        if not seq:
            return False
        upper = seq.upper()
        n = len(upper)
        cutoff = n * ratio
        return upper.count("A") >= cutoff or upper.count("T") >= cutoff

    @staticmethod
    def is_dominant_polyA(seq: str, ratio: float = 0.75) -> bool:
        """Check if A or T dominates *seq* above *ratio*.

        对应 xTea PolyA.is_dominant_polyA。
        """
        if not seq:
            return False
        upper = seq.upper()
        n = len(upper)
        return (upper.count("A") / n > ratio) or (upper.count("T") / n > ratio)

    @staticmethod
    def is_dominant_A_only(seq: str, ratio: float = 0.75) -> bool:
        """Check if only A (not T) dominates *seq*."""
        if not seq:
            return False
        upper = seq.upper()
        return upper.count("A") / len(upper) > ratio

    # ------------------------------------------------------------------
    # Consecutive polyA/T stretch detection
    # 连续 polyA/T 段检测
    # ------------------------------------------------------------------

    @staticmethod
    def has_consecutive_polyA_T(seq: str) -> bool:
        """Return True if *seq* contains AAAAA, TTTTT, AATAA, or TTATT.

        对应 xTea PolyA.is_consecutive_polyA_T。
        """
        upper = seq.upper()
        return any(p in upper for p in ("AAAAA", "TTTTT", "AATAA", "TTATT"))

    @staticmethod
    def has_consecutive_polyA_T_strict(seq: str) -> bool:
        """Stricter version requiring 6-mers (对应 xTea is_consecutive_polyA_T2)."""
        upper = seq.upper()
        return any(p in upper for p in ("AAAAAA", "TTTTTT", "AAATAA", "TTTATT"))

    @staticmethod
    def has_consecutive_polyA_T_with_orientation(seq: str) -> Tuple[bool, Optional[bool]]:
        """Detect consecutive polyA/T with orientation information.

        Returns:
            (found, is_polyA) -- is_polyA=True means A-rich (sense),
            False means T-rich (antisense).  None if not found.

        对应 xTea PolyA.is_consecutive_polyA_T_with_ori:
          左侧 clip => 预期 polyA (正向)
          右侧 clip => 预期 polyT (反向互补)
        """
        upper = seq.upper()
        if "AAAAA" in upper or "AATAA" in upper:
            return True, True
        if "TTTTT" in upper or "TTATT" in upper:
            return True, False
        return False, None

    @staticmethod
    def has_oriented_polyA_T(seq: str, is_left_clip: bool) -> bool:
        """Orientation-aware check: left-clip expects polyA, right-clip expects polyT.

        对应 xTea PolyA.is_consecutive_polyA_T_with_oritation -- 方向感知的检测:
        左侧 clipped reads 应只含 polyA, 右侧 clipped reads 应只含 polyT。
        """
        upper = seq.upper()
        if is_left_clip:
            return "AAAAA" in upper or "AATAA" in upper
        else:
            return "TTTTT" in upper or "TTATT" in upper

    # ------------------------------------------------------------------
    # Consecutive stretch length measurement
    # 最长连续 A/T 段长度
    # ------------------------------------------------------------------

    @staticmethod
    def max_consecutive_A_or_T(seq: str) -> Tuple[int, int]:
        """Return (max_consecutive_A, max_consecutive_T) in *seq*."""
        max_a = max_t = cur_a = cur_t = 0
        for ch in seq.upper():
            if ch == "A":
                cur_a += 1
            else:
                max_a = max(max_a, cur_a)
                cur_a = 0
            if ch == "T":
                cur_t += 1
            else:
                max_t = max(max_t, cur_t)
                cur_t = 0
        return max(max_a, cur_a), max(max_t, cur_t)

    def contains_poly_A_T(self, seq: str, min_consecutive: int = 5) -> bool:
        """True if *seq* has >= *min_consecutive* consecutive A or T bases.

        对应 xTea PolyA.contain_poly_A_T。
        """
        ma, mt = self.max_consecutive_A_or_T(seq)
        return ma >= min_consecutive or mt >= min_consecutive

    @staticmethod
    def contains_enough_A_or_T(seq: str, min_count: int = 5) -> bool:
        """True if total A or T count >= *min_count* (not necessarily consecutive)."""
        upper = seq.upper()
        return upper.count("A") >= min_count or upper.count("T") >= min_count

    # ------------------------------------------------------------------
    # Canonical polyA signal search
    # 典型 polyA 信号搜索 (AATAAA / ATTAAA 及其反向互补)
    # ------------------------------------------------------------------

    @staticmethod
    def contains_canonical_signal(seq: str, reverse_complement: bool = False) -> bool:
        """Check for canonical polyA hexamer signals.

        对应 xTea PolyA.contain_polyA_T(s_seq, b_rc)。
        """
        upper = seq.upper()
        if reverse_complement:
            return any(s in upper for s in ("TTTATT", "TTTTTT", "TTTAAT"))
        return any(s in upper for s in ("AATAAA", "AAAAAA", "ATTAAA"))

    @staticmethod
    def search_polyA_signal_positions(seq: str, reverse_complement: bool = False) -> List[int]:
        """Find all positions of canonical polyA signals in *seq*.

        Returns 0-based start positions.
        对应 xTea PolyA.search_multi_polyA_locations。
        """
        upper = seq.upper()
        positions: List[int] = []
        if reverse_complement:
            patterns = ("TTTATT", "TTTTTT", "TTTAAT")
        else:
            patterns = ("AATAAA", "AAAAAA", "ATTAAA")
        for pat in patterns:
            positions.extend(m.start() for m in re.finditer(pat, upper))
        return sorted(positions)

    # ------------------------------------------------------------------
    # TE-type-aware polyA detection
    # 根据 TE 类型确定 polyA 检测区域
    # ------------------------------------------------------------------

    def get_polya_start_for_te(self, te_type: str) -> int:
        """Return the expected polyA-region start position for a given TE type.

        不同 TE 类型的 consensus 序列中 polyA 区域起始位置不同:
          SVA   -> 1900 bp
          ALU   -> 255 bp
          LINE1 -> 5950 bp
        """
        te_upper = te_type.upper()
        if "SVA" in te_upper:
            return self.p.POLYA_START_SVA
        if "ALU" in te_upper:
            return self.p.POLYA_START_ALU
        if "LINE" in te_upper or "L1" in te_upper:
            return self.p.POLYA_START_LINE1
        # default: use SVA (most conservative)
        return self.p.POLYA_START_SVA

    def is_in_polya_region(self, hit_pos: int, te_type: str, tolerance: int = 100) -> bool:
        """Check whether a consensus hit position falls in the expected polyA region.

        用于判断 clipped reads 在 consensus 上的比对位置是否落在 polyA 区域,
        若是则该 clip 证据可能为 polyA 假阳性。
        """
        start = self.get_polya_start_for_te(te_type)
        return hit_pos >= (start - tolerance)

    # ------------------------------------------------------------------
    # Pre-defined rmsk polyA/T patterns (来自 xTea)
    # ------------------------------------------------------------------

    @staticmethod
    def get_predefined_polyA_rmsk() -> Set[str]:
        """Return pre-defined polyA repeat labels from RepeatMasker.

        对应 xTea PolyA.get_pre_defined_polyA_in_rmsk -- 包含所有旋转变体。
        """
        seeds = [
            "(AAAAAC)n", "(AAAAC)n", "(AAAATAA)n", "(AAAAT)n",
            "(AAATA)n", "(AAAT)n", "(AACAAA)n", "(AAC)n",
            "(AATAAAA)n", "(AATAAA)n", "(AATAA)n", "(AATA)n",
            "(AAT)n", "(ACA)n", "(A)n",
        ]
        result: Set[str] = set()
        for seed in seeds:
            result.add(seed)
            core = seed[1:-2]  # strip parens and 'n'
            for i in range(len(core)):
                rotated = core[i:] + core[:i]
                result.add(f"({rotated})n")
        return result

    @staticmethod
    def get_predefined_polyT_rmsk() -> Set[str]:
        """Return pre-defined polyT repeat labels from RepeatMasker."""
        seeds = [
            "(GTTTTT)n", "(GTTTT)n", "(TTATTTT)n", "(ATTTT)n",
            "(TATTT)n", "(ATTT)n", "(TTTGTT)n", "(GTT)n",
            "(TTTTGTT)n", "(TTTATT)n", "(TTATT)n", "(TATT)n",
            "(ATT)n", "(TGT)n", "(T)n",
        ]
        result: Set[str] = set()
        for seed in seeds:
            result.add(seed)
            core = seed[1:-2]
            for i in range(len(core)):
                rotated = core[i:] + core[:i]
                result.add(f"({rotated})n")
        return result


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Quick self-test."""
    det = PolyADetector()

    test_seq = "ACGTAATAAAGGGGGAAAAAA"
    print(f"Sequence: {test_seq}")
    print(f"  is_poly_A_or_T        : {det.is_poly_A_or_T(test_seq)}")
    print(f"  has_consecutive       : {det.has_consecutive_polyA_T(test_seq)}")
    print(f"  canonical signal (fwd): {det.contains_canonical_signal(test_seq, False)}")
    print(f"  signal positions (fwd): {det.search_polyA_signal_positions(test_seq, False)}")
    print(f"  max consec A/T        : {det.max_consecutive_A_or_T(test_seq)}")
    print(f"  oriented (left clip)  : {det.has_oriented_polyA_T(test_seq, True)}")
    print()
    print(f"PolyA start SVA  : {det.get_polya_start_for_te('SVA')}")
    print(f"PolyA start ALU  : {det.get_polya_start_for_te('ALU')}")
    print(f"PolyA start LINE1: {det.get_polya_start_for_te('LINE1')}")
    print()
    print(f"Predefined rmsk polyA patterns (sample): {list(det.get_predefined_polyA_rmsk())[:5]}")


if __name__ == "__main__":
    _cli_test()
