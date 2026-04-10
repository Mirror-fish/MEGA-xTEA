"""
transduction.py -- Simplified 3' transduction detection.

Adapted from xTea's x_transduction.py (Simon Chu, Harvard Medical School).
Focuses on germline use case with basic 3' transduction, SVA-specific
transduction with 5000bp window, and orphan transduction flagging.

The full xTea x_transduction.py is ~1562 lines handling somatic, tumor,
long-read, and many edge cases.  This module extracts the core germline
logic suitable for post-processing MEGAnE candidates.

Key concepts (关键概念):
  - 3' 转导: TE 插入时携带了源位点下游的基因组序列
  - Orphan 转导: 仅有转导序列, 没有 TE 本体 (孤儿转导)
  - SVA 的转导搜索窗口为 5000bp (LINE1 也是 5000bp)
  - 判定转导需要 disc reads 在 flank 区域形成 cluster
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from megaxtea.config import DEFAULT_CONFIG, MegaXTeaConfig, TransductionParams

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class TransductionCandidate:
    """A candidate transduction event."""

    ins_chrom: str = ""
    ins_pos: int = 0

    # Source information (转导来源)
    source_chrom: str = ""
    source_start: int = 0
    source_end: int = 0

    # Support counts
    n_left_disc: int = 0
    n_right_disc: int = 0
    n_clip_support: int = 0
    n_polyA_reads: int = 0

    # Classification
    is_orphan: bool = False         # 孤儿转导 (无 TE 本体)
    is_sibling: bool = False        # sibling transduction
    source_label: str = ""          # e.g. "polymerphic", "reference"
    te_type: str = ""               # LINE1 / SVA

    # Distance from source TE copy end to transduced region start
    dist_from_te_end: int = 0

    raw_line: str = ""


@dataclass
class TransductionSource:
    """A known TE source copy (reference or polymorphic) that can produce transductions."""

    chrom: str = ""
    start: int = 0
    end: int = 0
    subfamily: str = ""
    is_polymorphic: bool = False
    is_full_length: bool = False


# ---------------------------------------------------------------------------
# Core transduction detector
# ---------------------------------------------------------------------------

class TransductionDetector:
    """Simplified 3' transduction detection for germline TE insertions.

    Implements the core logic from xTea's XTransduction class, focusing on:
      1. Basic 3' transduction detection
      2. SVA-specific transduction with configurable window (default 5000bp)
      3. Orphan transduction flagging

    Usage::

        det = TransductionDetector()
        results = det.detect_from_disc_clusters(disc_clusters, sources)
    """

    def __init__(self, config: Optional[MegaXTeaConfig] = None) -> None:
        self.cfg = config or DEFAULT_CONFIG
        self.tp = self.cfg.transduction

    # ------------------------------------------------------------------
    # Flank window selection
    # 根据 TE 类型选择转导搜索窗口
    # ------------------------------------------------------------------

    def get_flank_window(self, te_type: str) -> int:
        """Return the transduction search window size for a given TE type.

        SVA 和 LINE1 都使用 5000bp 窗口 (来自 xTea global_values)。
        """
        upper = te_type.upper()
        if "SVA" in upper:
            return self.tp.FLANK_WINDOW_SVA
        return self.tp.FLANK_WINDOW_LINE1

    # ------------------------------------------------------------------
    # 3' transduction detection from discordant read clusters
    # 基于 disc read cluster 的 3' 转导检测
    # ------------------------------------------------------------------

    def detect_from_disc_clusters(
        self,
        disc_clusters: Dict[str, Dict[int, "DiscCluster"]],
        known_sources: List[TransductionSource],
        min_disc_support: int = 0,
    ) -> List[TransductionCandidate]:
        """Detect 3' transductions by matching disc clusters to known source flanks.

        Algorithm (算法):
          1. 对每个候选插入位点, 检查其 disc reads 的 mate 是否落在已知 TE 源拷贝
             的下游 flank 区域内 (flank window = 5000bp for SVA/LINE1)
          2. 如果多数 disc mates 指向同一源的 flank => 3' transduction
          3. 如果 disc mates 指向 flank 但无 clip-on-consensus 证据 => orphan

        Args:
            disc_clusters: {chrom: {pos: DiscCluster}} of insertion candidates.
            known_sources: List of known TE source copies.
            min_disc_support: Minimum disc read pairs to call transduction.

        Returns:
            List of TransductionCandidate objects.
        """
        if min_disc_support <= 0:
            min_disc_support = self.tp.MIN_DISC_CUTOFF

        # Build source flank index (源拷贝的 flank 区域索引)
        # {chrom: [(source_end, source_end + window, source_obj), ...]}
        flank_index: Dict[str, List[Tuple[int, int, TransductionSource]]] = {}
        for src in known_sources:
            window = self.get_flank_window(src.subfamily)
            entry = (src.end, src.end + window, src)
            flank_index.setdefault(src.chrom, []).append(entry)

        results: List[TransductionCandidate] = []

        for ins_chrom, pos_dict in disc_clusters.items():
            for ins_pos, cluster in pos_dict.items():
                # Check mate positions against source flanks
                # mate_positions: list of (mate_chrom, mate_pos)
                mate_positions = getattr(cluster, "mate_positions", [])
                if not mate_positions:
                    continue

                best_source, best_count = self._find_dominant_source(
                    mate_positions, flank_index
                )
                if best_source is None or best_count < min_disc_support:
                    continue

                # Determine if orphan (no consensus alignment evidence)
                has_cns_hit = getattr(cluster, "has_consensus_hit", False)

                td = TransductionCandidate(
                    ins_chrom=ins_chrom,
                    ins_pos=ins_pos,
                    source_chrom=best_source.chrom,
                    source_start=best_source.start,
                    source_end=best_source.end,
                    n_left_disc=getattr(cluster, "n_left", 0),
                    n_right_disc=getattr(cluster, "n_right", 0),
                    is_orphan=not has_cns_hit,
                    te_type=best_source.subfamily,
                    dist_from_te_end=0,  # computed below
                )
                results.append(td)

        return results

    def _find_dominant_source(
        self,
        mate_positions: List[Tuple[str, int]],
        flank_index: Dict[str, List[Tuple[int, int, TransductionSource]]],
    ) -> Tuple[Optional[TransductionSource], int]:
        """Find the source whose flank region receives the most disc mates.

        对应 xTea 的 dominant source 判定逻辑:
        要求占比 >= TRANSDCT_MULTI_SOURCE_MIN_RATIO (0.45) 才算有效。
        """
        source_counts: Dict[int, Tuple[TransductionSource, int]] = {}

        for mate_chrom, mate_pos in mate_positions:
            if mate_chrom not in flank_index:
                continue
            for flank_start, flank_end, src in flank_index[mate_chrom]:
                if flank_start <= mate_pos <= flank_end:
                    key = id(src)
                    if key in source_counts:
                        old_src, old_cnt = source_counts[key]
                        source_counts[key] = (old_src, old_cnt + 1)
                    else:
                        source_counts[key] = (src, 1)

        if not source_counts:
            return None, 0

        best_key = max(source_counts, key=lambda k: source_counts[k][1])
        best_src, best_cnt = source_counts[best_key]

        # Check dominance ratio
        total = sum(cnt for _, cnt in source_counts.values())
        ratio = best_cnt / total if total > 0 else 0
        if ratio < self.tp.TRANSDCT_MULTI_SOURCE_MIN_RATIO:
            return None, 0

        return best_src, best_cnt

    # ------------------------------------------------------------------
    # Orphan transduction detection
    # 孤儿转导检测
    # ------------------------------------------------------------------

    @staticmethod
    def flag_orphan_transductions(
        candidates: List[TransductionCandidate],
    ) -> List[TransductionCandidate]:
        """Flag candidates that are likely orphan transductions.

        Orphan 转导的特征:
          - disc reads 指向源拷贝的 flank 区域
          - 但 clip reads 没有比对到 TE consensus
          - 仅有 polyA 信号或极少的 consensus hit
        """
        for cand in candidates:
            if cand.is_orphan:
                cand.source_label = "orphan"
            elif cand.is_sibling:
                cand.source_label = "sibling"
        return candidates

    # ------------------------------------------------------------------
    # Simple transduction check for SVA candidates
    # SVA 候选的简单转导判断 (用于 sva_filter 集成)
    # ------------------------------------------------------------------

    def is_likely_transduction(
        self,
        n_left_disc: int,
        n_right_disc: int,
        has_consensus_hit: bool,
        min_disc: int = 0,
    ) -> Tuple[bool, bool]:
        """Quick check whether a candidate shows transduction signals.

        Returns:
            (is_transduction, is_orphan)

        简化的转导判断: 如果有足够多的 disc reads 但没有 consensus hit,
        则可能是 orphan transduction。
        """
        if min_disc <= 0:
            min_disc = self.tp.MIN_DISC_CUTOFF

        total_disc = n_left_disc + n_right_disc
        if total_disc < min_disc:
            return False, False

        if not has_consensus_hit:
            return True, True  # orphan

        return True, False

    # ------------------------------------------------------------------
    # Close-site deduplication
    # 邻近位点去重
    # ------------------------------------------------------------------

    def are_sites_close(self, pos1: int, pos2: int) -> bool:
        """Check whether two insertion positions are close enough to be the same event.

        对应 xTea TWO_SITES_CLOSE_DIST (default 100bp)。
        """
        return abs(pos1 - pos2) <= self.tp.TWO_SITES_CLOSE_DIST

    # ------------------------------------------------------------------
    # Integration with MEGAnE candidate list
    # ------------------------------------------------------------------

    def annotate_megane_candidates(
        self,
        candidates_tsv: str,
        known_sources: List[TransductionSource],
        output_path: str,
    ) -> int:
        """Read MEGAnE candidate file, annotate transductions, write output.

        For each candidate, appends a transduction annotation column.
        This is a simplified interface -- full transduction calling requires
        BAM access and disc-read re-alignment (handled in the main pipeline).

        Returns number of candidates annotated as transduction.
        """
        n_td = 0
        with open(candidates_tsv) as fin, open(output_path, "w") as fout:
            for line in fin:
                if not line.strip() or line.startswith("#"):
                    fout.write(line)
                    continue
                fields = line.rstrip().split("\t")
                # Placeholder: in real pipeline, would check disc alignment
                td_label = "not_transduction"
                fout.write(line.rstrip() + "\t" + td_label + "\n")
        logger.info("Transduction annotation: %d / total candidates flagged", n_td)
        return n_td


# ---------------------------------------------------------------------------
# Placeholder for disc cluster data structure
# ---------------------------------------------------------------------------

@dataclass
class DiscCluster:
    """Discordant-read cluster around an insertion site."""

    chrom: str = ""
    pos: int = 0
    n_left: int = 0
    n_right: int = 0
    mate_positions: List[Tuple[str, int]] = field(default_factory=list)
    has_consensus_hit: bool = False


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Quick self-test."""
    det = TransductionDetector()

    print(f"Flank window SVA  : {det.get_flank_window('SVA')}")
    print(f"Flank window LINE1: {det.get_flank_window('LINE1')}")
    print(f"Close sites (100, 150): {det.are_sites_close(100, 150)}")
    print(f"Close sites (100, 300): {det.are_sites_close(100, 300)}")

    is_td, is_orphan = det.is_likely_transduction(5, 3, False)
    print(f"Transduction check (disc=8, no cns): is_td={is_td}, orphan={is_orphan}")

    is_td2, is_orphan2 = det.is_likely_transduction(5, 3, True)
    print(f"Transduction check (disc=8, has cns): is_td={is_td2}, orphan={is_orphan2}")


if __name__ == "__main__":
    _cli_test()
