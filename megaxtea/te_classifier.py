"""
te_classifier.py -- Enhanced TE classification.

Combines:
  - MEGAnE's BLAST-based classification (parse_blastn_result.py)
  - xTea's binary encoding for output types (x_rep_type.py)
  - SVA subfamily awareness

Origin:
  - BLAST parsing logic: MEGAnE (Shohei Kojima, RIKEN)
  - Binary type encoding & VCF labels: xTea (Simon Chu, Harvard)

TE 分类逻辑说明:
  MEGAnE 使用 BLAST 将候选序列比对到 consensus, 按 e-value 取最佳 hit。
  xTea 使用二进制位掩码编码 TE 类型 (LINE1=1, ALU=2, SVA=4, HERV=8 ...)。
  本模块将两者统一, 支持从 BLAST 结果或特征向量判定 TE 类型, 并输出
  VCF 兼容的类型标签。
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from megaxtea.config import DEFAULT_CONFIG, MegaXTeaConfig, TEClassParams

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Binary encoding (来自 xTea x_rep_type.py)
# ---------------------------------------------------------------------------

# Bitmask constants
TE_LINE1 = 1
TE_ALU = 2
TE_SVA = 4
TE_HERV = 8
TE_MIT = 16
TE_MSTA = 32
TE_PSEUDOGENE = 64


# ---------------------------------------------------------------------------
# Classification result
# ---------------------------------------------------------------------------

@dataclass
class TEClassification:
    """Result of TE classification for a single candidate."""

    # Primary type
    te_type: str = "unknown"         # e.g. "LINE1", "ALU", "SVA", "HERV"
    te_bitmask: int = 0              # xTea-style binary encoding

    # BLAST hit details (from MEGAnE)
    blast_subject: str = ""          # best-hit subject name
    blast_evalue: float = 1.0
    blast_start: int = 0             # 0-based start on subject
    blast_end: int = 0
    blast_strand: str = "+"

    # SVA subfamily (if applicable)
    sva_subfamily: str = ""          # e.g. "SVA_A", "SVA_D", "SVA_E", "SVA_F"

    # VCF representation
    vcf_alt: str = ""                # e.g. "INS:ME:SVA"


# ---------------------------------------------------------------------------
# TEClassifier
# ---------------------------------------------------------------------------

class TEClassifier:
    """Enhanced TE classification engine.

    Usage::

        clf = TEClassifier()
        result = clf.classify_from_blast_hit("SVA_D", 100, 1800, 1e-50)
        print(result.te_type, result.vcf_alt)
    """

    def __init__(self, config: Optional[MegaXTeaConfig] = None) -> None:
        self.cfg = config or DEFAULT_CONFIG
        self.tc = self.cfg.te_class

        # Subject-name to type mapping (covers common RepeatMasker / Dfam names)
        # 科/超科名称到 TE 类型的映射
        self._name_map: Dict[str, str] = {}
        self._build_name_map()

    def _build_name_map(self) -> None:
        """Build mapping from consensus/subfamily names to TE type labels."""
        # LINE1 variants
        for prefix in ("L1", "LINE1", "L1HS", "L1PA", "L1PB", "L1MA", "L1MB", "L1MC", "L1MD"):
            self._name_map[prefix] = self.tc.LABEL_LINE1

        # ALU variants
        for prefix in ("ALU", "ALUJ", "ALUS", "ALUY", "ALUSC", "ALUSQ", "ALUSP", "ALUSN"):
            self._name_map[prefix] = self.tc.LABEL_ALU

        # SVA variants
        for prefix in ("SVA", "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"):
            self._name_map[prefix] = self.tc.LABEL_SVA

        # HERV-K
        for prefix in ("HERV", "HERVK", "HERV-K", "HERVK11", "LTR5_HS"):
            self._name_map[prefix] = self.tc.LABEL_HERV

    # ------------------------------------------------------------------
    # Name-based classification
    # ------------------------------------------------------------------

    def _resolve_type_from_name(self, subject_name: str) -> str:
        """Map a BLAST subject name / subfamily name to a TE type label.

        首先尝试精确匹配, 再尝试大写前缀匹配。
        """
        upper = subject_name.upper().replace("-", "")

        # Exact match
        if upper in self._name_map:
            return self._name_map[upper]

        # Prefix match (e.g. "L1PA2" -> "L1PA" -> "LINE1")
        for prefix, label in sorted(self._name_map.items(), key=lambda x: -len(x[0])):
            if upper.startswith(prefix):
                return label

        # Keyword fallback
        if "LINE" in upper or "L1" in upper:
            return self.tc.LABEL_LINE1
        if "ALU" in upper:
            return self.tc.LABEL_ALU
        if "SVA" in upper:
            return self.tc.LABEL_SVA
        if "HERV" in upper or "LTR" in upper:
            return self.tc.LABEL_HERV

        return "unknown"

    # ------------------------------------------------------------------
    # Bitmask encoding (来自 xTea RepType)
    # ------------------------------------------------------------------

    def type_to_bitmask(self, te_type: str) -> int:
        """Convert a type label to xTea binary encoding.

        xTea 使用位掩码编码, 可以组合多种类型。
        """
        upper = te_type.upper()
        if "LINE" in upper or "L1" in upper:
            return TE_LINE1
        if "ALU" in upper:
            return TE_ALU
        if "SVA" in upper:
            return TE_SVA
        if "HERV" in upper:
            return TE_HERV
        return 0

    def bitmask_to_type(self, bitmask: int) -> str:
        """Convert xTea bitmask back to primary type label.

        对应 xTea RepType.get_rep_type -- 返回第一个匹配的类型。
        """
        if bitmask & TE_LINE1:
            return self.tc.LABEL_LINE1
        if bitmask & TE_ALU:
            return self.tc.LABEL_ALU
        if bitmask & TE_SVA:
            return self.tc.LABEL_SVA
        if bitmask & TE_HERV:
            return self.tc.LABEL_HERV
        if bitmask & TE_MIT:
            return "Mitochondria"
        if bitmask & TE_MSTA:
            return "Msta"
        if bitmask & TE_PSEUDOGENE:
            return "Pseudogene"
        return "unknown"

    def bitmask_to_all_types(self, bitmask: int) -> List[str]:
        """Return all types encoded in a bitmask.

        对应 xTea RepType.get_all_rep_types。
        """
        types: List[str] = []
        for flag, label in [
            (TE_LINE1, self.tc.LABEL_LINE1),
            (TE_ALU, self.tc.LABEL_ALU),
            (TE_SVA, self.tc.LABEL_SVA),
            (TE_HERV, self.tc.LABEL_HERV),
        ]:
            if bitmask & flag:
                types.append(label)
        return types

    # ------------------------------------------------------------------
    # VCF ALT field
    # ------------------------------------------------------------------

    def type_to_vcf_alt(self, te_type: str) -> str:
        """Convert type label to VCF-style ALT field.

        对应 xTea RepType 的 VCF 标签 (INS:ME:LINE1 等)。
        """
        upper = te_type.upper()
        if "LINE" in upper or "L1" in upper:
            return self.tc.VCF_LINE1
        if "ALU" in upper:
            return self.tc.VCF_ALU
        if "SVA" in upper:
            return self.tc.VCF_SVA
        if "HERV" in upper:
            return self.tc.VCF_HERV
        return f"INS:ME:{te_type}"

    # ------------------------------------------------------------------
    # SVA subfamily detection
    # SVA 亚家族识别
    # ------------------------------------------------------------------

    @staticmethod
    def detect_sva_subfamily(subject_name: str) -> str:
        """Extract SVA subfamily from subject/consensus name.

        SVA has subfamilies A-F with different evolutionary ages:
          SVA_A (oldest) -> SVA_F (youngest, most active).
        """
        upper = subject_name.upper().replace("-", "_")
        for sub in ("SVA_F", "SVA_E", "SVA_D", "SVA_C", "SVA_B", "SVA_A"):
            if sub in upper:
                return sub
        if "SVA" in upper:
            return "SVA"
        return ""

    # ------------------------------------------------------------------
    # Main classification from BLAST hit
    # ------------------------------------------------------------------

    def classify_from_blast_hit(
        self,
        subject_name: str,
        hit_start: int,
        hit_end: int,
        evalue: float,
        evalue_threshold: float = 0.0,
    ) -> TEClassification:
        """Classify a candidate from its best BLAST hit.

        Args:
            subject_name: BLAST subject (consensus) name.
            hit_start: 0-based start on subject.
            hit_end: End on subject.
            evalue: BLAST e-value.
            evalue_threshold: Maximum e-value to accept.  0 = use config default.

        MEGAnE 的 parse_blastn_result.py 提取最佳 BLAST hit, 本函数在此基础上
        确定 TE 类型、亚家族和 VCF 标签。
        """
        if evalue_threshold <= 0:
            evalue_threshold = self.cfg.megane.overhang_evalue_threshold

        result = TEClassification(
            blast_subject=subject_name,
            blast_evalue=evalue,
            blast_start=hit_start,
            blast_end=hit_end,
        )

        if evalue > evalue_threshold:
            result.te_type = "unknown"
            return result

        result.te_type = self._resolve_type_from_name(subject_name)
        result.te_bitmask = self.type_to_bitmask(result.te_type)
        result.vcf_alt = self.type_to_vcf_alt(result.te_type)

        if result.te_type == self.tc.LABEL_SVA:
            result.sva_subfamily = self.detect_sva_subfamily(subject_name)

        return result

    # ------------------------------------------------------------------
    # Batch classification from MEGAnE BLAST output
    # ------------------------------------------------------------------

    def classify_blast_output(
        self,
        blast_parsed_path: str,
    ) -> Dict[str, TEClassification]:
        """Classify all candidates from MEGAnE's parsed BLAST result file.

        File format (from parse_blastn_result.py):
            query_name<TAB>subject,start,end,strand[;...]<TAB>evalue

        Returns:
            {query_name: TEClassification}
        """
        results: Dict[str, TEClassification] = {}
        try:
            with open(blast_parsed_path) as fh:
                for line in fh:
                    parts = line.rstrip().split("\t")
                    if len(parts) < 3 or parts[0] == "NA":
                        continue
                    qname = parts[0]
                    hit_info = parts[1].split(";")[0]  # take first hit
                    evalue = float(parts[2])

                    hit_parts = hit_info.split(",")
                    if len(hit_parts) < 4:
                        continue
                    subject = hit_parts[0]
                    start = int(hit_parts[1])
                    end = int(hit_parts[2])

                    results[qname] = self.classify_from_blast_hit(
                        subject, start, end, evalue
                    )
        except FileNotFoundError:
            logger.warning("BLAST parsed file not found: %s", blast_parsed_path)

        return results


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Quick self-test."""
    clf = TEClassifier()

    # Test classification
    for name, start, end, ev in [
        ("SVA_D", 100, 1800, 1e-50),
        ("L1HS", 5000, 6000, 1e-100),
        ("AluYa5", 10, 280, 1e-80),
        ("HERVK11", 200, 5000, 1e-30),
        ("unknown_seq", 0, 50, 0.5),
    ]:
        r = clf.classify_from_blast_hit(name, start, end, ev)
        print(f"  {name:15s} => type={r.te_type:8s}  vcf={r.vcf_alt:18s}  "
              f"subfamily={r.sva_subfamily or '-':6s}  mask={r.te_bitmask}")

    # Bitmask round-trip
    print("\nBitmask round-trips:")
    for t in ("LINE1", "ALU", "SVA", "HERV"):
        bm = clf.type_to_bitmask(t)
        print(f"  {t} -> {bm} -> {clf.bitmask_to_type(bm)}")

    # Combined bitmask
    combined = TE_LINE1 | TE_SVA
    print(f"\nCombined mask {combined}: {clf.bitmask_to_all_types(combined)}")


if __name__ == "__main__":
    _cli_test()
