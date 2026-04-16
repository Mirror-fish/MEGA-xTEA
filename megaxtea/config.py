"""
MEGA-xTEA central configuration.

Consolidates parameters from:
  - MEGAnE: coverage-adaptive thresholds, read-length estimation, BLAST e-values
  - xTea:   SVA-specific cutoffs, polyA constants, ML genotyping settings

All tunables live here so that every other module imports from one place.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict


# ---------------------------------------------------------------------------
# 来自 xTea global_values.py 的 SVA 相关常量
# (SVA structure-aware parameters ported from xTea)
# ---------------------------------------------------------------------------

@dataclass
class SVAParams:
    """SVA-specific detection parameters (源自 xTea x_post_filter / global_values)."""

    # Consensus region lengths
    REP_SVA_CNS_HEAD: int = 400        # SVA consensus 5' head length
    REP_SVA_POLYA_START: int = 1900    # polyA region start position in SVA consensus
    REP_SVA_MIN_LEN: int = 300         # minimum SVA insertion length

    # Cluster / filter cutoffs
    sva_clip_cluster_diff_cutoff: int = 200  # relaxed cluster diff for SVA VNTR
    SVA_ANNOTATION_EXTND: int = 200          # annotation boundary extension for SVA

    # Divergence
    REP_DIVERGENT_CUTOFF: int = 15     # max div-rate to keep as ref-copy overlap
    REP_LOW_DIVERGENT_CUTOFF: int = 7  # stricter cutoff for LINE1

    # Two-clip cluster default (non-SVA) -- restored after SVA filtering
    TWO_CLIP_CLUSTER_DIFF_CUTOFF: int = 300


# ---------------------------------------------------------------------------
# 来自 xTea 的 polyA 检测参数
# ---------------------------------------------------------------------------

@dataclass
class PolyAParams:
    """PolyA detection parameters (源自 xTea x_polyA / global_values)."""

    # Canonical polyA signals (正向/反向互补)
    CANONICAL_POLYA: tuple = ("AATAAA", "ATTAAA")
    CANONICAL_POLYT: tuple = ("TTTATT", "TTTAAT")

    # TE-type-specific polyA region start positions
    # 不同 TE 类型的 polyA 起始位置
    POLYA_START_SVA: int = 1900
    POLYA_START_ALU: int = 255
    POLYA_START_LINE1: int = 5950

    # Detection windows
    CK_POLYA_CLIP_WIN: int = 25
    CK_POLYA_SEQ_MAX: int = 20
    POLYA_RATIO: float = 0.4
    DOMINANT_POLYA_RATIO: float = 0.75
    ONE_SIDE_POLYA_CUTOFF: float = 0.75
    MAX_POLYA_RATIO: float = 0.85

    # Minimum consecutive A/T bases
    N_MIN_A_T: int = 5


# ---------------------------------------------------------------------------
# 来自 xTea 的 ML 基因分型参数
# ---------------------------------------------------------------------------

@dataclass
class MLGenotypeParams:
    """Machine-learning genotyping parameters (源自 xTea sklearn classifier)."""

    n_features: int = 15
    n_estimators: int = 20
    test_size: float = 0.3
    random_state: int = 0

    # BAM re-scan constants (源自 xTea global_values.py)
    # 用于精确特征提取的窗口和阈值参数
    TSD_CUTOFF: int = 100              # raw clip counting window (bp)
    CLIP_EXACT_CLIP_SLACK: int = 3     # effective (AF) clip counting window (bp)
    DFT_IS: int = 550                  # default insert size for disc/concordant
    DISC_THRESHOLD: int = 2000         # distance threshold for discordant pairs
    BWA_HALF_READ_MIN_SCORE: int = 45  # half read length for left/right full-map

    # Feature column indices in xTea output (0-based)
    # 特征向量列索引 (对应 xTea 的输出格式)
    IDX_LCLIP_CNS: int = 5
    IDX_RCLIP_CNS: int = 6
    IDX_LDISC_CNS: int = 7
    IDX_RDISC_CNS: int = 8
    IDX_LPOLYA: int = 9
    IDX_RPOLYA: int = 10
    IDX_LCOV: int = 11
    IDX_RCOV: int = 12
    IDX_EFF_CLIP: int = 35
    IDX_EFF_FMAP: int = 36
    IDX_RAW_LCLIP: int = 37
    IDX_RAW_RCLIP: int = 38
    IDX_DISC: int = 39
    IDX_CONC: int = 40


# ---------------------------------------------------------------------------
# 来自 MEGAnE 的覆盖度自适应参数
# ---------------------------------------------------------------------------

@dataclass
class MEGAnEParams:
    """Coverage-adaptive parameters from MEGAnE pipeline."""

    # BLAST classification
    overhang_evalue_threshold: float = 1e-5

    # Coverage estimation
    n_depth_samples: int = 200         # reads sampled for depth estimation
    max_cov_times: int = 3             # abnormal coverage multiplier

    # Genotyping (Gaussian-fit fallback)
    gaussian_min_evidence: int = 3     # min allele-evidence reads for Gaussian fit
    gaussian_n_components: int = 3     # GMM components (0/0, 0/1, 1/1)


# ---------------------------------------------------------------------------
# Transduction parameters (源自 xTea x_transduction / global_values)
# ---------------------------------------------------------------------------

@dataclass
class TransductionParams:
    """Transduction detection parameters (简化自 xTea)."""

    FLANK_WINDOW_SVA: int = 5000       # SVA-specific transduction search window
    FLANK_WINDOW_LINE1: int = 5000
    MIN_DISC_CUTOFF: int = 2
    MIN_POLYMORPHIC_SOURCE_DIST: int = 1000
    TRANSDCT_MULTI_SOURCE_MIN_RATIO: float = 0.45
    F_MIN_TRSDCT_DISC_MAP_RATIO: float = 0.65
    TWO_SITES_CLOSE_DIST: int = 100
    TD_REP_DIVERGENT_CUTOFF: int = 5


# ---------------------------------------------------------------------------
# 来自 xTea 的 TE 分类编码
# ---------------------------------------------------------------------------

@dataclass
class TEClassParams:
    """TE classification binary encoding (源自 xTea x_rep_type)."""

    # Binary flags (bitmask)
    LINE1: int = 1
    ALU: int = 2
    SVA: int = 4
    HERV: int = 8
    MIT: int = 16
    MSTA: int = 32
    PSEUDOGENE: int = 64

    # String labels
    LABEL_LINE1: str = "LINE1"
    LABEL_ALU: str = "ALU"
    LABEL_SVA: str = "SVA"
    LABEL_HERV: str = "HERV"

    # VCF ALT representations
    VCF_LINE1: str = "INS:ME:LINE1"
    VCF_ALU: str = "INS:ME:ALU"
    VCF_SVA: str = "INS:ME:SVA"
    VCF_HERV: str = "INS:ME:HERV-K"


# ---------------------------------------------------------------------------
# Clip / disc read filter constants (来自 xTea global_values)
# ---------------------------------------------------------------------------

@dataclass
class ClipDiscParams:
    """Clip and discordant-read filter constants."""

    INITIAL_MIN_CLIP_CUTOFF: int = 2
    MINIMUM_CLIP_MAPQ: int = 12
    MINIMUM_DISC_MAPQ: int = 20
    MIN_CLIP_MAPPED_RATIO: float = 0.65
    MIN_DISC_MAPPED_RATIO: float = 0.7
    TWO_SIDE_CLIP_MIN_LEN: int = 8
    NCLIP_HALF_CUTOFF: int = 5
    INDEL_READS_MAX_RATIO: float = 0.3
    ABNORMAL_COV_TIMES: int = 2


# ---------------------------------------------------------------------------
# Top-level config aggregator
# ---------------------------------------------------------------------------

@dataclass
class MegaXTeaConfig:
    """Top-level configuration aggregating all sub-configs.

    Usage::

        cfg = MegaXTeaConfig()
        print(cfg.sva.REP_SVA_CNS_HEAD)   # 400
        print(cfg.ml.n_features)           # 15
    """

    sva: SVAParams = field(default_factory=SVAParams)
    polya: PolyAParams = field(default_factory=PolyAParams)
    ml: MLGenotypeParams = field(default_factory=MLGenotypeParams)
    megane: MEGAnEParams = field(default_factory=MEGAnEParams)
    transduction: TransductionParams = field(default_factory=TransductionParams)
    te_class: TEClassParams = field(default_factory=TEClassParams)
    clip_disc: ClipDiscParams = field(default_factory=ClipDiscParams)


# Singleton default config
DEFAULT_CONFIG = MegaXTeaConfig()


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Print default configuration for verification."""
    import json
    from dataclasses import asdict

    cfg = MegaXTeaConfig()
    print(json.dumps(asdict(cfg), indent=2, default=str))


if __name__ == "__main__":
    _cli_test()
