"""
ml_genotype.py -- ML-based genotype classification.

Ported from xTea's x_genotype_classify_sklearn.py (Simon Chu, Harvard).
Uses a RandomForestClassifier to predict insertion genotypes (0/1, 1/1)
from 15-dimensional feature vectors extracted from candidate sites.

Falls back to MEGAnE's Gaussian-mixture fitting when no trained model is
available.

Feature vector (15维特征向量, 与 xTea 一致):
  0  lclipcns   -- left  clip aligned to consensus / left  coverage
  1  rclipcns   -- right clip aligned to consensus / right coverage
  2  ldisccns   -- left  disc aligned to consensus / left  coverage
  3  rdisccns   -- right disc aligned to consensus / right coverage
  4  polyA      -- (left_polyA/lcov + right_polyA/rcov)
  5  lcov       -- left  local coverage (raw)
  6  rcov       -- right local coverage (raw)
  7  clip       -- effective clip reads / total coverage
  8  fullmap    -- effective fully-mapped / total coverage
  9  clipratio  -- clip / (clip + fullmap)
  10 discratio  -- discordant / (discordant + concordant)
  11 rawlclip   -- raw left  clip / lcov
  12 rawrclip   -- raw right clip / rcov
  13 discordant -- discordant pairs / total coverage
  14 concordant -- concordant pairs / total coverage

Origin:
  - Feature extraction & RF classifier: xTea (x_genotype_classify_sklearn.py)
  - Gaussian fallback concept: MEGAnE (merge_allele_evidence_ins.py)
"""

from __future__ import annotations

import logging
import os
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)

# Feature names in order (与 xTea arff 头一致)
FEATURE_NAMES: Tuple[str, ...] = (
    "lclipcns", "rclipcns", "ldisccns", "rdisccns",
    "polyA",
    "lcov", "rcov",
    "clip", "fullmap", "clipratio", "discratio",
    "rawlclip", "rawrclip",
    "discordant", "concordant",
)

N_FEATURES = len(FEATURE_NAMES)  # 15

# Genotype label mapping (xTea convention)
# xTea 将 0/0 编码为 0, 0/1 编码为 1, 1/1 编码为 2
# RF 训练时只用 {1, 2} (het / hom), 因为 0/0 候选已在上游过滤掉
GENOTYPE_MAP = {0: "0/0", 1: "0/1", 2: "1/1"}


# ---------------------------------------------------------------------------
# Feature extraction
# 特征提取 (对应 xTea _parser_features)
# ---------------------------------------------------------------------------

def extract_features_from_xtea_fields(fields: List[str]) -> List[float]:
    """Extract 15-dimensional feature vector from xTea-format output fields.

    *fields* is the whitespace-split line from xTea's candidate output.
    Column indices follow xTea's convention (see config.MLGenotypeParams).

    关键计算逻辑:
      - 所有 clip/disc 计数都除以本地覆盖度进行归一化
      - clipratio = clip / (clip + fullmap), 反映断裂读段占比
      - discratio = disc / (disc + concordant), 反映不和谐配对占比
    """
    f_lcov = float(fields[11]) or 1e-10
    f_rcov = float(fields[12]) or 1e-10
    total_cov = f_lcov + f_rcov

    features: List[float] = []

    # 0-3: clip/disc aligned to consensus, normalized by local coverage
    features.append(float(fields[5]) / f_lcov)    # lclipcns
    features.append(float(fields[6]) / f_rcov)    # rclipcns
    features.append(float(fields[7]) / f_lcov)    # ldisccns
    features.append(float(fields[8]) / f_rcov)    # rdisccns

    # 4: polyA (combined left + right, normalized)
    features.append(float(fields[9]) / f_lcov + float(fields[10]) / f_rcov)

    # 5-6: raw local coverage
    features.append(f_lcov)
    features.append(f_rcov)

    # 7-8: effective clip and fullmap, normalized
    eff_clip = float(fields[35])
    eff_fmap = float(fields[36])
    features.append(eff_clip / total_cov)
    features.append(eff_fmap / total_cov)

    # 9: clip ratio
    denom_clip = eff_clip + eff_fmap
    features.append(eff_clip / denom_clip if denom_clip > 0 else 0.0)

    # 10: disc ratio
    disc = float(fields[39])
    conc = float(fields[40])
    denom_disc = disc + conc
    features.append(disc / denom_disc if denom_disc > 0 else 0.0)

    # 11-14: raw clip & disc/concordant, normalized
    features.append(float(fields[37]) / f_lcov)    # rawlclip
    features.append(float(fields[38]) / f_rcov)    # rawrclip
    features.append(disc / total_cov)               # discordant
    features.append(conc / total_cov)               # concordant

    return features


def extract_features_from_megane_evidence(evidence: Dict[str, Any]) -> List[float]:
    """Extract features from MEGAnE's allele-evidence dictionary.

    MEGAnE stores evidence differently from xTea; this function maps MEGAnE
    keys to the 15-dim vector.  Missing keys default to 0.

    MEGAnE 的 allele evidence 字典键名可能不同, 此函数做映射。
    """
    def _g(key: str) -> float:
        return float(evidence.get(key, 0))

    lcov = _g("left_coverage") or 1e-10
    rcov = _g("right_coverage") or 1e-10
    total_cov = lcov + rcov

    eff_clip = _g("clipped_reads")
    eff_fmap = _g("fully_mapped_reads")
    disc = _g("discordant_pairs")
    conc = _g("concordant_pairs")

    denom_clip = eff_clip + eff_fmap
    denom_disc = disc + conc

    return [
        _g("left_clip_consensus") / lcov,
        _g("right_clip_consensus") / rcov,
        _g("left_disc_consensus") / lcov,
        _g("right_disc_consensus") / rcov,
        _g("left_polyA") / lcov + _g("right_polyA") / rcov,
        lcov,
        rcov,
        eff_clip / total_cov,
        eff_fmap / total_cov,
        eff_clip / denom_clip if denom_clip > 0 else 0.0,
        disc / denom_disc if denom_disc > 0 else 0.0,
        _g("raw_left_clip") / lcov,
        _g("raw_right_clip") / rcov,
        disc / total_cov,
        conc / total_cov,
    ]


def extract_features_from_megane_genotyped_bed(fields: List[str]) -> List[float]:
    """Extract 15-dim feature vector from MEGAnE's genotyped BED columns.

    MEGAnE genotyped BED (15 columns):
      col4:  MEI_left:ref_pos=X,chimeric=Y,hybrid=Z
      col5:  MEI_right:ref_pos=X,chimeric=Y,hybrid=Z
      col12: tsd_depth: "allele;ratio;structure"
      col13: spanning:  "allele;count"
      col14: disc:      "allele;count"

    Maps available MEGAnE evidence to xTea's 15 features with best-effort
    approximation where exact values are not available.

    注意: MEGAnE 不提供 LCOV/RCOV, polyA counts, concordant pairs 等原始值。
    这里使用近似值。对于精确的 ML 基因分型, 未来需要修改 MEGAnE 的
    输出管线以保存原始计数。
    """
    import re

    # --- Parse chimeric/hybrid from col4/col5 ---
    left_chimeric = 0
    left_hybrid = 0
    right_chimeric = 0
    right_hybrid = 0

    if len(fields) > 4:
        m = re.search(r'chimeric=(\d+)', fields[4])
        if m:
            left_chimeric = int(m.group(1))
        m = re.search(r'hybrid=(\d+)', fields[4])
        if m:
            left_hybrid = int(m.group(1))
    if len(fields) > 5:
        m = re.search(r'chimeric=(\d+)', fields[5])
        if m:
            right_chimeric = int(m.group(1))
        m = re.search(r'hybrid=(\d+)', fields[5])
        if m:
            right_hybrid = int(m.group(1))

    # --- Parse spanning and disc counts from col12/col13 ---
    spanning_count = 0
    disc_count = 0

    if len(fields) > 12:
        parts = fields[12].split(";")
        if len(parts) >= 2:
            try:
                spanning_count = int(float(parts[1]))
            except (ValueError, IndexError):
                pass
    if len(fields) > 13:
        parts = fields[13].split(";")
        if len(parts) >= 2:
            try:
                disc_count = int(float(parts[1]))
            except (ValueError, IndexError):
                pass

    # --- Estimate coverage ---
    # Without direct LCOV/RCOV values, estimate from total evidence.
    # A 30x genome typically has ~30 reads at each position.
    # Use disc + spanning + chimeric as a rough total.
    total_support = left_chimeric + right_chimeric + left_hybrid + right_hybrid
    eff_clip = float(left_chimeric + right_chimeric)
    eff_disc = float(disc_count)
    eff_fmap = float(spanning_count)

    # Rough coverage estimate: total support reads + spanning + disc
    # This is imprecise but provides a reasonable normalization factor.
    est_total = max(eff_clip + eff_fmap + eff_disc + 1, 10)
    est_lcov = est_total / 2.0
    est_rcov = est_total / 2.0
    total_cov = est_lcov + est_rcov

    # Concordant estimate: total_reads - disc - clip
    est_conc = max(est_total - eff_disc - eff_clip, 0)

    denom_clip = eff_clip + eff_fmap
    denom_disc = eff_disc + est_conc

    features = [
        left_chimeric / est_lcov,                             # lclipcns
        right_chimeric / est_rcov,                            # rclipcns
        left_hybrid / est_lcov,                               # ldisccns
        right_hybrid / est_rcov,                              # rdisccns
        0.0,                                                  # polyA (not available)
        est_lcov,                                             # lcov
        est_rcov,                                             # rcov
        eff_clip / total_cov,                                 # clip
        eff_fmap / total_cov,                                 # fullmap
        eff_clip / denom_clip if denom_clip > 0 else 0.0,     # clipratio
        eff_disc / denom_disc if denom_disc > 0 else 0.0,     # discratio
        left_chimeric / est_lcov,                             # rawlclip
        right_chimeric / est_rcov,                            # rawrclip
        eff_disc / total_cov,                                 # discordant
        est_conc / total_cov,                                 # concordant
    ]
    return features


# ---------------------------------------------------------------------------
# ML Genotyper
# ---------------------------------------------------------------------------

class MLGenotyper:
    """ML-based genotype classifier supporting both sklearn RandomForest
    and Deep Forest (CascadeForestClassifier) models.

    Wraps either model type with xTea-compatible feature extraction,
    training, saving/loading, and prediction interfaces.

    Usage::

        gtyper = MLGenotyper()
        gtyper.load_model("DF21_model_1_2")  # directory = Deep Forest
        genotypes = gtyper.predict(feature_matrix)
    """

    def __init__(self, n_estimators: int = 20, random_state: int = 0) -> None:
        self.n_estimators = n_estimators
        self.random_state = random_state
        self.model: Any = None
        self._model_type: str = "none"  # "sklearn", "deepforest", or "none"

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def train(
        self,
        X: np.ndarray,
        y: np.ndarray,
        test_size: float = 0.3,
        verbose: bool = True,
    ) -> float:
        """Train a RandomForestClassifier on feature matrix *X* and labels *y*.

        Labels: 1 = heterozygous (0/1), 2 = homozygous (1/1).

        Returns test-set accuracy.

        训练逻辑与 xTea GntpClassifier_sklearn.train_model 一致:
          使用 sklearn train_test_split 划分, RF n_estimators=20。
        """
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import accuracy_score
        from sklearn.model_selection import train_test_split

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=self.random_state
        )
        clf = RandomForestClassifier(
            n_jobs=-1,
            random_state=self.random_state,
            n_estimators=self.n_estimators,
        )
        clf.fit(X_train, y_train)
        self.model = clf

        preds = clf.predict(X_test)
        acc = accuracy_score(y_test, preds)
        if verbose:
            logger.info("RF genotyper trained -- test accuracy: %.4f", acc)
        return float(acc)

    def train_from_xtea_files(
        self,
        sf_het_list: str,
        sf_hom_list: str,
        model_out: str,
        balance: bool = False,
        test_size: float = 0.3,
    ) -> float:
        """Train from xTea-format output file lists (one file path per line).

        *sf_het_list*: file listing paths to 0/1 sample outputs.
        *sf_hom_list*: file listing paths to 1/1 sample outputs.

        对应 xTea gnrt_training_arff_from_xTEA_output, 但跳过 arff 中间格式,
        直接读取特征并训练。
        """
        all_X: List[List[float]] = []
        all_y: List[int] = []

        # Load hom (label=2)
        n_hom = 0
        with open(sf_hom_list) as fh:
            for path_line in fh:
                feats = self._load_features_from_xtea_output(path_line.strip(), True)
                for f in feats:
                    all_X.append(f)
                    all_y.append(2)
                n_hom += 1

        # Load het (label=1)
        n_het = 0
        with open(sf_het_list) as fh:
            for path_line in fh:
                if balance and n_het >= n_hom:
                    break
                feats = self._load_features_from_xtea_output(path_line.strip(), True)
                for f in feats:
                    all_X.append(f)
                    all_y.append(1)
                n_het += 1

        X = np.array(all_X, dtype=float)
        y = np.array(all_y, dtype=int)
        acc = self.train(X, y, test_size=test_size)
        self.save_model(model_out)
        return acc

    # ------------------------------------------------------------------
    # Prediction
    # ------------------------------------------------------------------

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict genotypes for feature matrix *X*.

        Returns array of labels (1=het, 2=hom).
        """
        if self.model is None:
            raise RuntimeError("No model loaded. Call load_model() or train() first.")
        return self.model.predict(X)

    def predict_genotype_strings(self, X: np.ndarray) -> List[str]:
        """Predict and return VCF-style genotype strings ("0/1" or "1/1")."""
        preds = self.predict(X)
        return [GENOTYPE_MAP.get(int(p), "0/0") for p in preds]

    def predict_for_xtea_file(
        self, sf_xtea: str, sf_output: str
    ) -> None:
        """Read xTea output, predict genotypes, write augmented output.

        对应 xTea predict_for_site -- 逐行读取, 追加预测的基因型。
        """
        features = self._load_features_from_xtea_output(sf_xtea, b_train=False)
        if not features:
            logger.warning("No features loaded from %s", sf_xtea)
            return

        X = np.array(features, dtype=float)
        preds = self.predict(X)

        with open(sf_xtea) as fin, open(sf_output, "w") as fout:
            for i, line in enumerate(fin):
                gt = GENOTYPE_MAP.get(int(preds[i]), "0/0")
                fout.write(line.rstrip() + "\t" + gt + "\n")

    # ------------------------------------------------------------------
    # Model persistence
    # ------------------------------------------------------------------

    def save_model(self, path: str) -> None:
        """Save trained model to pickle file."""
        if self.model is None:
            raise RuntimeError("No model to save.")
        with open(path, "wb") as fh:
            pickle.dump(self.model, fh)
        logger.info("Model saved to %s", path)

    def load_model(self, path: str) -> None:
        """Load model from file or directory.

        Supports two formats:
          - Directory (e.g. DF21_model_1_2/) → Deep Forest CascadeForestClassifier
          - .pkl file → sklearn RandomForestClassifier (pickle)

        对应 xTea 的两种分类器:
          GntpClassifier_DF21 使用 CascadeForestClassifier.load(目录)
          GntpClassifier_sklearn 使用 pickle.load(文件)
        """
        if os.path.isdir(path):
            # Deep Forest model (directory containing param.pkl + estimator/)
            try:
                from deepforest import CascadeForestClassifier
                self.model = CascadeForestClassifier()
                self.model.load(path)
                self._model_type = "deepforest"
                logger.info("Deep Forest model loaded from %s", path)
            except ImportError:
                logger.error(
                    "deepforest package not installed. "
                    "Install with: pip install deep-forest"
                )
                raise
        else:
            # sklearn pickle model
            try:
                import joblib
                self.model = joblib.load(path)
                self._model_type = "sklearn"
                logger.info("sklearn model loaded via joblib from %s", path)
            except Exception:
                # Fallback to raw pickle with Python 2 compatibility
                with open(path, "rb") as fh:
                    if sys.version_info >= (3, 0):
                        self.model = pickle.load(fh, encoding="latin1")
                    else:
                        self.model = pickle.load(fh)
                self._model_type = "sklearn"
                logger.info("sklearn model loaded via pickle from %s", path)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _load_features_from_xtea_output(
        sf_xtea: str, b_train: bool = True
    ) -> List[List[float]]:
        """Load feature vectors from xTea output file.

        对应 xTea load_in_feature_from_xTEA_output。
        """
        all_features: List[List[float]] = []
        if not os.path.isfile(sf_xtea):
            return all_features
        with open(sf_xtea) as fh:
            for line in fh:
                fields = line.split()
                if len(fields) < 41:
                    continue
                feats = extract_features_from_xtea_fields(fields)
                all_features.append(feats)
        return all_features


# ---------------------------------------------------------------------------
# Gaussian-mixture fallback (来自 MEGAnE 的高斯拟合基因分型)
# ---------------------------------------------------------------------------

class GaussianGenotypeFallback:
    """Gaussian-mixture-model fallback genotyper.

    When no trained RF model is available, use a simple GMM on the
    allele-evidence ratio to call 0/0, 0/1, 1/1.

    This mirrors MEGAnE's approach in merge_allele_evidence_ins.py where
    a Gaussian fit is applied to the distribution of supporting-read ratios.

    基本思路: 对 clip ratio (支持插入的读段占总读段的比例) 拟合三个高斯分量,
    分别对应 0/0, 0/1, 1/1。
    """

    def __init__(self, n_components: int = 3, min_evidence: int = 3) -> None:
        self.n_components = n_components
        self.min_evidence = min_evidence

    def genotype_by_ratio(
        self,
        clip_reads: int,
        total_reads: int,
        het_low: float = 0.15,
        het_high: float = 0.75,
    ) -> str:
        """Simple ratio-based genotyping (no GMM, just thresholds).

        当样本量不足以拟合 GMM 时使用阈值法:
          ratio < het_low  => 0/0
          het_low <= ratio <= het_high => 0/1
          ratio > het_high => 1/1
        """
        if total_reads < self.min_evidence:
            return "./."
        ratio = clip_reads / total_reads
        if ratio < het_low:
            return "0/0"
        elif ratio <= het_high:
            return "0/1"
        else:
            return "1/1"

    def genotype_batch_gmm(
        self, ratios: np.ndarray
    ) -> List[str]:
        """Fit a GMM to *ratios* and assign genotypes.

        ratios: 1-D array of supporting-read ratios for all candidates.

        使用 sklearn GMM 对比值分布拟合, 根据组分均值排序分配基因型。
        """
        try:
            from sklearn.mixture import GaussianMixture
        except ImportError:
            logger.warning("sklearn not available; using threshold fallback.")
            return [self.genotype_by_ratio(int(r * 100), 100) for r in ratios]

        if len(ratios) < self.n_components * 2:
            return [self.genotype_by_ratio(int(r * 100), 100) for r in ratios]

        gmm = GaussianMixture(
            n_components=self.n_components, random_state=0
        )
        X = ratios.reshape(-1, 1)
        gmm.fit(X)

        # Sort components by mean to assign labels
        means = gmm.means_.flatten()
        order = np.argsort(means)
        label_map = {order[0]: "0/0", order[1]: "0/1", order[2]: "1/1"}

        predictions = gmm.predict(X)
        return [label_map[p] for p in predictions]


# ---------------------------------------------------------------------------
# Unified genotyping interface
# 统一基因分型接口
# ---------------------------------------------------------------------------

class UnifiedGenotyper:
    """Unified genotyper: tries ML model first, falls back to Gaussian.

    Usage::

        gt = UnifiedGenotyper(model_path="rf_model.pkl")
        genotypes = gt.genotype_candidates(feature_matrix)
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        n_estimators: int = 20,
    ) -> None:
        self.ml = MLGenotyper(n_estimators=n_estimators)
        self.fallback = GaussianGenotypeFallback()
        self._has_model = False

        if model_path and (os.path.isfile(model_path) or os.path.isdir(model_path)):
            self.ml.load_model(model_path)
            self._has_model = True

    def genotype_candidates(
        self,
        features: np.ndarray,
    ) -> List[str]:
        """Genotype candidates using ML model or fallback.

        Args:
            features: (N, 15) feature matrix.

        Returns:
            List of genotype strings ("0/1", "1/1", etc.).
        """
        if self._has_model:
            try:
                return self.ml.predict_genotype_strings(features)
            except Exception as exc:
                logger.warning("ML prediction failed (%s); using fallback.", exc)

        # Fallback: use clip-ratio column (index 9) for Gaussian fitting
        if features.shape[1] > 9:
            ratios = features[:, 9]
        else:
            ratios = np.zeros(features.shape[0])
        return self.fallback.genotype_batch_gmm(ratios)


# ---------------------------------------------------------------------------
# CLI test
# ---------------------------------------------------------------------------

def _cli_test() -> None:
    """Quick self-test with synthetic data."""
    print("Feature names:", FEATURE_NAMES)
    print("N_FEATURES:", N_FEATURES)

    # Synthetic features
    rng = np.random.RandomState(42)
    X_het = rng.rand(20, N_FEATURES) * 0.5
    X_hom = rng.rand(20, N_FEATURES) * 0.5 + 0.5
    X = np.vstack([X_het, X_hom])
    y = np.array([1] * 20 + [2] * 20)

    gt = MLGenotyper()
    try:
        acc = gt.train(X, y, test_size=0.3)
        print(f"Training accuracy: {acc:.3f}")
        preds = gt.predict_genotype_strings(X[:5])
        print(f"Sample predictions: {preds}")
    except ImportError:
        print("sklearn not installed; skipping ML test.")

    # Fallback test
    fb = GaussianGenotypeFallback()
    print(f"Ratio-based genotype (10/100): {fb.genotype_by_ratio(10, 100)}")
    print(f"Ratio-based genotype (50/100): {fb.genotype_by_ratio(50, 100)}")
    print(f"Ratio-based genotype (90/100): {fb.genotype_by_ratio(90, 100)}")


if __name__ == "__main__":
    _cli_test()
