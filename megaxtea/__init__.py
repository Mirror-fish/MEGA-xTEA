"""
MEGA-xTEA: Merged MEGAnE + xTea integration package.

Combines MEGAnE's fast k-mer-based TE detection with xTea's SVA-specific
filtering, ML genotyping, and transduction detection for improved accuracy.

Modules:
    sva_filter    - SVA post-processing filter (from xTea x_post_filter logic)
    ml_genotype   - ML-based genotype classification (from xTea sklearn classifier)
    polyA_detector - Advanced polyA signal detection (from xTea x_polyA)
    transduction  - Simplified 3' transduction detection (from xTea x_transduction)
    te_classifier - Enhanced TE classification combining MEGAnE + xTea
    config        - Central configuration parameters
"""

__version__ = "0.1.0"
__author__ = "MEGA-xTEA contributors"
__all__ = [
    "sva_filter",
    "ml_genotype",
    "polyA_detector",
    "transduction",
    "te_classifier",
    "config",
]
