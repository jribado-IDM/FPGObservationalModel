"""FPG Observational Model - Genomic surveillance sampling for EMOD simulations.

This package provides tools for converting EMOD's Full Parasite Genetics (FPG)
simulation output into recapitulative sampling for genomic surveillance analysis.
"""

# Import main functions from modules
from .run_observational_model import (
    run_observational_model,
    get_default_config,
    update_matrix_indices,
    extract_sampled_infections,
)

from .unified_sampling import (
    run_sampling_model,
    subset_randomly,
    subset_by_seasons,
    subset_by_age,
)

from .unified_metric_calculations import (
    register_matrix,
    get_matrix,
    run_time_summaries,
    comprehensive_group_summary,
)

# Define public API
__all__ = [
    # Main workflow
    "run_observational_model",
    "get_default_config",
    "update_matrix_indices",
    "extract_sampled_infections",

    # Sampling functions
    "run_sampling_model",
    "subset_randomly",
    "subset_by_seasons",
    "subset_by_age",

    # Metric calculations
    "register_matrix",
    "get_matrix",
    "run_time_summaries",
    "comprehensive_group_summary",
]