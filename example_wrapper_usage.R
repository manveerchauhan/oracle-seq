# ORACLE-seq Wrapper Function Usage Example
# 
# This example demonstrates how to use the run_oracle_seq() wrapper function
# for organized, user-friendly analysis.

# Load the wrapper function
source("run_oracle_seq.R")

# Example 1: Basic usage
# ======================
# This will create a directory structure like:
# oracle_results/
# ├── lung_cancer_study/
# │   ├── analysis_config.yaml
# │   ├── logs/
# │   ├── statistical_analysis/
# │   ├── plots/
# │   └── reports/

run_oracle_seq(
  oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
  data_file = "data/rarefaction_curves.rds",
  analysis_name = "lung_cancer_study",
  output_dir = "oracle_results",
  thresholds = c(10, 1, 0.1),
  bootstrap_samples = 1000,
  use_bayesian_model_averaging = TRUE,
  confidence_level = 0.95
)

# Example 2: Custom thresholds
# ============================
# Use different marginal returns thresholds

run_oracle_seq(
  oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
  data_file = "data/rarefaction_curves.rds",
  analysis_name = "high_resolution_analysis",
  output_dir = "oracle_results",
  thresholds = c(20, 10, 5, 1, 0.5, 0.1),
  bootstrap_samples = 1000,
  use_bayesian_model_averaging = TRUE,
  confidence_level = 0.95
)

# Example 3: Quick analysis with fewer bootstrap samples
# =====================================================
# For faster analysis during development

run_oracle_seq(
  oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
  data_file = "data/rarefaction_curves.rds",
  analysis_name = "quick_test",
  output_dir = "oracle_results",
  thresholds = c(10, 1, 0.1),
  bootstrap_samples = 100,  # Fewer samples for speed
  use_bayesian_model_averaging = FALSE,          # Disable BMA for speed
  confidence_level = 0.95
)

# Example 4: Running multiple analyses
# ===================================
# The wrapper automatically handles duplicate names

# First run
run_oracle_seq(
  oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
  data_file = "data/protocol_comparison.rds",
  analysis_name = "protocol_comparison",
  output_dir = "oracle_results",
  thresholds = c(10, 1, 0.1),
  bootstrap_samples = 1000,
  use_bayesian_model_averaging = TRUE,
  confidence_level = 0.95
)

# Second run with same name - will create "protocol_comparison_2"
run_oracle_seq(
  oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
  data_file = "data/protocol_comparison_updated.rds",
  analysis_name = "protocol_comparison",  # Same name
  output_dir = "oracle_results",
  thresholds = c(15, 3, 0.5),           # Different parameters
  bootstrap_samples = 1000,
  use_bayesian_model_averaging = TRUE,
  confidence_level = 0.95
)

# What happens if parameters are missing?
# =======================================
# The function will stop with helpful error messages:

# This will fail with clear error message
# run_oracle_seq(
#   oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
#   data_file = "data/rarefaction_curves.rds",
#   analysis_name = "incomplete_analysis",
#   output_dir = "oracle_results"
#   # Missing: thresholds, bootstrap_samples, use_bayesian_model_averaging, confidence_level
# )

# Error output would be:
# ERROR: Missing required parameter(s):
#   - thresholds
#   - bootstrap_samples
#   - use_bayesian_model_averaging
#   - confidence_level
# 
# Please provide all required parameters and try again.

# Directory structure created:
# ===========================
# oracle_results/
# ├── lung_cancer_study/
# │   ├── analysis_config.yaml          # Parameters used
# │   ├── logs/
# │   │   └── analysis_log.txt         # Full analysis log
# │   ├── statistical_analysis/
# │   │   ├── fitted_models.rds        # Statistical results
# │   │   ├── model_competition_report.txt
# │   │   └── threshold_predictions.csv
# │   ├── plots/
# │   │   ├── rarefaction_curves.png   # All plots
# │   │   ├── marginal_returns.png
# │   │   └── model_comparison.png
# │   └── reports/
# │       ├── analysis_summary.html
# │       ├── model_equations.txt
# │       └── threshold_recommendations.txt
# ├── high_resolution_analysis/
# │   └── ... (same structure with different results)
# ├── quick_test/
# │   └── ... (same structure)
# ├── protocol_comparison/
# │   └── ... (first run results)
# └── protocol_comparison_2/
#     └── ... (second run results) 