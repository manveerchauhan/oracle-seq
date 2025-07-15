# ORACLE-seq Wrapper Function
# 
# Simple wrapper interface for running ORACLE-seq analysis with organized output
# directories and parameter validation.

#' Run ORACLE-seq Analysis with Organized Output
#'
#' Main wrapper function that provides a simple interface for running complete
#' ORACLE-seq rarefaction curve analysis with automatic output organization.
#'
#' @param oracle_seq_dir Character string specifying path to oracle-seq directory containing R/ and main.R
#' @param data_file Character string specifying path to rarefaction data file (RDS format)
#' @param analysis_name Character string specifying name for the analysis (used for output directory)
#' @param output_dir Character string specifying base output directory path
#' @param thresholds Numeric vector of marginal returns thresholds (features/million reads)
#' @param bootstrap_samples Integer number of bootstrap samples for uncertainty quantification
#' @param use_bayesian_model_averaging Logical indicating whether to use Bayesian Model Averaging
#' @param confidence_level Numeric confidence level for intervals (0-1)
#'
#' @return NULL (side effect: creates analysis outputs in organized directory structure)
#'
#' @details
#' **Directory Structure Created:**
#' ```
#' output_dir/
#' ├── analysis_name/           # or analysis_name_2, analysis_name_3, etc.
#' │   ├── analysis_config.yaml # Parameters used for this analysis
#' │   ├── logs/               # Analysis logs and warnings
#' │   ├── statistical_analysis/ # Model results and reports
#' │   ├── plots/              # All visualization outputs
#' │   └── reports/            # Text reports and summaries
#' ```
#'
#' **Parameter Validation:**
#' All required parameters must be provided. The function will stop with
#' helpful error messages if any required parameters are missing.
#'
#' **Directory Naming:**
#' If the specified analysis_name already exists in output_dir, the function
#' automatically appends "_2", "_3", etc. until an available name is found.
#'
#' @examples
#' # Run complete ORACLE-seq analysis
#' run_oracle_seq(
#'   oracle_seq_dir = "/path/to/oracle-seq-dev",
#'   data_file = "data/rarefaction_curves.rds",
#'   analysis_name = "lung_cancer_study",
#'   output_dir = "oracle_results",
#'   thresholds = c(10, 1, 0.1),
#'   bootstrap_samples = 1000,
#'   use_bayesian_model_averaging = TRUE,
#'   confidence_level = 0.95
#' )
run_oracle_seq <- function(oracle_seq_dir, data_file, analysis_name, output_dir, thresholds, 
                          bootstrap_samples, use_bayesian_model_averaging, confidence_level) {
  
  # === PARAMETER VALIDATION ===
  cat("ORACLE-seq Analysis Wrapper\n")
  cat("===========================\n\n")
  
  # Check for missing required parameters
  missing_params <- c()
  
  if (missing(oracle_seq_dir) || is.null(oracle_seq_dir)) {
    missing_params <- c(missing_params, "oracle_seq_dir")
  }
  if (missing(data_file) || is.null(data_file)) {
    missing_params <- c(missing_params, "data_file")
  }
  if (missing(analysis_name) || is.null(analysis_name)) {
    missing_params <- c(missing_params, "analysis_name")
  }
  if (missing(output_dir) || is.null(output_dir)) {
    missing_params <- c(missing_params, "output_dir")
  }
  if (missing(thresholds) || is.null(thresholds)) {
    missing_params <- c(missing_params, "thresholds")
  }
  if (missing(bootstrap_samples) || is.null(bootstrap_samples)) {
    missing_params <- c(missing_params, "bootstrap_samples")
  }
  if (missing(use_bayesian_model_averaging) || is.null(use_bayesian_model_averaging)) {
    missing_params <- c(missing_params, "use_bayesian_model_averaging")
  }
  if (missing(confidence_level) || is.null(confidence_level)) {
    missing_params <- c(missing_params, "confidence_level")
  }
  
  # Stop if any parameters are missing
  if (length(missing_params) > 0) {
    cat("ERROR: Missing required parameter(s):\n")
    for (param in missing_params) {
      cat("  -", param, "\n")
    }
    cat("\nPlease provide all required parameters and try again.\n")
    stop("Missing required parameters")
  }
  
  # === INPUT VALIDATION ===
  
  # Check if oracle-seq directory exists and contains required files
  if (!dir.exists(oracle_seq_dir)) {
    cat("ERROR: ORACLE-seq directory not found:", oracle_seq_dir, "\n")
    stop("ORACLE-seq directory not found")
  }
  
  required_files <- c("R/utils.R", "R/config.R", "R/models.R", "R/analysis.R", "R/plots.R", "main.R")
  missing_files <- c()
  for (file in required_files) {
    if (!file.exists(file.path(oracle_seq_dir, file))) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    cat("ERROR: Missing required ORACLE-seq files in", oracle_seq_dir, ":\n")
    for (file in missing_files) {
      cat("  -", file, "\n")
    }
    stop("Missing ORACLE-seq files")
  }
  
  # Check if data file exists
  if (!file.exists(data_file)) {
    cat("ERROR: Data file not found:", data_file, "\n")
    stop("Data file not found")
  }
  
  # Validate parameter types and ranges
  if (!is.character(analysis_name) || length(analysis_name) != 1) {
    cat("ERROR: analysis_name must be a single character string\n")
    stop("Invalid analysis_name")
  }
  
  if (!is.numeric(thresholds) || any(thresholds <= 0)) {
    cat("ERROR: thresholds must be a numeric vector with positive values\n")
    stop("Invalid thresholds")
  }
  
  if (!is.numeric(bootstrap_samples) || bootstrap_samples <= 0 || bootstrap_samples != round(bootstrap_samples)) {
    cat("ERROR: bootstrap_samples must be a positive integer\n")
    stop("Invalid bootstrap_samples")
  }
  
  if (!is.logical(use_bayesian_model_averaging) || length(use_bayesian_model_averaging) != 1) {
    cat("ERROR: use_bayesian_model_averaging must be TRUE or FALSE\n")
    stop("Invalid use_bayesian_model_averaging")
  }
  
  if (!is.numeric(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    cat("ERROR: confidence_level must be between 0 and 1\n")
    stop("Invalid confidence_level")
  }
  
  # === DIRECTORY SETUP ===
  
  # Create base output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created base output directory:", output_dir, "\n")
  }
  
  # Find available directory name (handle duplicates)
  final_analysis_name <- find_available_directory_name(output_dir, analysis_name)
  analysis_dir <- file.path(output_dir, final_analysis_name)
  
  # Create analysis directory structure
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(analysis_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(analysis_dir, "statistical_analysis"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(analysis_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(analysis_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
  
  cat("Analysis directory created:", analysis_dir, "\n\n")
  
  # === CONFIGURATION SETUP ===
  
  # Save analysis configuration
  config_file <- file.path(analysis_dir, "analysis_config.yaml")
  save_analysis_config(config_file, data_file, analysis_name, thresholds, 
                      bootstrap_samples, use_bayesian_model_averaging, confidence_level)
  
  # Update global configuration for the pipeline
  update_pipeline_config(analysis_dir, data_file, thresholds, bootstrap_samples, 
                        use_bayesian_model_averaging, confidence_level)
  
  # === RUN ANALYSIS ===
  
  cat("Starting ORACLE-seq analysis...\n")
  cat("Parameters:\n")
  cat("  Data file:", data_file, "\n")
  cat("  Analysis name:", final_analysis_name, "\n")
  cat("  Thresholds:", paste(thresholds, collapse = ", "), "f/M\n")
  cat("  Bootstrap samples:", bootstrap_samples, "\n")
  cat("  Bayesian Model Averaging:", ifelse(use_bayesian_model_averaging, "ENABLED", "DISABLED"), "\n")
  cat("  Confidence level:", confidence_level, "\n\n")
  
  # Start logging
  log_file <- file.path(analysis_dir, "logs", "analysis_log.txt")
  sink(log_file, append = FALSE, split = TRUE)
  
  tryCatch({
    # Use the provided oracle-seq directory path
    cat("Loading ORACLE-seq modules from:", oracle_seq_dir, "\n")
    
    # Set global variable for main.R to use
    assign("ORACLE_SEQ_DIR", oracle_seq_dir, envir = .GlobalEnv)
    
    # Load required libraries
    source(file.path(oracle_seq_dir, "R", "utils.R"))
    load_libraries()
    
    # Load configuration
    source(file.path(oracle_seq_dir, "R", "config.R"))
    
    # Run the main analysis pipeline
    source(file.path(oracle_seq_dir, "main.R"))
    
    cat("\n", paste(rep("=", 50), collapse = ""), "\n")
    cat("ORACLE-seq analysis completed successfully!\n")
    cat("Results saved to:", analysis_dir, "\n")
    cat(paste(rep("=", 50), collapse = ""), "\n")
    
  }, error = function(e) {
    cat("\nERROR: Analysis failed with error:\n")
    cat(e$message, "\n")
    stop("Analysis failed")
  }, finally = {
    # Stop logging
    sink()
  })
}

#' Find Available Directory Name with Incremental Suffix
#'
#' Finds an available directory name by appending incremental numbers
#' if the base name already exists.
#'
#' @param base_dir Character string specifying base directory path
#' @param desired_name Character string specifying desired directory name
#'
#' @return Character string of available directory name
#'
#' @details
#' If desired_name already exists in base_dir, the function appends
#' "_2", "_3", etc. until an available name is found.
#'
#' @examples
#' # Returns "study_1" if "study_1" doesn't exist
#' # Returns "study_1_2" if "study_1" exists but "study_1_2" doesn't
#' name <- find_available_directory_name("/output", "study_1")
find_available_directory_name <- function(base_dir, desired_name) {
  # Start with the desired name
  candidate_name <- desired_name
  counter <- 2
  
  # Keep incrementing until we find an available name
  while (dir.exists(file.path(base_dir, candidate_name))) {
    candidate_name <- paste0(desired_name, "_", counter)
    counter <- counter + 1
  }
  
  # Inform user if name was changed
  if (candidate_name != desired_name) {
    cat("Directory '", desired_name, "' already exists.\n", sep = "")
    cat("Using '", candidate_name, "' instead.\n", sep = "")
  }
  
  return(candidate_name)
}

#' Save Analysis Configuration to YAML File
#'
#' Saves the analysis parameters to a YAML configuration file for
#' reproducibility and reference.
#'
#' @param config_file Character string specifying output file path
#' @param data_file Character string specifying data file path
#' @param analysis_name Character string specifying analysis name
#' @param thresholds Numeric vector of thresholds
#' @param bootstrap_samples Integer number of bootstrap samples
#' @param use_bayesian_model_averaging Logical indicating BMA usage
#' @param confidence_level Numeric confidence level
#'
#' @return NULL (side effect: creates YAML file)
save_analysis_config <- function(config_file, data_file, analysis_name, 
                                thresholds, bootstrap_samples, use_bayesian_model_averaging, 
                                confidence_level) {
  
  config_content <- paste0(
    "# ORACLE-seq Analysis Configuration\n",
    "# Generated: ", Sys.time(), "\n\n",
    "analysis:\n",
    "  name: \"", analysis_name, "\"\n",
    "  timestamp: \"", Sys.time(), "\"\n\n",
    "data:\n",
    "  rarefaction_file: \"", data_file, "\"\n\n",
    "parameters:\n",
    "  marginal_thresholds: [", paste(thresholds, collapse = ", "), "]\n",
    "  bootstrap_samples: ", bootstrap_samples, "\n",
    "  confidence_level: ", confidence_level, "\n",
    "  use_bayesian_averaging: ", tolower(use_bayesian_model_averaging), "\n\n",
    "r_session:\n",
    "  r_version: \"", R.Version()$version.string, "\"\n",
    "  platform: \"", R.Version()$platform, "\"\n"
  )
  
  writeLines(config_content, config_file)
  cat("Configuration saved to:", config_file, "\n")
}

#' Update Pipeline Configuration Variables
#'
#' Updates the global configuration variables used by the ORACLE-seq
#' pipeline to use the wrapper-provided parameters.
#'
#' @param analysis_dir Character string specifying analysis directory
#' @param data_file Character string specifying data file path
#' @param thresholds Numeric vector of thresholds
#' @param bootstrap_samples Integer number of bootstrap samples
#' @param use_bayesian_model_averaging Logical indicating BMA usage
#' @param confidence_level Numeric confidence level
#'
#' @return NULL (side effect: updates global variables)
update_pipeline_config <- function(analysis_dir, data_file, thresholds, 
                                  bootstrap_samples, use_bayesian_model_averaging, confidence_level) {
  
  # Update directory paths
  assign("BASE_DIR", analysis_dir, envir = .GlobalEnv)
  assign("STATS_DIR", file.path(analysis_dir, "statistical_analysis"), envir = .GlobalEnv)
  assign("PLOTS_DIR", file.path(analysis_dir, "plots"), envir = .GlobalEnv)
  assign("REPORTS_DIR", file.path(analysis_dir, "reports"), envir = .GlobalEnv)
  assign("LOGS_DIR", file.path(analysis_dir, "logs"), envir = .GlobalEnv)
  
  # Update analysis parameters
  assign("MARGINAL_THRESHOLDS", thresholds, envir = .GlobalEnv)
  assign("N_BOOTSTRAP", bootstrap_samples, envir = .GlobalEnv)
  assign("BMA_BOOTSTRAP_SAMPLES", bootstrap_samples, envir = .GlobalEnv)
  assign("USE_BAYESIAN_MODEL_AVERAGING", use_bayesian_model_averaging, envir = .GlobalEnv)
  assign("CONFIDENCE_LEVEL", confidence_level, envir = .GlobalEnv)
  
  # Update data file path (the pipeline will need to load this)
  assign("RAREFACTION_DATA_FILE", data_file, envir = .GlobalEnv)
  
  cat("Pipeline configuration updated for wrapper analysis\n")
}

# Example usage (commented out):
# run_oracle_seq(
#   oracle_seq_dir = "/data/gpfs/projects/punim2251/Aim1_LongBench/oracle-seq-dev",
#   data_file = "data/rarefaction_curves.rds",
#   analysis_name = "lung_cancer_study",
#   output_dir = "oracle_results",
#   thresholds = c(10, 1, 0.1),
#   bootstrap_samples = 1000,
#   use_bayesian_model_averaging = TRUE,
#   confidence_level = 0.95
# ) 