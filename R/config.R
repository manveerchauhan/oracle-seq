# Configuration for ORACLE-seq Rarefaction Analysis
# 
# This configuration file contains constants and helper functions that support
# the ORACLE-seq pipeline. All analysis parameters are provided by the wrapper
# function run_oracle_seq().

# === STATIC CONFIGURATION ===
# These values don't change between analyses

# Bayesian Model Averaging parameters
BMA_INCLUSION_THRESHOLD <- 4  # ΔAIC threshold for including models in ensemble
BMA_MIN_MODELS <- 2           # Minimum number of models required for ensemble

# Model types supported by ORACLE-seq
MODEL_TYPES <- c("linear", "michaelis_menten", "asymptotic_exp", "power_law", 
                 "logarithmic", "shifted_logarithmic", "hill")

# Color schemes for plots
COLOR_SCHEMES <- list(
  standard = c(PacBio = "#DA1884", ONT = "#003F5C"),
  primary_overlap = c(PacBio = "#E31A1C", ONT = "#1F78B4")
)

# Line and shape mappings for plots
LINETYPE_MAP <- c(Bulk = "solid", SC = "dashed", SN = "dotted")
SHAPE_MAP <- c(SC = 16, SN = 17, Bulk = 15)

# Expected curves for validation (used only when data contains these)
EXPECTED_CURVES <- c(
  "ONT_Bulk_Genes Discovery", "ONT_Bulk_Isoforms Discovery",
  "ONT_SC_Genes Discovery", "ONT_SC_Isoforms Discovery", 
  "ONT_SN_Genes Discovery", "ONT_SN_Isoforms Discovery",
  "PacBio_Bulk_Genes Discovery", "PacBio_Bulk_Isoforms Discovery",
  "PacBio_SC_Genes Discovery", "PacBio_SC_Isoforms Discovery",
  "PacBio_SN_Genes Discovery", "PacBio_SN_Isoforms Discovery"
)

# Random seed for reproducibility
SEED <- 42

# === DYNAMIC THRESHOLD CONFIGURATION ===
# These functions work with any threshold values provided by the wrapper

#' Get Threshold Configuration for Dynamic Styling
#'
#' Creates styling configuration for any set of threshold values.
#' Used by plotting functions to ensure consistent visual styling.
#'
#' @param thresholds Numeric vector of threshold values
#' @return List containing styling configuration
get_threshold_config <- function(thresholds = MARGINAL_THRESHOLDS) {
  n_thresh <- length(thresholds)
  
  config <- list(
    values = thresholds,
    
    # Labels for different contexts
    labels = list(
      short = paste0(thresholds, " f/M"),
      descriptive = paste0("Threshold: ", thresholds, " features/million reads")
    ),
    
    # Visual styling for plots
    styling = list(
      line_weights = seq(1.2, 0.4, length.out = n_thresh),
      line_types = rep(c("solid", "dashed", "dotted"), length.out = n_thresh),
      colors = if (n_thresh <= 3) {
        c("#d73027", "#4575b4", "#91bfdb")[1:n_thresh]
      } else {
        rainbow(n_thresh, start = 0, end = 0.7)
      },
      alphas = seq(0.8, 0.5, length.out = n_thresh)
    ),
    
    # Default threshold for plot extensions (use median)
    default_extension_threshold = median(thresholds),
    
    # Ordering for display (high to low)
    display_order = "desc"
  )
  
  return(config)
}

#' Get Ordered Thresholds for Display
#'
#' Returns thresholds in specified order for consistent display.
#'
#' @param thresholds Numeric vector of threshold values
#' @param order Character string: "desc" (default), "asc", or "natural"
#' @return Ordered numeric vector of thresholds
get_ordered_thresholds <- function(thresholds = MARGINAL_THRESHOLDS, order = "desc") {
  if (order == "desc") {
    return(sort(thresholds, decreasing = TRUE))
  } else if (order == "asc") {
    return(sort(thresholds, decreasing = FALSE))
  } else {
    return(thresholds)
  }
}

#' Create Dynamic Threshold Labels
#'
#' Generates labels for threshold values in different formats.
#'
#' @param feature_type Character string specifying feature type
#' @param thresholds Numeric vector of threshold values
#' @param label_type Character string: "full", "short", or "descriptive"
#' @return Character vector of threshold labels
create_threshold_labels <- function(feature_type, thresholds = MARGINAL_THRESHOLDS, label_type = "full") {
  feature_name <- str_remove(feature_type, " Discovery")
  
  switch(label_type,
    "full" = paste0(thresholds, " ", feature_name, " per Million Reads"),
    "short" = paste0(thresholds, " f/M"),
    "descriptive" = paste0(thresholds, " ", tolower(feature_name), " per million reads"),
    paste0(thresholds, " f/M")  # default
  )
}

#' Get Color Map for Plots
#'
#' Returns color mapping for different plot types.
#'
#' @param plot_type Character string specifying plot type
#' @return Named vector of colors
get_color_map <- function(plot_type = "standard") {
  COLOR_SCHEMES[[plot_type]] %||% COLOR_SCHEMES[["standard"]]
}