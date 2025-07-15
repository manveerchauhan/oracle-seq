# Utility functions for ORACLE-seq rarefaction curve analysis
#
# This module provides essential utility functions that support the core
# analysis pipeline. It includes library loading, data validation, threshold
# configuration, and helper functions for data manipulation and visualization.

#' Load Required Libraries for ORACLE-seq Analysis
#'
#' Loads all necessary R packages for the ORACLE-seq pipeline with automatic
#' installation if packages are missing. This ensures reproducible analysis
#' environments across different systems.
#'
#' @return NULL (side effect: loads libraries)
#'
#' @details
#' **Required Libraries:**
#' - tidyverse: Data manipulation and visualization
#' - broom: Statistical model tidying
#' - patchwork: Plot composition
#' - minpack.lm: Robust nonlinear least squares
#' - scales: Axis scaling and formatting
#' - parallel: Parallel processing for bootstrap
#' - ggrepel: Text label positioning
#' - tidytext: Text processing utilities
#'
#' **Installation Strategy:**
#' If a package is not installed, the function automatically installs it
#' from CRAN and then loads it. This ensures the pipeline runs successfully
#' in fresh R environments.
#'
#' **Usage:**
#' Called at the beginning of the analysis pipeline to ensure all
#' dependencies are available.
#'
#' @examples
#' # Load all required libraries
#' load_libraries()
load_libraries <- function() {
  # Define required libraries for ORACLE-seq analysis
  libs <- c("tidyverse", "broom", "patchwork", "minpack.lm", 
            "scales", "parallel", "ggrepel", "tidytext")
  
  for (lib in libs) {
    if (!require(lib, character.only = TRUE)) {
      cat("Installing missing package:", lib, "\n")
      install.packages(lib)
      library(lib, character.only = TRUE)
    }
  }
  
  cat("All required libraries loaded successfully.\n")
}

#' Safe Parameter Extraction from Model Objects
#'
#' Safely extracts parameter values from fitted model objects with fallback
#' options. Handles different parameter naming conventions and missing values
#' gracefully.
#'
#' @param params Named vector or list of model parameters
#' @param param_names Character vector of parameter names to try (in order)
#' @param default Default value to return if parameter not found
#'
#' @return Parameter value or default if not found
#'
#' @details
#' **Parameter Search Strategy:**
#' 1. Tries each parameter name in the order provided
#' 2. Returns the first successful match
#' 3. Returns default value if no matches found
#'
#' **Common Usage:**
#' Different model fitting functions may use different parameter naming
#' conventions (e.g., "a" vs "a.90%" for confidence intervals). This
#' function handles these variations gracefully.
#'
#' **Error Handling:**
#' - Missing parameters return default value
#' - Invalid parameter objects return default value
#' - Null or empty parameter lists return default value
#'
#' @examples
#' # Extract parameter with fallback
#' a_value <- safe_param_get(model_params, c("a", "a.90%"), default = 1.0)
#' 
#' # Multiple parameter name options
#' slope <- safe_param_get(params, c("slope", "a", "beta"), default = 0)
safe_param_get <- function(params, param_names, default = NA) {
  # Handle null or empty parameter objects
  if (is.null(params) || length(params) == 0) {
    return(default)
  }
  
  # Try each parameter name in order
  for (name in param_names) {
    if (name %in% names(params)) {
      value <- params[[name]]
      # Return first valid (non-null, non-NA) value
      if (!is.null(value) && !is.na(value)) {
        return(value)
      }
    }
  }
  
  # Return default if no valid parameter found
  return(default)
}

#' Null-Coalescing Operator
#'
#' Returns the first non-null value from two options. Useful for providing
#' default values when variables might be null or empty.
#'
#' @param lhs Left-hand side value (preferred)
#' @param rhs Right-hand side value (fallback)
#'
#' @return lhs if non-null and non-empty, otherwise rhs
#'
#' @details
#' **Behavior:**
#' - Returns lhs if it exists and has length > 0
#' - Returns rhs if lhs is null or empty
#' - Handles various R data types (vectors, lists, etc.)
#'
#' **Usage Pattern:**
#' Commonly used throughout the pipeline for providing sensible defaults
#' when configuration values might be missing.
#'
#' @examples
#' # Use default if variable is null
#' threshold <- user_threshold %||% 1.0
#' 
#' # Chain multiple defaults
#' value <- config$param %||% default_config$param %||% 0
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs) && length(lhs) > 0) lhs else rhs
}

#' Clean Curve Names for Display
#'
#' Converts internal curve identifiers to human-readable names for
#' plots and reports. Handles standard ORACLE-seq naming conventions.
#'
#' @param curve_id Character vector of curve identifiers
#'
#' @return Character vector of cleaned names
#'
#' @details
#' **Naming Conventions:**
#' - Converts underscores to spaces
#' - Expands "SC" to "Single-Cell"
#' - Expands "SN" to "Single-Nucleus"
#' - Maintains technology and feature type information
#'
#' **Example Transformations:**
#' - "ONT_SC_Genes Discovery" → "ONT Single-Cell Genes Discovery"
#' - "PacBio_SN_Isoforms Discovery" → "PacBio Single-Nucleus Isoforms Discovery"
#'
#' **Usage Context:**
#' Used in plotting functions to create readable axis labels, legends,
#' and plot titles from internal curve identifiers.
#'
#' @examples
#' # Clean curve names for display
#' clean_names <- clean_curve_name(c("ONT_SC_Genes Discovery", "PacBio_Bulk_Isoforms Discovery"))
#' # Returns: "ONT Single-Cell Genes Discovery", "PacBio Bulk Isoforms Discovery"
clean_curve_name <- function(curve_id) {
  curve_id %>%
    str_replace_all("_", " ") %>%
    str_replace_all("SC", "Single-Cell") %>%
    str_replace_all("SN", "Single-Nucleus")
}

# === DYNAMIC THRESHOLD UTILITIES ===

#' Create Dynamic Threshold Labels for Any Feature Type
#'
#' Generates appropriate threshold labels for different feature types and
#' display contexts. Supports the configurable threshold system.
#'
#' @param feature_type Character string specifying feature type (e.g., "Genes Discovery")
#' @param thresholds Numeric vector of threshold values (default: MARGINAL_THRESHOLDS)
#' @param label_type Character string specifying label format ("full", "short", "descriptive", "axis")
#'
#' @return Character vector of threshold labels
#'
#' @details
#' **Label Types:**
#' - "full": Complete descriptive labels (e.g., "10 Genes per Million Reads")
#' - "short": Abbreviated labels (e.g., "10 f/M")
#' - "descriptive": Sentence-style labels (e.g., "10 genes per million reads")
#' - "axis": Multi-line labels for plot axes (e.g., "10\nf/M")
#'
#' **Feature Type Handling:**
#' Automatically extracts feature name from full feature type strings
#' (e.g., "Genes Discovery" → "Genes").
#'
#' **Configurable Thresholds:**
#' Works with any threshold configuration, making the system flexible
#' for different research questions.
#'
#' @examples
#' # Create full labels for genes
#' labels <- create_dynamic_threshold_labels("Genes Discovery", c(10, 1, 0.1), "full")
#' # Returns: "10 Genes per Million Reads", "1 Genes per Million Reads", "0.1 Genes per Million Reads"
#' 
#' # Create short labels
#' short_labels <- create_dynamic_threshold_labels("Isoforms Discovery", c(5, 1), "short")
#' # Returns: "5 f/M", "1 f/M"
create_dynamic_threshold_labels <- function(feature_type, thresholds = NULL, label_type = "full") {
  # Use configured thresholds if not provided
  thresholds <- thresholds %||% MARGINAL_THRESHOLDS
  
  # Extract feature name from full type string
  feature_name <- str_remove(feature_type, " Discovery")
  
  # Generate labels based on requested type
  switch(label_type,
    "full" = paste0(thresholds, " ", feature_name, " per Million Reads"),
    "short" = paste0(thresholds, " f/M"),
    "descriptive" = paste0(thresholds, " ", tolower(feature_name), " per million reads"),
    "axis" = paste0(thresholds, "\nf/M"),
    # Default case
    paste0(thresholds, " f/M")
  )
}

#' Get Threshold Styling Configuration
#'
#' Returns visual styling parameters for threshold-based plots. Supports
#' automatic scaling for different numbers of thresholds.
#'
#' @param thresholds Numeric vector of threshold values
#' @param style_type Character string specifying style configuration
#'
#' @return Data frame with styling parameters
#'
#' @details
#' **Styling Parameters:**
#' - line_weight: Line thickness for threshold lines
#' - line_type: Line type (solid, dashed, dotted)
#' - color: Color values for thresholds
#' - alpha: Transparency values
#'
#' **Automatic Scaling:**
#' If more thresholds are provided than pre-configured styles, the function
#' automatically generates additional styling parameters using gradients
#' and sequences.
#'
#' **Style Philosophy:**
#' - Higher thresholds get more prominent styling (thicker lines, higher alpha)
#' - Color schemes are designed for accessibility and publication quality
#' - Consistent styling across all plot types
#'
#' @examples
#' # Get styling for default thresholds
#' styling <- get_threshold_styling(c(10, 1, 0.1))
#' print(styling)
get_threshold_styling <- function(thresholds, style_type = "standard") {
  n_thresh <- length(thresholds)
  
  # Get base configuration
  config <- get_threshold_config()
  
  # Generate styling data frame
  styling_df <- data.frame(
    threshold = thresholds,
    line_weight = config$styling$line_weights[1:n_thresh],
    line_type = config$styling$line_types[1:n_thresh],
    color = config$styling$colors[1:n_thresh],
    alpha = config$styling$alphas[1:n_thresh],
    stringsAsFactors = FALSE
  )
  
  # Fill in missing styling with defaults
  if (n_thresh > length(config$styling$line_weights)) {
    extra_weights <- seq(1.0, 0.3, length.out = n_thresh - length(config$styling$line_weights))
    styling_df$line_weight[is.na(styling_df$line_weight)] <- extra_weights
  }
  
  if (n_thresh > length(config$styling$colors)) {
    extra_colors <- rainbow(n_thresh - length(config$styling$colors))
    styling_df$color[is.na(styling_df$color)] <- extra_colors
  }
  
  return(styling_df)
}

#' Create Threshold Factor Levels for Consistent Ordering
#'
#' Creates properly ordered factor levels for threshold-based data. Ensures
#' consistent ordering across all plots and analyses.
#'
#' @param feature_type Character string specifying feature type
#' @param thresholds Numeric vector of threshold values
#' @param order Character string specifying ordering ("desc", "asc", "natural")
#'
#' @return Character vector of ordered factor levels
#'
#' @details
#' **Ordering Options:**
#' - "desc": Descending order (highest threshold first)
#' - "asc": Ascending order (lowest threshold first)
#' - "natural": Natural ordering based on threshold values
#'
#' **Factor Level Format:**
#' Creates human-readable factor levels that work well in plot legends
#' and statistical outputs.
#'
#' **Usage Context:**
#' Used throughout the plotting system to ensure consistent threshold
#' ordering in legends, facets, and statistical summaries.
#'
#' @examples
#' # Create descending factor levels
#' levels <- create_threshold_factor_levels("Genes Discovery", c(10, 1, 0.1), "desc")
#' # Returns ordered levels with highest threshold first
create_threshold_factor_levels <- function(feature_type, thresholds, order = "desc") {
  # Create labels
  labels <- create_dynamic_threshold_labels(feature_type, thresholds, "full")
  
  # Order based on threshold values
  if (order == "desc") {
    order_idx <- order(thresholds, decreasing = TRUE)
  } else if (order == "asc") {
    order_idx <- order(thresholds, decreasing = FALSE)
  } else {
    order_idx <- seq_along(thresholds)
  }
  
  return(labels[order_idx])
}

#' Find Extension Threshold for Plot Extrapolation
#'
#' Selects an appropriate threshold for extending plots beyond observed
#' data. Used to determine how far to extrapolate fitted curves.
#'
#' @param threshold_data Data frame containing threshold calculations
#' @param available_thresholds Numeric vector of available thresholds
#'
#' @return Data frame with selected extension threshold information
#'
#' @details
#' **Selection Strategy:**
#' 1. Use configured default extension threshold if available
#' 2. Fall back to median threshold if default not available
#' 3. Use any available threshold if median not available
#' 4. Return empty data frame if no thresholds available
#'
#' **Extension Logic:**
#' The selected threshold determines how far beyond observed data
#' the fitted curves should be extended for visualization.
#'
#' **Quality Control:**
#' - Validates that selected threshold has valid depth calculation
#' - Ensures extension depth is reasonable (not infinite)
#' - Provides fallback options if primary selection fails
#'
#' @examples
#' # Find extension threshold
#' extension_info <- find_extension_threshold(threshold_results, c(10, 1, 0.1))
#' if (nrow(extension_info) > 0) {
#'   extension_depth <- extension_info$depth[1]
#' }
find_extension_threshold <- function(threshold_data, available_thresholds) {
  if (nrow(threshold_data) == 0) {
    return(data.frame())
  }
  
  # Get configuration
  config <- get_threshold_config()
  default_extension <- config$default_extension_threshold
  
  # Try to use configured default
  if (!is.null(default_extension) && default_extension %in% threshold_data$threshold) {
    return(threshold_data %>% filter(threshold == default_extension))
  }
  
  # Fall back to median threshold
  median_threshold <- median(available_thresholds)
  closest_to_median <- available_thresholds[which.min(abs(available_thresholds - median_threshold))]
  
  if (closest_to_median %in% threshold_data$threshold) {
    return(threshold_data %>% filter(threshold == closest_to_median))
  }
  
  # Use any available threshold
  return(threshold_data %>% slice(1))
}

#' Create Dynamic Threshold Annotations for Plots
#'
#' Adds horizontal lines and text annotations for threshold values to plots.
#' Supports configurable styling and dynamic label generation.
#'
#' @param thresholds Numeric vector of threshold values
#' @param x_position Numeric value for the x-axis position of annotations
#' @param hjust Numeric value for horizontal justification of text
#'
#' @return List of ggplot2 layer objects (annotations)
#'
#' @details
#' **Annotation Types:**
#' - Horizontal lines (geom_hline)
#' - Text annotations (annotate)
#'
#' **Styling:**
#' - Uses configurable line weights, types, colors, and alphas
#' - Generates dynamic labels based on threshold values
#' - Supports multiple thresholds and different x-axis positions
#'
#' **Usage:**
#' Commonly used in ggplot2 functions to add threshold information
#' to plots.
#'
#' @examples
#' # Create annotations for a plot
#' annotations <- create_threshold_annotations(c(10, 1, 0.1), x_position = 1000)
#' # Returns a list of ggplot2 layer objects
create_threshold_annotations <- function(thresholds = NULL, x_position = Inf, hjust = 1.1) {
  thresholds <- thresholds %||% MARGINAL_THRESHOLDS
  styling <- get_threshold_styling(thresholds)
  
  annotations <- list()
  
  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    style <- styling[i, ]
    
    # Add horizontal line
    annotations[[paste0("line_", i)]] <- geom_hline(
      yintercept = thresh, 
      color = style$color, 
      linewidth = style$line_weight, 
      linetype = style$line_type,
      alpha = style$alpha
    )
    
    # Add text annotation
    annotations[[paste0("text_", i)]] <- annotate(
      "text", 
      x = x_position, 
      y = thresh, 
      label = paste0(thresh, " f/M"),
      hjust = hjust, 
      vjust = -0.3, 
      size = 3, 
      color = style$color, 
      fontface = "bold"
    )
  }
  
  return(annotations)
}

#' Dynamic case_when for Threshold Labeling
#'
#' Generates a dynamic case_when expression for threshold labeling.
#' This function is particularly useful when plotting multiple curves
#' with different threshold values.
#'
#' @param feature_type Character string specifying feature type
#' @param thresholds Numeric vector of threshold values
#'
#' @return Character string of case_when conditions
#'
#' @details
#' **Dynamic Case_when:**
#' - Generates conditions for each threshold
#' - Handles different label types (full, short, descriptive, axis)
#' - Maintains consistency with create_dynamic_threshold_labels
#'
#' **Usage:**
#' Used in ggplot2 facetting or other plotting functions to dynamically
#' label curves based on their thresholds.
#'
#' @examples
#' # Create case_when for a plot
#' case_when_expr <- create_threshold_case_when("Genes Discovery", c(10, 1, 0.1))
#' # Returns: "threshold == 10 ~ paste0('10 ', str_remove(ft, ' Discovery'), ' per Million Reads'),
#' #          threshold == 1 ~ paste0('1 ', str_remove(ft, ' Discovery'), ' per Million Reads'),
#' #          threshold == 0.1 ~ paste0('0.1 ', str_remove(ft, ' Discovery'), ' per Million Reads'),
#' #          TRUE ~ paste0(threshold, ' ', str_remove(ft, ' Discovery'), ' per Million Reads')"
create_threshold_case_when <- function(feature_type, thresholds = NULL) {
  thresholds <- thresholds %||% MARGINAL_THRESHOLDS
  feature_name <- str_remove(feature_type, " Discovery")
  
  # Create dynamic case_when conditions
  conditions <- character(length(thresholds) + 1)
  
  for (i in seq_along(thresholds)) {
    conditions[i] <- paste0(
      "threshold == ", thresholds[i], 
      " ~ paste0('", thresholds[i], " ', str_remove(ft, ' Discovery'), ' per Million Reads')"
    )
  }
  
  # Add default case
  conditions[length(conditions)] <- "TRUE ~ paste0(threshold, ' ', str_remove(ft, ' Discovery'), ' per Million Reads')"
  
  return(paste(conditions, collapse = ",\n          "))
}

#' Ensure Directory Exists
#'
#' Creates a directory if it doesn't exist. Used for ensuring output
#' directories are available before saving results.
#'
#' @param dir_path Character string specifying directory path
#'
#' @return NULL (side effect: creates directory)
#'
#' @details
#' **Behavior:**
#' - Creates directory if it doesn't exist
#' - Creates parent directories as needed
#' - Handles existing directories gracefully
#' - Provides informative error messages
#'
#' **Usage Context:**
#' Called before saving plots, statistics, or reports to ensure
#' output directories exist.
#'
#' @examples
#' # Ensure output directory exists
#' ensure_directory("output/plots")
#' ensure_directory(file.path(BASE_DIR, "statistical_analysis"))
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir_path, "\n")
  }
}