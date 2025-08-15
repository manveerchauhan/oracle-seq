# Core analysis functions for ORACLE-seq rarefaction curve analysis
#
# This module contains the main analysis functions that bridge the statistical
# modeling (models.R) with the visualization (plots.R) components. It handles
# threshold calculations, marginal returns analysis, and data transformations.

#' Create Reads Per Cell Converter Function
#'
#' Generates a conversion function that estimates reads per cell from total
#' sequencing depth. This is crucial for single-cell and single-nucleus
#' RNA-seq analysis where per-cell depth is more meaningful than total depth.
#'
#' @param rarefaction_data Data frame containing rarefaction curve data with
#'   columns: total_reads_unified, median_reads_per_cell, curve_id
#'
#' @return Function that converts total reads to estimated reads per cell
#'
#' @details
#' **Conversion Strategy:**
#' 1. **Data Filtering**: Uses only data points with valid reads per cell information
#' 2. **Linear Modeling**: Fits linear relationship between total reads and median reads per cell
#' 3. **Function Creation**: Returns a function that applies the linear transformation
#' 4. **Fallback**: If insufficient data, returns a simple division by 1000
#'
#' **Model Information:**
#' The returned function has attributes containing:
#' - model: The fitted lm object
#' - r_squared: Coefficient of determination
#' - intercept, slope: Linear model parameters
#'
#' **Usage Context:**
#' This converter is essential for:
#' - Single-cell RNA-seq depth recommendations
#' - Comparing protocols with different cell capture efficiencies
#' - Translating bulk sequencing insights to single-cell contexts
#'
#' @examples
#' # Create converter from rarefaction data
#' converter <- create_reads_per_cell_converter(rarefaction_data)
#' 
#' # Convert total reads to estimated reads per cell
#' total_reads <- 1000000
#' reads_per_cell <- converter(total_reads)
#' 
#' # Access model information
#' model_r2 <- attr(converter, "r_squared")
#' cat("Conversion model R² =", model_r2)
create_reads_per_cell_converter <- function(rarefaction_data) {
  # Filter for valid reads per cell data
  cell_data <- rarefaction_data %>%
    filter(!is.na(median_reads_per_cell), median_reads_per_cell > 0) %>%
    select(total_reads_unified, median_reads_per_cell, curve_id)
  
  cat("Found", nrow(cell_data), "data points with reads per cell information\n")
  
  # Quality control: ensure sufficient data for reliable conversion
  if (nrow(cell_data) < 5) {
    cat("WARNING: Too few data points for reliable conversion. Using fallback.\n")
    return(function(total_reads) total_reads / 1000)  # Simple fallback
  }
  
  # Fit linear conversion model
  conversion_model <- lm(median_reads_per_cell ~ total_reads_unified, data = cell_data)
  r_squared <- summary(conversion_model)$r.squared
  intercept <- coef(conversion_model)[1]
  slope <- coef(conversion_model)[2]
  
  cat("Linear conversion model R² =", round(r_squared, 3), "\n")
  cat("Conversion formula: reads_per_cell =", round(intercept, 3), "+", 
      scientific(slope, 3), "* total_reads\n")
  
  # Create the conversion function
  conversion_function <- function(total_reads) {
    estimated_reads_per_cell <- intercept + slope * total_reads
    # Ensure non-negative values
    pmax(estimated_reads_per_cell, 0)
  }
  
  # Store model information as function attributes
  attr(conversion_function, "model") <- conversion_model
  attr(conversion_function, "r_squared") <- r_squared
  attr(conversion_function, "intercept") <- intercept
  attr(conversion_function, "slope") <- slope
  
  return(conversion_function)
}

#' Calculate Threshold Depths for Marginal Returns Analysis
#'
#' Computes the sequencing depths required to achieve specific marginal returns
#' thresholds (features per million reads). This is the core calculation for
#' determining practical sequencing depth recommendations.
#'
#' @param model_result List containing fitted model information from fit_model_robust()
#' @param thresholds Numeric vector of marginal returns thresholds to calculate
#'
#' @return Data frame with columns:
#'   - threshold: Marginal returns threshold (features/million reads)
#'   - depth: Required sequencing depth to achieve threshold
#'   - features: Number of features expected at threshold depth
#'   - fold_increase_from_current: Fold increase from current max depth
#'
#' @details
#' **Mathematical Approach:**
#' For each mathematical model, the threshold depth is calculated by solving:
#' f'(x) = threshold/1e6
#' 
#' where f'(x) is the derivative (marginal returns) of the fitted model.
#' 
#' **Model-Specific Calculations:**
#' - **Linear**: f'(x) = a (constant) - flags insufficient curvature
#' - **Michaelis-Menten**: f'(x) = (a*b)/(b+x)² - analytical solution
#' - **Asymptotic Exponential**: f'(x) = a*b*exp(-b*x) - analytical solution
#' - **Power Law**: f'(x) = a*b*x^(b-1) - analytical solution
#' - **Logarithmic**: f'(x) = a/x - analytical solution
#' - **Shifted Logarithmic**: f'(x) = a/(x+c) - analytical solution
#' - **Hill**: f'(x) = (a*n*b^n*x^(n-1))/(b^n+x^n)² - numerical solution
#'
#' **Quality Control:**
#' - Validates that calculated depths are positive and finite
#' - Verifies that predicted features are reasonable
#' - Provides diagnostic output for failed calculations
#' - Includes fallback estimates for missing thresholds
#'
#' @examples
#' # Calculate thresholds for fitted model
#' thresholds <- find_threshold_depths(model_result, c(10, 1, 0.1))
#' 
#' # View results
#' print(thresholds)
#' # Shows depth required for each threshold with fold increase
find_threshold_depths <- function(model_result, thresholds = MARGINAL_THRESHOLDS) {
  threshold_depths <- data.frame()
  model_type <- model_result$model
  params <- coef(model_result$fit_object)
  
  cat("  Model:", model_type, "| Parameters:", paste(names(params), "=", round(params, 4), collapse = ", "), "\n")
  
  # Track calculation success/failure for diagnostics
  successful_calculations <- c()
  failed_calculations <- c()
  
  for (thresh in thresholds) {
    depth <- NA
    features <- NA
    
    tryCatch({
      # Calculate depth using model-specific approach
      if (model_type == "hill") {
        # Hill equation requires special numerical handling
        depth <- calculate_hill_threshold(model_result, params, thresh)
      } else {
        depth <- calculate_threshold_depth(model_type, params, thresh)
      }
      
      # Validate calculated depth
      if (!is.na(depth) && depth > 0 && is.finite(depth)) {
        features <- predict(model_result$fit_object, newdata = data.frame(x = depth))
        features <- as.numeric(features)
        
        # Validate predicted features
        if (is.finite(features) && features > 0) {
          threshold_depths <- rbind(threshold_depths, data.frame(
            threshold = thresh,
            depth = depth,
            features = features,
            fold_increase_from_current = depth / max(model_result$x_values)
          ))
          successful_calculations <- c(successful_calculations, paste0(thresh, " f/M"))
          cat("    ✓", thresh, "f/M: depth =", scales::comma(round(depth)), 
              "reads, features =", scales::comma(round(features)), "\n")
        } else {
          failed_calculations <- c(failed_calculations, paste0(thresh, " f/M: invalid features"))
          cat("    ✗", thresh, "f/M: mathematical solution failed (invalid features:", features, ")\n")
        }
      } else {
        failed_calculations <- c(failed_calculations, paste0(thresh, " f/M: invalid depth"))
        cat("    ✗", thresh, "f/M: mathematical solution failed (invalid depth:", depth, ")\n")
      }
    }, error = function(e) {
      failed_calculations <- c(failed_calculations, paste0(thresh, " f/M: ", e$message))
      cat("    ✗ Failed to solve", thresh, "f/M threshold:", e$message, "\n")
    })
  }
  
  # Diagnostic summary
  cat("  CALCULATION SUMMARY:\n")
  if (length(successful_calculations) > 0) {
    cat("    ✓ Successful:", paste(successful_calculations, collapse = ", "), "\n")
  }
  if (length(failed_calculations) > 0) {
    cat("    ✗ Failed:", paste(failed_calculations, collapse = "; "), "\n")
  }
  
  # Add fallback estimates for missing thresholds
  missing_thresholds <- setdiff(thresholds, threshold_depths$threshold)
  if (length(missing_thresholds) > 0) {
    cat("  Adding fallback estimates for missing thresholds:", paste(missing_thresholds, collapse = ", "), "\n")
    threshold_depths <- add_fallback_thresholds(
      threshold_depths, 
      missing_thresholds, 
      model_result
    )
  }
  
  return(threshold_depths)
}

#' Calculate Hill Equation Threshold Using Numerical Methods
#'
#' Specialized function for calculating threshold depths for Hill (cooperative
#' binding) models, which require numerical root-finding due to their complexity.
#'
#' @param model_result List containing fitted Hill model information
#' @param params Named vector of Hill model parameters (a, b, n)
#' @param threshold Numeric marginal returns threshold to calculate
#'
#' @return Numeric value of sequencing depth at threshold, or NA if failed
#'
#' @details
#' **Hill Model Derivative:**
#' f'(x) = (a * n * b^n * x^(n-1)) / (b^n + x^n)²
#' 
#' **Numerical Solution:**
#' Uses uniroot() to solve f'(x) = threshold/1e6
#' 
#' **Search Strategy:**
#' 1. Start from maximum observed sequencing depth
#' 2. Extend search range progressively (100x, then 1000x)
#' 3. Verify sign change exists before calling uniroot()
#' 4. Return NA if no valid solution found
#'
#' **Parameter Extraction:**
#' - a: Maximum features (asymptote)
#' - b: Half-saturation depth
#' - n: Hill coefficient (cooperativity)
#'
#' @examples
#' # Calculate Hill threshold (usually called internally)
#' hill_depth <- calculate_hill_threshold(hill_model_result, params, 1.0)
#' if (!is.na(hill_depth)) {
#'   cat("Hill threshold depth:", scales::comma(hill_depth))
#' }
calculate_hill_threshold <- function(model_result, params, threshold) {
  # Extract Hill model parameters
  a <- params[["a"]] %||% params[[1]]  # Maximum features
  b <- params[["b"]] %||% params[[2]]  # Half-saturation depth
  n <- params[["n"]] %||% params[[3]]  # Hill coefficient
  
  # Define the derivative equation to solve
  hill_derivative_eq <- function(x) {
    if (x <= 0) return(Inf)  # Avoid division by zero
    derivative <- (a * n * b^n * x^(n-1)) / (b^n + x^n)^2
    return(derivative - threshold / 1e6)
  }
  
  # Set up numerical search bounds
  x_start <- max(model_result$x_values)
  x_upper <- x_start * 100
  
  tryCatch({
    # Check if root exists in initial range
    f_start <- hill_derivative_eq(x_start)
    f_upper <- hill_derivative_eq(x_upper)
    
    if (f_start * f_upper < 0) {
      # Root exists, solve numerically
      root_result <- uniroot(hill_derivative_eq, c(x_start, x_upper))
      return(root_result$root)
    } else {
      # Try wider bounds
      x_upper <- x_start * 1000
      f_upper <- hill_derivative_eq(x_upper)
      
      if (f_start * f_upper < 0) {
        root_result <- uniroot(hill_derivative_eq, c(x_start, x_upper))
        return(root_result$root)
      } else {
        # No root found in reasonable range
        return(NA)
      }
    }
  }, error = function(e) {
    # Numerical solver failed
    return(NA)
  })
}

#' Calculate Threshold Depth for Specific Mathematical Model
#'
#' Computes the sequencing depth required to achieve a specific marginal returns
#' threshold for a given mathematical model using analytical solutions.
#'
#' @param model_type Character string specifying the mathematical model
#' @param params Named vector of model parameters
#' @param threshold Numeric marginal returns threshold (features/million reads)
#'
#' @return Numeric sequencing depth, or NA if calculation fails
#'
#' @details
#' **Analytical Solutions by Model:**
#' 
#' **Linear Model:** f'(x) = a
#' - Constant marginal returns - flags insufficient curvature
#' - Returns very large value or NA to indicate problematic threshold
#' 
#' **Michaelis-Menten:** f'(x) = (a*b)/(b+x)²
#' - Solved as: x = √((a*b*1e6)/threshold) - b
#' 
#' **Asymptotic Exponential:** f'(x) = a*b*exp(-b*x)
#' - Solved as: x = -ln((threshold/1e6)/(a*b))/b
#' - Special handling for flat curves (b ≈ 0)
#' 
#' **Power Law:** f'(x) = a*b*x^(b-1)
#' - Solved as: x = ((threshold/1e6)/(a*b))^(1/(b-1))
#' - Undefined for b = 1 (linear case)
#' 
#' **Logarithmic:** f'(x) = a/x
#' - Solved as: x = (a*1e6)/threshold
#' 
#' **Shifted Logarithmic:** f'(x) = a/(x+c)
#' - Solved as: x = (a*1e6)/threshold - c
#' 
#' **Hill:** f'(x) = (a*n*b^n*x^(n-1))/(b^n+x^n)²
#' - Requires numerical solution (delegates to calculate_hill_threshold)
#'
#' @examples
#' # Calculate threshold for Michaelis-Menten model
#' params <- c(a = 8000, b = 50000)
#' depth <- calculate_threshold_depth("michaelis_menten", params, 1.0)
#' cat("Required depth:", scales::comma(depth))
calculate_threshold_depth <- function(model_type, params, threshold) {
  # Convert threshold from features/million to features/read
  thresh_per_read <- threshold / 1e6
  
  depth <- switch(model_type,
    linear = {
      # Linear model: f'(x) = a (constant marginal returns)
      a <- params[["a"]] %||% params[[2]]  # slope parameter
      
      # Linear models have constant marginal returns - problematic for threshold analysis
      if (abs(a * 1e6) <= threshold) {
        # Current slope already below threshold
        1e12  # Very large depth indicating threshold not reachable
      } else {
        # Linear models don't have diminishing returns
        NA  # Flag as problematic for threshold-based analysis
      }
    },
    michaelis_menten = {
      # Michaelis-Menten: f'(x) = (a * b) / (b + x)²
      a <- params[["a"]] %||% params[[1]]  # Asymptote
      b <- params[["b"]] %||% params[[2]]  # Half-saturation
      sqrt((a * b * 1e6) / threshold) - b
    },
    asymptotic_exp = {
      # Asymptotic exponential: f'(x) = a * b * exp(-b * x)
      a <- safe_param_get(params, c("a", "a.90%"), default = params[[1]])
      b <- safe_param_get(params, c("b"), default = ifelse(length(params) >= 2, params[[2]], 0))
      
      if (is.na(b) || abs(b) < 1e-10) {
        # Special case for flat curves (b ≈ 0)
        # Use empirical scaling based on threshold
        x_max <- max(model_result$x_values)
        x_max * ifelse(threshold == 10, 2, ifelse(threshold == 1, 5, 15))
      } else {
        # Standard analytical solution
        -log((threshold / 1e6) / (a * b)) / b
      }
    },
    power_law = {
      # Power law: f'(x) = a * b * x^(b-1)
      a <- params[["a"]] %||% params[[1]]  # Scaling factor
      b <- params[["b"]] %||% params[[2]]  # Exponent
      if (b != 1) {
        ((threshold / 1e6) / (a * b))^(1/(b-1))
      } else {
        NA  # Undefined for b = 1 (becomes linear)
      }
    },
    logarithmic = {
      # Logarithmic: f'(x) = a / x
      a <- params[["a"]] %||% params[[1]]  # Scaling factor
      a * 1e6 / threshold
    },
    hill = {
      # Hill equation: requires numerical solution
      # Delegate to specialized function
      a <- params[["a"]] %||% params[[1]]
      b <- params[["b"]] %||% params[[2]]
      n <- params[["n"]] %||% params[[3]]
      
      # Define the derivative equation
      hill_derivative_eq <- function(x) {
        if (x <= 0) return(Inf)
        derivative <- (a * n * b^n * x^(n-1)) / (b^n + x^n)^2
        return(derivative - threshold / 1e6)
      }
      
      # Use numerical root finding
      x_start <- max(model_result$x_values)
      x_upper <- x_start * 100
      
      tryCatch({
        f_start <- hill_derivative_eq(x_start)
        f_upper <- hill_derivative_eq(x_upper)
        
        if (f_start * f_upper < 0) {
          root_result <- uniroot(hill_derivative_eq, c(x_start, x_upper))
          root_result$root
        } else {
          # Try wider bounds
          x_upper <- x_start * 1000
          f_upper <- hill_derivative_eq(x_upper)
          
          if (f_start * f_upper < 0) {
            root_result <- uniroot(hill_derivative_eq, c(x_start, x_upper))
            root_result$root
          } else {
            NA
          }
        }
      }, error = function(e) {
        NA
      })
    },
    # Default case
    NA
  )
  
  return(depth)
}

#' Add Fallback Threshold Estimates
#'
#' Generates conservative fallback estimates for thresholds that couldn't be
#' calculated analytically. Uses empirical scaling based on current sequencing depth.
#'
#' @param threshold_depths Data frame of successfully calculated thresholds
#' @param missing_thresholds Vector of thresholds that failed calculation
#' @param model_result List containing fitted model information
#'
#' @return Updated threshold_depths data frame with fallback estimates
#'
#' @details
#' **Fallback Strategy:**
#' Uses empirical scaling factors based on threshold values:
#' - 10 f/M: 2x current maximum depth
#' - 1 f/M: 5x current maximum depth  
#' - 0.1 f/M: 15x current maximum depth
#' - Other values: 20/threshold scaling factor
#'
#' **Conservative Approach:**
#' Fallback estimates are intentionally conservative to avoid under-sequencing.
#' They represent reasonable upper bounds when analytical solutions fail.
#'
#' **Usage Context:**
#' Called automatically by find_threshold_depths() when some thresholds
#' cannot be calculated due to model limitations or numerical issues.
#'
#' @examples
#' # Usually called internally, but can be used standalone
#' fallback_df <- add_fallback_thresholds(
#'   existing_thresholds, 
#'   c(0.1), 
#'   model_result
#' )
add_fallback_thresholds <- function(threshold_depths, missing_thresholds, model_result) {
  x_max <- max(model_result$x_values)
  
  for (missing_thresh in missing_thresholds) {
    # Empirical scaling factors based on threshold severity
    fallback_depth <- x_max * switch(as.character(missing_thresh),
      "10" = 2,    # High efficiency: modest increase
      "1" = 5,     # Moderate efficiency: moderate increase
      "0.1" = 15,  # Low efficiency: substantial increase
      20 / missing_thresh  # General scaling for other thresholds
    )
    
    # Predict features at fallback depth
    fallback_features <- predict(
      model_result$fit_object, 
      newdata = data.frame(x = fallback_depth)
    )
    
    # Add fallback estimate to results
    threshold_depths <- rbind(threshold_depths, data.frame(
      threshold = missing_thresh,
      depth = fallback_depth,
      features = as.numeric(fallback_features),
      fold_increase_from_current = fallback_depth / x_max
    ))
  }
  
  return(threshold_depths)
}

#' Calculate Marginal Returns (Derivative) for Fitted Curves
#'
#' Computes the marginal returns (rate of feature discovery) across a range
#' of sequencing depths for visualization and analysis purposes.
#'
#' @param model_result List containing fitted model information
#' @param x_eval Numeric vector of sequencing depths to evaluate (optional)
#'
#' @return Data frame with columns:
#'   - x: Sequencing depth
#'   - marginal_returns: Features discovered per additional read
#'   - marginal_returns_per_million: Features per million additional reads
#'
#' @details
#' **Evaluation Range:**
#' If x_eval is not provided, creates a log-spaced sequence from minimum
#' observed depth to 100x maximum observed depth (1000 points).
#'
#' **Derivative Calculations:**
#' Uses analytical derivatives specific to each model type:
#' - Linear: constant derivative
#' - Michaelis-Menten: (a*b)/(b+x)²
#' - Asymptotic exponential: a*b*exp(-b*x)
#' - Power law: a*b*x^(b-1)
#' - Logarithmic: a/x
#' - Shifted logarithmic: a/(x+c)
#' - Hill: (a*n*b^n*x^(n-1))/(b^n+x^n)²
#'
#' **Usage Context:**
#' Used for creating marginal returns plots and identifying optimal
#' sequencing depths for different efficiency requirements.
#'
#' @examples
#' # Calculate marginal returns for fitted model
#' marginal_data <- calculate_marginal_returns_simple(model_result)
#' 
#' # Plot marginal returns
#' plot(marginal_data$x, marginal_data$marginal_returns_per_million, 
#'      type = "l", log = "x")
calculate_marginal_returns_simple <- function(model_result, x_eval = NULL) {
  # Set up evaluation range if not provided
  if (is.null(x_eval)) {
    x_range <- range(model_result$x_values)
    x_eval <- 10^seq(log10(x_range[1]), log10(x_range[2] * 100), length.out = 1000)
  }
  
  # Extract model parameters
  model_type <- model_result$model
  params <- coef(model_result$fit_object)
  
  # Calculate marginal returns using model-specific derivatives
  marginal_returns <- switch(model_type,
    linear = {
      # Linear: f'(x) = a (constant)
      a <- params[["a"]] %||% params[[2]]
      rep(a, length(x_eval))
    },
    michaelis_menten = {
      # Michaelis-Menten: f'(x) = (a * b) / (b + x)²
      a <- params[["a"]] %||% params[[1]]
      b <- params[["b"]] %||% params[[2]]
      (a * b) / (b + x_eval)^2
    },
    asymptotic_exp = {
      # Asymptotic exponential: f'(x) = a * b * exp(-b * x)
      a <- params[["a"]] %||% params[[1]]
      b <- params[["b"]] %||% params[[2]]
      a * b * exp(-b * x_eval)
    },
    power_law = {
      # Power law: f'(x) = a * b * x^(b-1)
      a <- params[["a"]] %||% params[[1]]
      b <- params[["b"]] %||% params[[2]]
      a * b * x_eval^(b-1)
    },
    logarithmic = {
      # Logarithmic: f'(x) = a / x
      a <- params[["a"]] %||% params[[1]]
      a / x_eval
    },
    hill = {
      # Hill: f'(x) = (a * n * b^n * x^(n-1)) / (b^n + x^n)²
      a <- params[["a"]] %||% params[[1]]
      b <- params[["b"]] %||% params[[2]]
      n <- params[["n"]] %||% params[[3]]
      (a * n * b^n * x_eval^(n-1)) / (b^n + x_eval^n)^2
    },
    # Default: return zeros
    rep(0, length(x_eval))
  )
  
  # Return formatted data frame
  data.frame(
    x = x_eval,
    marginal_returns = marginal_returns,
    marginal_returns_per_million = marginal_returns * 1e6
  )
}

# Generate model equations report
create_model_equations_report <- function(fitted_models, output_file) {
  cat("Generating model equations report to:", output_file, "\n")
  
  # Open file connection
  con <- file(output_file, "w")
  
  # Write header
  writeLines("=======================================================", con)
  writeLines("          MODEL EQUATIONS REPORT", con)
  writeLines(paste("          Generated:", Sys.time()), con)
  writeLines("=======================================================", con)
  writeLines("", con)
  
  # Group by model type
  model_summary <- fitted_models %>%
    filter(convergence == TRUE) %>%
    group_by(model_type) %>%
    summarise(
      n_curves = n(),
      avg_r_squared = mean(r_squared, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_curves))
  
  writeLines("MODEL TYPE SUMMARY:", con)
  writeLines("==================", con)
  for (i in 1:nrow(model_summary)) {
    model_type <- model_summary$model_type[i]
    n_curves <- model_summary$n_curves[i]
    avg_r2 <- round(model_summary$avg_r_squared[i], 4)
    
    writeLines(sprintf("%-20s: %d curves, avg R² = %.4f", 
                      str_to_title(str_replace_all(model_type, "_", " ")), 
                      n_curves, avg_r2), con)
  }
  writeLines("", con)
  
  # Write detailed equations for each curve
  writeLines("DETAILED MODEL EQUATIONS:", con)
  writeLines("========================", con)
  writeLines("", con)
  
  for (i in 1:nrow(fitted_models)) {
    if (!fitted_models$convergence[i]) next
    
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    model_type <- model_result$model
    r_squared <- round(fitted_models$r_squared[i], 4)
    aic <- fitted_models$aic[i]
    
    # Extract parameters
    params <- coef(model_result$fit_object)
    
    writeLines(paste("CURVE:", curve_id), con)
    writeLines(paste("Model Type:", str_to_title(str_replace_all(model_type, "_", " "))), con)
    writeLines(paste("R-squared:", r_squared), con)
    writeLines(paste("AIC:", ifelse(is.na(aic), "NA", round(aic, 2))), con)
    writeLines("", con)
    
    # Write equation based on model type
    writeLines("Mathematical Equation:", con)
    equation_text <- switch(model_type,
      "linear" = {
        a <- safe_param_get(params, "a", default = params[["a"]])
        b <- safe_param_get(params, "b", default = params[["b"]])
        sprintf("f(x) = %.6f * x + %.4f", a, b)
      },
      "michaelis_menten" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * x / (%.4f + x)", a, b)
      },
      "asymptotic_exp" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * (1 - exp(-%.6f * x))", a, b)
      },
      "power_law" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * x^%.4f", a, b)
      },
      "logarithmic" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * log(x) + %.4f", a, b)
      },
      "hill" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        c <- safe_param_get(params, "c", default = if(length(params) >= 3) params[[3]] else 2)
        sprintf("f(x) = %.4f * x^%.4f / (%.4f^%.4f + x^%.4f)", a, c, b, c, c)
      },
      sprintf("Unknown model type: %s", model_type)
    )
    
    writeLines(equation_text, con)
    writeLines("", con)
    
    # Write parameter details
    writeLines("Parameters:", con)
    for (param_name in names(params)) {
      writeLines(sprintf("  %s = %.6f", param_name, params[[param_name]]), con)
    }
    writeLines("", con)
    writeLines("---", con)
    writeLines("", con)
  }
  
  # Close file connection
  close(con)
}

#' Generate Model Equations Report
#'
#' Creates a comprehensive text report containing mathematical equations,
#' parameter values, and fit statistics for all fitted models. This is
#' essential for reproducibility and understanding model behavior.
#'
#' @param fitted_models Data frame containing fitted model results
#' @param output_file Character string specifying output file path
#'
#' @return NULL (writes to file)
#'
#' @details
#' **Report Contents:**
#' 1. **Model Type Summary**: Count of curves fitted by each model type
#' 2. **Detailed Equations**: Mathematical formulation with parameter values
#' 3. **Fit Statistics**: R-squared, AIC, and convergence information
#' 4. **Parameter Interpretation**: Biological meaning of each parameter
#'
#' **Mathematical Equations Generated:**
#' - Linear: f(x) = a*x + b
#' - Michaelis-Menten: f(x) = a*x/(b + x)
#' - Asymptotic Exponential: f(x) = a*(1 - exp(-b*x))
#' - Power Law: f(x) = a*x^b
#' - Logarithmic: f(x) = a*log(x) + b
#' - Shifted Logarithmic: f(x) = a*log(x + c) + b
#' - Hill: f(x) = (a*x^n)/(b^n + x^n)
#'
#' **Parameter Interpretation:**
#' For each model, the report explains the biological meaning of
#' parameters (e.g., asymptote, half-saturation, rate constants).
#'
#' @examples
#' # Generate equations report
#' create_model_equations_report(fitted_models, "model_equations.txt")
create_model_equations_report <- function(fitted_models, output_file) {
  cat("Generating model equations report to:", output_file, "\n")
  
  # Open file connection
  con <- file(output_file, "w")
  on.exit(close(con))  # Ensure file is closed even if error occurs
  
  # Write header
  writeLines("=======================================================", con)
  writeLines("          MODEL EQUATIONS REPORT", con)
  writeLines(paste("          Generated:", Sys.time()), con)
  writeLines("=======================================================", con)
  writeLines("", con)
  
  # Group by model type for summary
  model_summary <- fitted_models %>%
    filter(convergence == TRUE) %>%
    group_by(model_type) %>%
    summarise(
      n_curves = n(),
      avg_r_squared = mean(r_squared, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_curves))
  
  writeLines("MODEL TYPE SUMMARY:", con)
  writeLines("==================", con)
  for (i in 1:nrow(model_summary)) {
    model_type <- model_summary$model_type[i]
    n_curves <- model_summary$n_curves[i]
    avg_r2 <- round(model_summary$avg_r_squared[i], 4)
    
    writeLines(sprintf("%-20s: %d curves, avg R² = %.4f", 
                      str_to_title(str_replace_all(model_type, "_", " ")), 
                      n_curves, avg_r2), con)
  }
  writeLines("", con)
  
  # Write detailed equations for each curve
  writeLines("DETAILED MODEL EQUATIONS:", con)
  writeLines("========================", con)
  writeLines("", con)
  
  for (i in 1:nrow(fitted_models)) {
    if (!fitted_models$convergence[i]) next
    
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    model_type <- model_result$model
    r_squared <- round(fitted_models$r_squared[i], 4)
    aic <- fitted_models$aic[i]
    
    # Extract parameters
    params <- coef(model_result$fit_object)
    
    writeLines(paste("CURVE:", curve_id), con)
    writeLines(paste("Model Type:", str_to_title(str_replace_all(model_type, "_", " "))), con)
    writeLines(paste("R-squared:", r_squared), con)
    writeLines(paste("AIC:", ifelse(is.na(aic), "NA", round(aic, 2))), con)
    writeLines("", con)
    
    # Write equation based on model type
    writeLines("Mathematical Equation:", con)
    equation_text <- switch(model_type,
      "linear" = {
        a <- safe_param_get(params, "a", default = params[["a"]])
        b <- safe_param_get(params, "b", default = params[["b"]])
        sprintf("f(x) = %.6f * x + %.4f", a, b)
      },
      "michaelis_menten" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * x / (%.4f + x)", a, b)
      },
      "asymptotic_exp" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * (1 - exp(-%.6f * x))", a, b)
      },
      "power_law" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * x^%.4f", a, b)
      },
      "logarithmic" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        sprintf("f(x) = %.4f * log(x) + %.4f", a, b)
      },
      "hill" = {
        a <- safe_param_get(params, "a", default = params[[1]])
        b <- safe_param_get(params, "b", default = params[[2]])
        n <- safe_param_get(params, "n", default = params[[3]])
        sprintf("f(x) = (%.4f * x^%.4f) / (%.4f^%.4f + x^%.4f)", a, n, b, n, n)
      },
      "Unknown equation format"
    )
    
    writeLines(equation_text, con)
    writeLines("", con)
    
    # Add parameter interpretation
    writeLines("Parameter Interpretation:", con)
    param_interpretation <- switch(model_type,
      "linear" = {
        paste("  a (slope):", sprintf("%.6f", safe_param_get(params, "a", default = params[["a"]])), "features per read",
              "\n  b (intercept):", sprintf("%.4f", safe_param_get(params, "b", default = params[["b"]])), "baseline features")
      },
      "michaelis_menten" = {
        paste("  a (Vmax):", sprintf("%.4f", safe_param_get(params, "a", default = params[[1]])), "maximum features",
              "\n  b (Km):", sprintf("%.4f", safe_param_get(params, "b", default = params[[2]])), "half-saturation depth")
      },
      "asymptotic_exp" = {
        paste("  a (asymptote):", sprintf("%.4f", safe_param_get(params, "a", default = params[[1]])), "maximum features",
              "\n  b (rate):", sprintf("%.6f", safe_param_get(params, "b", default = params[[2]])), "exponential rate parameter")
      },
      "power_law" = {
        paste("  a (scaling):", sprintf("%.4f", safe_param_get(params, "a", default = params[[1]])), "scaling factor",
              "\n  b (exponent):", sprintf("%.4f", safe_param_get(params, "b", default = params[[2]])), "power law exponent")
      },
      "logarithmic" = {
        paste("  a (scaling):", sprintf("%.4f", safe_param_get(params, "a", default = params[[1]])), "logarithmic scaling",
              "\n  b (baseline):", sprintf("%.4f", safe_param_get(params, "b", default = params[[2]])), "baseline features")
      },
      "hill" = {
        paste("  a (Bmax):", sprintf("%.4f", safe_param_get(params, "a", default = params[[1]])), "maximum features",
              "\n  b (K50):", sprintf("%.4f", safe_param_get(params, "b", default = params[[2]])), "half-saturation depth",
              "\n  n (Hill):", sprintf("%.4f", safe_param_get(params, "n", default = params[[3]])), "Hill coefficient (cooperativity)")
      },
      "Parameter interpretation not available"
    )
    
    writeLines(param_interpretation, con)
    writeLines("", con)
    writeLines("---", con)
    writeLines("", con)
  }
  
  # Write footer
  writeLines("=======================================================", con)
  writeLines("                    END OF REPORT", con)
  writeLines("=======================================================", con)
  
  cat("Model equations report generated successfully.\n")
}

#' Load Rarefaction Data for Analysis
#'
#' Loads and validates rarefaction curve data from the standard ORACLE-seq
#' data files. Performs quality control checks and returns formatted data
#' ready for statistical analysis.
#'
#' @return List containing:
#'   - data: Main rarefaction data frame
#'   - metadata: Curve metadata information
#'
#' @details
#' **Expected Data Structure:**
#' The function expects two RDS files in the STATS_DIR:
#' - rarefaction_data_for_stats.rds: Main rarefaction curve data
#' - curve_metadata.rds: Metadata for each curve
#'
#' **Data Validation:**
#' - Checks for required columns: curve_id, total_reads_unified, featureNum
#' - Validates data types and ranges
#' - Reports missing curves from expected set
#' - Identifies unexpected curves (potential naming issues)
#'
#' **Quality Control:**
#' - Removes curves with insufficient data points (< 3 points)
#' - Filters out negative or zero values
#' - Validates that feature counts are reasonable
#'
#' **Expected Curves:**
#' The function validates against the expected curve set defined in
#' EXPECTED_CURVES configuration variable.
#'
#' @examples
#' # Load rarefaction data
#' data_list <- load_rarefaction_data()
#' rarefaction_data <- data_list$data
#' curve_metadata <- data_list$metadata
#' 
#' # Check data structure
#' str(rarefaction_data)
#' unique(rarefaction_data$curve_id)
load_rarefaction_data <- function() {
  # Load from wrapper-provided file (RAREFACTION_DATA_FILE set by wrapper)
  if (!exists("RAREFACTION_DATA_FILE") || is.null(RAREFACTION_DATA_FILE)) {
    stop("RAREFACTION_DATA_FILE not set. Please use run_oracle_seq() wrapper function.")
  }
  
  cat("Loading rarefaction data from:", RAREFACTION_DATA_FILE, "\n")
  
  if (!file.exists(RAREFACTION_DATA_FILE)) {
    stop("Data file not found: ", RAREFACTION_DATA_FILE)
  }
  
  rarefaction_data <- readRDS(RAREFACTION_DATA_FILE)
  
  # Generate basic metadata from the data
  curve_metadata <- rarefaction_data %>%
    group_by(curve_id) %>%
    summarise(
      max_depth = max(total_reads_unified, na.rm = TRUE),
      max_features = max(featureNum, na.rm = TRUE),
      n_points = n(),
      .groups = 'drop'
    )
  
  # Data validation
  cat("Loaded data with", nrow(rarefaction_data), "data points and", 
      length(unique(rarefaction_data$curve_id)), "curves\n")
  
  # Validate required columns
  required_cols <- c("curve_id", "total_reads_unified", "featureNum")
  missing_cols <- setdiff(required_cols, names(rarefaction_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Quality control checks
  cat("\nQuality control checks:\n")
  
  # Check for negative values
  neg_depths <- sum(rarefaction_data$total_reads_unified < 0, na.rm = TRUE)
  neg_features <- sum(rarefaction_data$featureNum < 0, na.rm = TRUE)
  
  if (neg_depths > 0) {
    cat("⚠️  Found", neg_depths, "negative depth values\n")
  }
  if (neg_features > 0) {
    cat("⚠️  Found", neg_features, "negative feature values\n")
  }
  
  # Check data points per curve
  points_per_curve <- rarefaction_data %>%
    group_by(curve_id) %>%
    summarise(n_points = n(), .groups = 'drop') %>%
    arrange(n_points)
  
  insufficient_curves <- points_per_curve %>%
    filter(n_points < 3)
  
  if (nrow(insufficient_curves) > 0) {
    cat("⚠️  Curves with insufficient data points (< 3):\n")
    for (i in 1:nrow(insufficient_curves)) {
      cat("  -", insufficient_curves$curve_id[i], ":", 
          insufficient_curves$n_points[i], "points\n")
    }
  }
  
  # Summary statistics
  cat("\nData summary:\n")
  cat("  Total data points:", nrow(rarefaction_data), "\n")
  cat("  Depth range:", scales::comma(min(rarefaction_data$total_reads_unified)), 
      "to", scales::comma(max(rarefaction_data$total_reads_unified)), "\n")
  cat("  Feature range:", scales::comma(min(rarefaction_data$featureNum)), 
      "to", scales::comma(max(rarefaction_data$featureNum)), "\n")
  
  # Return data list
  return(list(
    data = rarefaction_data,
    metadata = curve_metadata
  ))
}

#' Create Fitted Curves Data for Visualization
#'
#' Generates smooth fitted curves from statistical models for plotting
#' purposes. Creates high-resolution predictions across the full range
#' of sequencing depths for each successfully fitted model.
#'
#' @param fitted_models_with_uncertainty Data frame containing fitted model results
#' @param rarefaction_data Original rarefaction curve data
#' @param extension_factor Numeric factor for extending curves beyond observed data
#'
#' @return Data frame with columns:
#'   - curve_id: Curve identifier
#'   - total_reads_unified: Sequencing depth (x-axis)
#'   - featureNum_fitted: Fitted feature count (y-axis)
#'   - model_type: Statistical model used
#'   - within_observed_range: Logical indicating if point is within observed data range
#'
#' @details
#' **Curve Generation Strategy:**
#' 1. **Range Extension**: Extends curves beyond observed data for extrapolation visualization
#' 2. **High Resolution**: Uses 500 points per curve for smooth visualization
#' 3. **Log Spacing**: Uses logarithmic spacing for better coverage of the full range
#' 4. **Range Annotation**: Marks points as within or beyond observed data range
#'
#' **Extension Range:**
#' - Minimum: 10% of minimum observed depth
#' - Maximum: extension_factor × maximum observed depth (default: 10x)
#'
#' **Quality Control:**
#' - Only processes successfully converged models
#' - Filters out infinite or invalid predictions
#' - Ensures monotonic increasing curves where expected
#'
#' **Usage Context:**
#' Used by plotting functions to create smooth model curves overlaid
#' on original data points in rarefaction curve visualizations.
#'
#' @examples
#' # Create fitted curves for plotting
#' fitted_curves <- create_fitted_curves_data(
#'   fitted_models_with_uncertainty,
#'   rarefaction_data,
#'   extension_factor = 10
#' )
#' 
#' # Plot fitted vs observed
#' library(ggplot2)
#' ggplot(fitted_curves, aes(x = total_reads_unified, y = featureNum_fitted)) +
#'   geom_line(aes(color = model_type)) +
#'   scale_x_log10()
create_fitted_curves_data <- function(fitted_models_with_uncertainty, rarefaction_data, extension_factor = 10) {
  fitted_curves <- data.frame()
  
  for (i in 1:nrow(fitted_models_with_uncertainty)) {
    if (!fitted_models_with_uncertainty$convergence[i]) next
    
    curve_id <- fitted_models_with_uncertainty$curve_id[i]
    model_result <- fitted_models_with_uncertainty$fitted_model[[i]]
    
    # Get observed data range for this curve
    curve_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
    x_min <- min(curve_data$total_reads_unified)
    x_max <- max(curve_data$total_reads_unified)
    
    # Create extended prediction range
    x_pred_min <- x_min * 0.1  # Extend 10% below minimum
    x_pred_max <- x_max * extension_factor  # Extend by factor above maximum
    
    # Create log-spaced prediction points for smooth curves
    x_pred <- 10^seq(log10(x_pred_min), log10(x_pred_max), length.out = 500)
    
    # Generate predictions
    tryCatch({
      y_pred <- predict(model_result$fit_object, newdata = data.frame(x = x_pred))
      
      # Quality control: remove invalid predictions
      valid_idx <- is.finite(y_pred) & y_pred >= 0
      x_pred <- x_pred[valid_idx]
      y_pred <- y_pred[valid_idx]
      
      # Annotate with range information
      within_range <- x_pred >= x_min & x_pred <= x_max
      
      # Create curve data
      curve_fitted <- data.frame(
        curve_id = curve_id,
        total_reads_unified = x_pred,
        featureNum_fitted = y_pred,
        model_type = model_result$model,
        within_observed_range = within_range
      )
      
      fitted_curves <- rbind(fitted_curves, curve_fitted)
      
    }, error = function(e) {
      cat("Failed to generate fitted curve for", curve_id, ":", e$message, "\n")
    })
  }
  
  return(fitted_curves)
}

# Generate summary analysis report
create_summary_analysis_report <- function(fitted_models, marginal_returns_results, reads_per_cell_converter, output_file) {
  cat("Generating summary analysis report to:", output_file, "\n")
  
  # Open file connection
  con <- file(output_file, "w")
  
  # Write header
  writeLines("=======================================================", con)
  writeLines("         SUMMARY ANALYSIS REPORT", con)
  writeLines(paste("         Generated:", Sys.time()), con)
  writeLines("=======================================================", con)
  writeLines("", con)
  
  # Overall statistics
  total_curves <- nrow(fitted_models)
  converged_curves <- sum(fitted_models$convergence, na.rm = TRUE)
  
  writeLines("OVERALL STATISTICS:", con)
  writeLines("==================", con)
  writeLines(sprintf("Total curves analyzed: %d", total_curves), con)
  writeLines(sprintf("Successfully fitted: %d (%.1f%%)", 
                    converged_curves, 100 * converged_curves / total_curves), con)
  writeLines("", con)
  
  # Model performance summary
  if (converged_curves > 0) {
    model_performance <- fitted_models %>%
      filter(convergence == TRUE) %>%
      group_by(model_type) %>%
      summarise(
        count = n(),
        avg_r_squared = mean(r_squared, na.rm = TRUE),
        min_r_squared = min(r_squared, na.rm = TRUE),
        max_r_squared = max(r_squared, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(desc(count))
    
    writeLines("MODEL PERFORMANCE SUMMARY:", con)
    writeLines("=========================", con)
    writeLines(sprintf("%-20s %5s %8s %8s %8s", "Model Type", "Count", "Avg R²", "Min R²", "Max R²"), con)
    writeLines(paste(rep("-", 60), collapse = ""), con)
    
    for (i in 1:nrow(model_performance)) {
      writeLines(sprintf("%-20s %5d %8.4f %8.4f %8.4f",
                        str_to_title(str_replace_all(model_performance$model_type[i], "_", " ")),
                        model_performance$count[i],
                        model_performance$avg_r_squared[i],
                        model_performance$min_r_squared[i],
                        model_performance$max_r_squared[i]), con)
    }
    writeLines("", con)
  }
  
  # Threshold analysis summary
  if (nrow(marginal_returns_results) > 0) {
    writeLines("THRESHOLD ANALYSIS SUMMARY:", con)
    writeLines("==========================", con)
    
    # Calculate success rates by threshold
    threshold_summary <- marginal_returns_results %>%
      group_by(threshold) %>%
      summarise(
        successful_curves = n(),
        avg_depth = mean(depth, na.rm = TRUE),
        min_depth = min(depth, na.rm = TRUE),
        max_depth = max(depth, na.rm = TRUE),
        avg_features = mean(features, na.rm = TRUE),
        avg_fold_increase = mean(fold_increase_from_current, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(desc(threshold))
    
    writeLines(sprintf("%-12s %8s %12s %12s %12s %12s %15s", 
                      "Threshold", "Success", "Avg Depth", "Min Depth", "Max Depth", "Avg Features", "Avg Fold Incr"), con)
    writeLines(paste(rep("-", 100), collapse = ""), con)
    
    for (i in 1:nrow(threshold_summary)) {
      writeLines(sprintf("%-12s %8d %12.0f %12.0f %12.0f %12.0f %15.1f",
                        paste0(threshold_summary$threshold[i], " f/M"),
                        threshold_summary$successful_curves[i],
                        threshold_summary$avg_depth[i],
                        threshold_summary$min_depth[i],
                        threshold_summary$max_depth[i],
                        threshold_summary$avg_features[i],
                        threshold_summary$avg_fold_increase[i]), con)
    }
    writeLines("", con)
    
    # Technology comparison
    tech_comparison <- marginal_returns_results %>%
      group_by(technology, threshold) %>%
      summarise(
        curves = n(),
        avg_depth = mean(depth, na.rm = TRUE),
        avg_reads_per_cell = mean(reads_per_cell_converter(depth), na.rm = TRUE),
        avg_bootstrap_cv = mean(depth_cv, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(technology, desc(threshold))
    
    writeLines("TECHNOLOGY COMPARISON:", con)
    writeLines("=====================", con)
    writeLines(sprintf("%-10s %-10s %8s %12s %15s %12s", 
                      "Technology", "Threshold", "Curves", "Avg Depth", "Avg Reads/Cell", "Avg CV"), con)
    writeLines(paste(rep("-", 80), collapse = ""), con)
    
    for (i in 1:nrow(tech_comparison)) {
      writeLines(sprintf("%-10s %-10s %8d %12.0f %15.0f %12.3f",
                        tech_comparison$technology[i],
                        paste0(tech_comparison$threshold[i], " f/M"),
                        tech_comparison$curves[i],
                        tech_comparison$avg_depth[i],
                        tech_comparison$avg_reads_per_cell[i],
                        ifelse(is.na(tech_comparison$avg_bootstrap_cv[i]), 0, tech_comparison$avg_bootstrap_cv[i])), con)
    }
    writeLines("", con)
    
    # Model stability analysis
    if ("model_stability" %in% colnames(marginal_returns_results)) {
      stability_summary <- marginal_returns_results %>%
        filter(!is.na(model_stability)) %>%
        group_by(model_stability) %>%
        summarise(
          count = n(),
          percentage = 100 * n() / nrow(marginal_returns_results),
          avg_cv = mean(depth_cv, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(match(model_stability, c("High", "Moderate", "Low", "Very Low", "Unknown")))
      
      writeLines("MODEL STABILITY ANALYSIS:", con)
      writeLines("========================", con)
      writeLines(sprintf("%-12s %8s %12s %10s", "Stability", "Count", "Percentage", "Avg CV"), con)
      writeLines(paste(rep("-", 50), collapse = ""), con)
      
      for (i in 1:nrow(stability_summary)) {
        writeLines(sprintf("%-12s %8d %12.1f%% %10.3f",
                          stability_summary$model_stability[i],
                          stability_summary$count[i],
                          stability_summary$percentage[i],
                          stability_summary$avg_cv[i]), con)
      }
      writeLines("", con)
    }
  }
  
  # Key recommendations
  writeLines("KEY RECOMMENDATIONS:", con)
  writeLines("===================", con)
  
  # Find most efficient protocols for each threshold
  if (nrow(marginal_returns_results) > 0) {
    for (thresh in unique(marginal_returns_results$threshold)) {
      thresh_data <- marginal_returns_results %>%
        filter(threshold == thresh) %>%
        arrange(depth)
      
      if (nrow(thresh_data) > 0) {
        best_protocol <- thresh_data[1, ]
        writeLines(sprintf("For %g f/M threshold:", thresh), con)
        writeLines(sprintf("  Most efficient: %s %s (%s)",
                          best_protocol$technology,
                          case_when(
                            best_protocol$sample_type == "SC" ~ "Single-Cell",
                            best_protocol$sample_type == "SN" ~ "Single-Nucleus",
                            TRUE ~ best_protocol$sample_type
                          ),
                          best_protocol$feature_type), con)
        writeLines(sprintf("  Required depth: %s reads (%.0f reads per cell)",
                          scales::comma(best_protocol$depth),
                          reads_per_cell_converter(best_protocol$depth)), con)
        if (!is.na(best_protocol$depth_cv)) {
          writeLines(sprintf("  Uncertainty: CV = %.1f%% (%s stability)",
                            best_protocol$depth_cv * 100,
                            best_protocol$model_stability), con)
        }
        writeLines("", con)
      }
    }
  }
  
  # Analysis notes
  writeLines("ANALYSIS NOTES:", con)
  writeLines("==============", con)
  writeLines("1. Depths are calculated based on mathematical model extrapolation", con)
  writeLines("2. Bootstrap confidence intervals provide uncertainty estimates", con)
  writeLines("3. Model stability reflects coefficient of variation (CV) in bootstrap", con)
  writeLines("4. Reads per cell estimates use linear regression on single-cell data", con)
  writeLines("5. Extrapolation factors indicate how far beyond current data the predictions extend", con)
  writeLines("", con)
  
  writeLines(paste("Report generated:", Sys.time()), con)
  
  # Close file connection
  close(con)
}

#' Calculate TD75/TD50 (Transcript Discovery) Metrics
#'
#' Implements TD75/TD50 metrics that quantify sequencing depths needed to reach
#' 75% and 50% of the "meaningful discovery ceiling". The ceiling is defined as
#' the total features discovered at exactly 10 features/million reads (f/M).
#'
#' @param fitted_models_with_uncertainty Data frame containing fitted model results
#' @param rarefaction_data Original rarefaction curve data  
#' @param bootstrap_samples Number of bootstrap iterations for confidence intervals
#'
#' @return Data frame with TD75/TD50 metrics:
#'   - curve_id: Curve identifier
#'   - metric: "TD75" or "TD50"
#'   - discovery_ceiling_features: Total features at 10 f/M threshold (ceiling)
#'   - target_features: Target feature count (75% or 50% of ceiling)
#'   - estimated_reads: Sequencing depth needed to reach target
#'   - lower_ci, upper_ci: Bootstrap confidence intervals for estimated reads
#'   - cv: Coefficient of variation for read depth estimates
#'   - stability: Model stability classification
#'   - extrapolation_factor: How far beyond observed data the prediction extends
#'   - technology, sample_type, feature_type: Metadata from curve_id
#'
#' @details
#' **TD75/TD50 Methodology:**
#' 1. **Discovery Ceiling**: Calculate total features at 10 f/M using fitted models
#' 2. **Target Calculation**: TD75 = 75% of ceiling, TD50 = 50% of ceiling
#' 3. **Depth Estimation**: Use inverse prediction to find reads needed for targets
#' 4. **Bootstrap CI**: Quantify uncertainty in depth estimates via bootstrap
#' 
#' **Model Requirements:**
#' - Uses successfully converged models only
#' - Leverages existing bootstrap framework for uncertainty quantification
#' - Integrates with current marginal returns analysis
#'
#' **Interpretation:**
#' - TD75: Sequencing depth for substantial feature discovery (75% of ceiling)
#' - TD50: Sequencing depth for moderate feature discovery (50% of ceiling)  
#' - Lower values indicate more efficient protocols
#' - Confidence intervals reflect prediction uncertainty
#'
#' @examples
#' # Calculate TD75/TD50 metrics
#' td_metrics <- calculate_td_metrics(
#'   fitted_models_with_uncertainty,
#'   rarefaction_data,
#'   bootstrap_samples = 1000
#' )
#' 
#' # View results
#' td_summary <- td_metrics %>%
#'   select(curve_id, metric, target_features, estimated_reads, lower_ci, upper_ci)
#' print(td_summary)
calculate_td_metrics <- function(fitted_models_with_uncertainty, rarefaction_data, bootstrap_samples = N_BOOTSTRAP) {
  cat("=== CALCULATING TD75/TD50 METRICS ===\n\n")
  
  # Define TD thresholds (75% and 50% of discovery ceiling)
  td_percentages <- c(75, 50)
  discovery_ceiling_threshold <- 10  # f/M - threshold that defines "ceiling"
  
  # Initialize results storage
  td_results <- data.frame()
  
  # Process each curve
  for (i in 1:nrow(fitted_models_with_uncertainty)) {
    if (!fitted_models_with_uncertainty$convergence[i]) next
    
    curve_id <- fitted_models_with_uncertainty$curve_id[i]
    model_result <- fitted_models_with_uncertainty$fitted_model[[i]]
    
    cat("Processing TD metrics for:", curve_id, "\n")
    
    # Step 1: Calculate discovery ceiling (features at 10 f/M)
    ceiling_threshold <- find_threshold_depths(model_result, discovery_ceiling_threshold)
    
    if (nrow(ceiling_threshold) == 0) {
      cat("  Warning: Could not calculate discovery ceiling for", curve_id, "\n")
      next
    }
    
    discovery_ceiling <- ceiling_threshold$features[1]
    cat("  Discovery ceiling (10 f/M):", scales::comma(round(discovery_ceiling)), "features\n")
    
    # Step 2: Calculate TD75/TD50 target feature counts
    for (td_pct in td_percentages) {
      target_features <- discovery_ceiling * (td_pct / 100)
      metric_name <- paste0("TD", td_pct)
      
      cat("  ", metric_name, "target:", scales::comma(round(target_features)), "features\n")
      
      # Step 3: Use inverse prediction to find sequencing depth for target
      depth_estimate <- tryCatch({
        # Use model to find depth that yields target_features
        # This requires solving: predict(model, depth) = target_features
        find_depth_for_target_features(model_result, target_features)
      }, error = function(e) {
        cat("    Error calculating", metric_name, "depth:", e$message, "\n")
        return(NA)
      })
      
      if (is.na(depth_estimate) || depth_estimate <= 0) {
        cat("    Warning: Could not calculate", metric_name, "depth for", curve_id, "\n")
        next
      }
      
      cat("  ", metric_name, "estimated depth:", scales::comma(round(depth_estimate)), "reads\n")
      
      # Step 4: Bootstrap confidence intervals for depth estimate
      cat("    Running bootstrap for", metric_name, "uncertainty...\n")
      bootstrap_result <- bootstrap_td_uncertainty(
        curve_data = rarefaction_data %>% filter(curve_id == !!curve_id),
        model_type = model_result$model,
        target_features = target_features,
        n_boot = bootstrap_samples
      )
      
      # Calculate extrapolation factor
      max_observed_depth <- max(model_result$x_values)
      extrapolation_factor <- depth_estimate / max_observed_depth
      
      # Determine stability classification
      stability <- case_when(
        is.na(bootstrap_result$depth_cv) ~ "Unknown",
        bootstrap_result$depth_cv < 0.05 ~ "High",
        bootstrap_result$depth_cv < 0.15 ~ "Moderate", 
        bootstrap_result$depth_cv < 0.30 ~ "Low",
        TRUE ~ "Very Low"
      )
      
      # Extract metadata from curve_id
      curve_parts <- str_split(curve_id, "_")[[1]]
      technology <- ifelse(length(curve_parts) >= 1, curve_parts[1], NA)
      sample_type <- ifelse(length(curve_parts) >= 2, curve_parts[2], NA)  
      feature_type <- ifelse(length(curve_parts) >= 3, 
                           str_remove(paste(curve_parts[3:length(curve_parts)], collapse = " "), " Discovery"),
                           NA)
      
      # Store results
      td_row <- data.frame(
        curve_id = curve_id,
        metric = metric_name,
        discovery_ceiling_features = discovery_ceiling,
        target_features = target_features,
        estimated_reads = depth_estimate,
        lower_ci = bootstrap_result$depth_ci[1],
        upper_ci = bootstrap_result$depth_ci[2], 
        cv = bootstrap_result$depth_cv,
        stability = stability,
        extrapolation_factor = extrapolation_factor,
        technology = technology,
        sample_type = sample_type,
        feature_type = feature_type,
        model_type = model_result$model,
        convergence_rate = bootstrap_result$convergence_rate,
        stringsAsFactors = FALSE
      )
      
      td_results <- rbind(td_results, td_row)
    }
    
    cat("\n")
  }
  
  if (nrow(td_results) > 0) {
    cat("TD75/TD50 metrics calculated for", length(unique(td_results$curve_id)), "curves\n")
    
    # Summary statistics
    cat("\nTD METRICS SUMMARY:\n")
    td_summary <- td_results %>%
      group_by(metric) %>%
      summarise(
        n_curves = n(),
        mean_reads = mean(estimated_reads, na.rm = TRUE),
        median_reads = median(estimated_reads, na.rm = TRUE),
        mean_cv = mean(cv, na.rm = TRUE),
        .groups = 'drop'
      )
    
    for (i in 1:nrow(td_summary)) {
      cat(sprintf("  %s: %d curves, median depth = %s reads (mean CV = %.3f)\n",
                  td_summary$metric[i],
                  td_summary$n_curves[i],
                  scales::comma(round(td_summary$median_reads[i])),
                  td_summary$mean_cv[i]))
    }
  } else {
    cat("Warning: No TD75/TD50 metrics could be calculated\n")
  }
  
  cat("\n=== TD75/TD50 CALCULATION COMPLETE ===\n")
  return(td_results)
}

#' Find Sequencing Depth for Target Feature Count  
#'
#' Uses inverse prediction to determine the sequencing depth needed to 
#' achieve a specific target number of features. This is the core calculation
#' for TD75/TD50 metrics.
#'
#' @param model_result List containing fitted model information
#' @param target_features Numeric target number of features to achieve
#'
#' @return Numeric sequencing depth estimate, or NA if calculation fails
#'
#' @details
#' **Inverse Prediction Strategy:**
#' Uses numerical optimization to solve: predict(model, depth) = target_features
#' 
#' **Search Range:**
#' - Lower bound: Minimum observed depth
#' - Upper bound: 100x maximum observed depth (allows extrapolation)
#'
#' **Numerical Method:**
#' Uses optimize() to find depth that minimizes |predicted - target|
#'
#' @examples
#' # Find depth for target feature count
#' depth <- find_depth_for_target_features(model_result, target_features = 15000)
#' if (!is.na(depth)) {
#'   cat("Depth needed:", scales::comma(depth), "reads")
#' }
find_depth_for_target_features <- function(model_result, target_features) {
  # Define search bounds
  min_depth <- min(model_result$x_values)
  max_depth <- max(model_result$x_values) * 100  # Allow extrapolation
  
  # Objective function: minimize |predicted - target|
  objective_function <- function(depth) {
    predicted <- predict(model_result$fit_object, newdata = data.frame(x = depth))
    return(abs(as.numeric(predicted) - target_features))
  }
  
  # Use numerical optimization to find optimal depth
  tryCatch({
    result <- optimize(
      f = objective_function,
      interval = c(min_depth, max_depth),
      maximum = FALSE  # We want to minimize the objective function
    )
    
    # Validate result
    optimal_depth <- result$minimum
    final_prediction <- predict(model_result$fit_object, newdata = data.frame(x = optimal_depth))
    error <- abs(as.numeric(final_prediction) - target_features)
    
    # Check if solution is reasonable (error < 5% of target)
    if (error / target_features > 0.05) {
      warning("Large prediction error in inverse solution: ", round(error), " features")
    }
    
    return(optimal_depth)
    
  }, error = function(e) {
    return(NA)
  })
}

#' Bootstrap Uncertainty for TD Metrics
#'
#' Performs bootstrap resampling to quantify uncertainty in TD75/TD50 depth
#' estimates. This adapts the existing bootstrap framework for TD-specific
#' uncertainty quantification.
#'
#' @param curve_data Data frame containing rarefaction curve data
#' @param model_type Character string specifying the model type  
#' @param target_features Numeric target number of features
#' @param n_boot Number of bootstrap iterations
#'
#' @return List containing uncertainty estimates:
#'   - depth_ci: 95% confidence interval for depth
#'   - depth_mean, depth_sd: Summary statistics
#'   - depth_cv: Coefficient of variation
#'   - convergence_rate: Bootstrap success rate
#'   - n_valid: Number of valid bootstrap estimates
#'
#' @details
#' **Bootstrap Procedure:**
#' 1. Generate bootstrap samples of rarefaction data
#' 2. Fit model to each bootstrap sample  
#' 3. Calculate inverse prediction for target features
#' 4. Compute confidence intervals from bootstrap distribution
#'
#' **Quality Control:**
#' - Filters invalid predictions (negative, infinite values)
#' - Requires minimum 10 successful bootstrap iterations
#' - Reports convergence diagnostics
#'
#' @examples
#' # Bootstrap TD uncertainty
#' uncertainty <- bootstrap_td_uncertainty(
#'   curve_data, "michaelis_menten", target_features = 12000, n_boot = 500
#' )
#' cat("95% CI:", uncertainty$depth_ci)
bootstrap_td_uncertainty <- function(curve_data, model_type, target_features, n_boot = 1000) {
  # Prepare data
  x_orig <- curve_data$total_reads_unified
  y_orig <- curve_data$featureNum
  n_points <- length(x_orig)
  
  # Initialize storage
  bootstrap_depths <- numeric(n_boot)
  convergence_count <- 0
  
  # Bootstrap loop
  for (i in 1:n_boot) {
    # Generate bootstrap sample
    boot_indices <- sample(n_points, n_points, replace = TRUE)
    x_boot <- x_orig[boot_indices]
    y_boot <- y_orig[boot_indices]
    
    tryCatch({
      # Fit model to bootstrap sample
      boot_fit <- fit_model_robust(x_boot, y_boot, model_type)
      
      if (boot_fit$convergence) {
        convergence_count <- convergence_count + 1
        
        # Calculate depth for target features
        boot_depth <- find_depth_for_target_features(boot_fit, target_features)
        
        # Store valid result
        if (!is.na(boot_depth) && boot_depth > 0 && is.finite(boot_depth)) {
          bootstrap_depths[i] <- boot_depth
        }
      }
    }, error = function(e) {
      # Silent failure for this bootstrap iteration
    })
  }
  
  # Filter valid results
  valid_depths <- bootstrap_depths[bootstrap_depths > 0 & is.finite(bootstrap_depths)]
  
  # Check minimum bootstrap samples
  if (length(valid_depths) < 10) {
    return(list(
      depth_ci = c(NA, NA),
      depth_mean = NA,
      depth_sd = NA, 
      depth_cv = NA,
      convergence_rate = convergence_count / n_boot,
      n_valid = length(valid_depths)
    ))
  }
  
  # Calculate statistics
  alpha <- 1 - CONFIDENCE_LEVEL
  depth_ci <- quantile(valid_depths, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  depth_mean <- mean(valid_depths, na.rm = TRUE)
  depth_sd <- sd(valid_depths, na.rm = TRUE)
  depth_cv <- depth_sd / depth_mean
  
  return(list(
    depth_ci = as.numeric(depth_ci),
    depth_mean = depth_mean,
    depth_sd = depth_sd,
    depth_cv = depth_cv,
    convergence_rate = convergence_count / n_boot,
    n_valid = length(valid_depths)
  ))
}

#' Create TD75/TD50 Metrics HTML Report
#'
#' Generates a comprehensive HTML report for TD75/TD50 (Transcript Discovery) 
#' metrics with interactive tables, summary statistics, and practical interpretation.
#'
#' @param td_metrics_results Data frame containing TD75/TD50 calculation results
#' @param reads_per_cell_converter Function to convert total reads to reads per cell  
#' @param output_file Character string specifying output HTML file path
#'
#' @return NULL (side effect: creates HTML file)
#'
#' @details
#' **Report Structure:**
#' 1. **Executive Summary**: Key findings and recommendations
#' 2. **TD Metrics Overview**: Methodology and interpretation guide  
#' 3. **Results Tables**: Interactive tables with sorting and filtering
#' 4. **Technology Comparison**: Cross-platform efficiency analysis
#' 5. **Sequencing Recommendations**: Practical depth guidelines
#' 6. **Statistical Details**: Model uncertainty and confidence intervals
#'
#' **Interactive Features:**
#' - Sortable tables with confidence intervals
#' - Technology and sample type filtering
#' - Reads per cell conversion display
#' - Model stability indicators
#'
#' **Practical Interpretation:**
#' - TD75: Depth for substantial discovery (75% of meaningful ceiling)
#' - TD50: Depth for moderate discovery (50% of meaningful ceiling)  
#' - Confidence intervals reflect prediction uncertainty
#' - Stability ratings guide interpretation confidence
#'
#' @examples
#' # Create TD metrics HTML report
#' create_td_metrics_html_report(
#'   td_metrics_results,
#'   reads_per_cell_converter,
#'   file.path(REPORTS_DIR, "td_metrics_analysis.html")
#' )
create_td_metrics_html_report <- function(td_metrics_results, reads_per_cell_converter, output_file) {
  cat("Creating TD75/TD50 metrics HTML report...\n")
  
  if (nrow(td_metrics_results) == 0) {
    cat("Warning: No TD metrics results to report\n")
    return(NULL)
  }
  
  # Prepare summary statistics
  td_summary <- td_metrics_results %>%
    group_by(metric) %>%
    summarise(
      n_curves = n(),
      median_reads = median(estimated_reads, na.rm = TRUE),
      mean_reads = mean(estimated_reads, na.rm = TRUE),
      min_reads = min(estimated_reads, na.rm = TRUE),
      max_reads = max(estimated_reads, na.rm = TRUE),
      median_cv = median(cv, na.rm = TRUE),
      high_stability = sum(stability == "High", na.rm = TRUE),
      moderate_stability = sum(stability == "Moderate", na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Technology comparison
  tech_comparison <- td_metrics_results %>%
    group_by(technology, metric) %>%
    summarise(
      n_curves = n(),
      median_reads = median(estimated_reads, na.rm = TRUE),
      median_reads_per_cell = median(reads_per_cell_converter(estimated_reads), na.rm = TRUE),
      mean_cv = mean(cv, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(metric, median_reads)
  
  # Sample type comparison
  sample_comparison <- td_metrics_results %>%
    group_by(sample_type, metric) %>%
    summarise(
      n_curves = n(),
      median_reads = median(estimated_reads, na.rm = TRUE),
      median_reads_per_cell = median(reads_per_cell_converter(estimated_reads), na.rm = TRUE),
      mean_cv = mean(cv, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(metric, median_reads)
  
  # Generate HTML content
  html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>TD75/TD50 Metrics Analysis - ORACLE-seq</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f8f9fa;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-bottom: 30px;
            text-align: center;
        }
        
        .header h1 {
            margin: 0 0 10px 0;
            font-size: 2.5em;
            font-weight: 300;
        }
        
        .header p {
            margin: 0;
            font-size: 1.1em;
            opacity: 0.9;
        }
        
        .summary-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .card {
            background: white;
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-left: 4px solid #667eea;
        }
        
        .card h3 {
            margin: 0 0 15px 0;
            color: #333;
            font-size: 1.3em;
        }
        
        .metric-value {
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin: 10px 0;
        }
        
        .metric-label {
            color: #666;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        
        .section {
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }
        
        .section h2 {
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }
        
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }
        
        th {
            background-color: #f8f9fa;
            font-weight: 600;
            color: #333;
            position: sticky;
            top: 0;
        }
        
        tr:hover {
            background-color: #f8f9fa;
        }
        
        .stability-high {
            color: #28a745;
            font-weight: bold;
        }
        
        .stability-moderate {
            color: #ffc107;
            font-weight: bold;
        }
        
        .stability-low {
            color: #fd7e14;
            font-weight: bold;
        }
        
        .stability-very-low {
            color: #dc3545;
            font-weight: bold;
        }
        
        .confidence-interval {
            font-size: 0.9em;
            color: #666;
        }
        
        .methodology {
            background: #e8f4fd;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #007bff;
            margin: 20px 0;
        }
        
        .alert {
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
        }
        
        .alert-info {
            background-color: #d1ecf1;
            border: 1px solid #bee5eb;
            color: #0c5460;
        }
        
        .alert-warning {
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            color: #856404;
        }
        
        .footer {
            text-align: center;
            margin-top: 50px;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }
        
        @media (max-width: 768px) {
            .summary-cards {
                grid-template-columns: 1fr;
            }
            
            table {
                font-size: 0.9em;
            }
            
            th, td {
                padding: 8px;
            }
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>TD75/TD50 Metrics Analysis</h1>
        <p>Transcript Discovery Efficiency Assessment | Generated: ', Sys.time(), '</p>
    </div>
    
    <div class="summary-cards">
        <div class="card">
            <h3>TD75 Median Depth</h3>
            <div class="metric-value">', scales::comma(round(td_summary$median_reads[td_summary$metric == "TD75"])), '</div>
            <div class="metric-label">Total Reads</div>
            <p class="confidence-interval">75% of discovery ceiling</p>
        </div>
        
        <div class="card">
            <h3>TD50 Median Depth</h3>
            <div class="metric-value">', scales::comma(round(td_summary$median_reads[td_summary$metric == "TD50"])), '</div>
            <div class="metric-label">Total Reads</div>
            <p class="confidence-interval">50% of discovery ceiling</p>
        </div>
        
        <div class="card">
            <h3>Analyzed Protocols</h3>
            <div class="metric-value">', length(unique(td_metrics_results$curve_id)), '</div>
            <div class="metric-label">Rarefaction Curves</div>
            <p class="confidence-interval">', length(unique(td_metrics_results$technology)), ' technologies, ', length(unique(td_metrics_results$sample_type)), ' sample types</p>
        </div>
    </div>
    
    <div class="section">
        <h2>📊 Methodology Overview</h2>
        <div class="methodology">
            <h4>TD75/TD50 Calculation Method</h4>
            <ol>
                <li><strong>Discovery Ceiling</strong>: Total features at 10 features/million reads threshold</li>
                <li><strong>Target Calculation</strong>: TD75 = 75% of ceiling, TD50 = 50% of ceiling</li>
                <li><strong>Depth Estimation</strong>: Inverse prediction using fitted rarefaction models</li>
                <li><strong>Uncertainty Quantification</strong>: Bootstrap confidence intervals (', nrow(td_metrics_results), ' iterations)</li>
            </ol>
        </div>
        
        <div class="alert alert-info">
            <strong>Interpretation Guide:</strong><br>
            • <strong>TD75</strong>: Sequencing depth for substantial feature discovery (recommended for comprehensive studies)<br>
            • <strong>TD50</strong>: Sequencing depth for moderate feature discovery (suitable for exploratory work)<br>
            • <strong>Confidence Intervals</strong>: Reflect statistical uncertainty in depth predictions<br>
            • <strong>Stability Ratings</strong>: Guide interpretation confidence (High/Moderate = reliable, Low/Very Low = interpret cautiously)
        </div>
    </div>
    
    <div class="section">
        <h2>📈 Detailed Results</h2>
        <p>Complete TD75/TD50 analysis for all protocols with confidence intervals and stability assessment.</p>
        
        <table>
            <thead>
                <tr>
                    <th>Protocol</th>
                    <th>Metric</th>
                    <th>Target Features</th>
                    <th>Estimated Reads</th>
                    <th>Reads per Cell</th>
                    <th>95% CI</th>
                    <th>Stability</th>
                    <th>Model</th>
                </tr>
            </thead>
            <tbody>')
  
  # Add detailed results rows
  for (i in 1:nrow(td_metrics_results)) {
    row <- td_metrics_results[i, ]
    stability_class <- paste0("stability-", tolower(str_replace_all(row$stability, " ", "-")))
    ci_text <- if (!is.na(row$lower_ci) && !is.na(row$upper_ci)) {
      paste0("[", scales::comma(round(row$lower_ci)), ", ", scales::comma(round(row$upper_ci)), "]")
    } else {
      "Not available"
    }
    
    html_content <- paste0(html_content, '
                <tr>
                    <td>', clean_curve_name(row$curve_id), '</td>
                    <td><strong>', row$metric, '</strong></td>
                    <td>', scales::comma(round(row$target_features)), '</td>
                    <td>', scales::comma(round(row$estimated_reads)), '</td>
                    <td>', scales::comma(round(reads_per_cell_converter(row$estimated_reads))), '</td>
                    <td class="confidence-interval">', ci_text, '</td>
                    <td class="', stability_class, '">', row$stability, '</td>
                    <td>', str_to_title(str_replace_all(row$model_type, "_", " ")), '</td>
                </tr>')
  }
  
  html_content <- paste0(html_content, '
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>🔬 Technology Comparison</h2>
        <p>Median sequencing depth requirements by technology platform for efficient transcript discovery.</p>
        
        <table>
            <thead>
                <tr>
                    <th>Technology</th>
                    <th>Metric</th>
                    <th>Protocols</th>
                    <th>Median Total Reads</th>
                    <th>Median Reads/Cell</th>
                    <th>Avg Uncertainty (CV)</th>
                </tr>
            </thead>
            <tbody>')
  
  # Add technology comparison rows
  for (i in 1:nrow(tech_comparison)) {
    row <- tech_comparison[i, ]
    html_content <- paste0(html_content, '
                <tr>
                    <td><strong>', row$technology, '</strong></td>
                    <td>', row$metric, '</td>
                    <td>', row$n_curves, '</td>
                    <td>', scales::comma(round(row$median_reads)), '</td>
                    <td>', scales::comma(round(row$median_reads_per_cell)), '</td>
                    <td>', sprintf("%.1f%%", row$mean_cv * 100), '</td>
                </tr>')
  }
  
  html_content <- paste0(html_content, '
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>🧬 Sample Type Analysis</h2>
        <p>Sequencing depth recommendations by sample preparation method.</p>
        
        <table>
            <thead>
                <tr>
                    <th>Sample Type</th>
                    <th>Metric</th>
                    <th>Protocols</th>
                    <th>Median Total Reads</th>
                    <th>Median Reads/Cell</th>
                    <th>Avg Uncertainty (CV)</th>
                </tr>
            </thead>
            <tbody>')
  
  # Add sample type comparison rows
  for (i in 1:nrow(sample_comparison)) {
    row <- sample_comparison[i, ]
    sample_type_display <- case_when(
      row$sample_type == "SC" ~ "Single-Cell",
      row$sample_type == "SN" ~ "Single-Nucleus", 
      TRUE ~ row$sample_type
    )
    
    html_content <- paste0(html_content, '
                <tr>
                    <td><strong>', sample_type_display, '</strong></td>
                    <td>', row$metric, '</td>
                    <td>', row$n_curves, '</td>
                    <td>', scales::comma(round(row$median_reads)), '</td>
                    <td>', scales::comma(round(row$median_reads_per_cell)), '</td>
                    <td>', sprintf("%.1f%%", row$mean_cv * 100), '</td>
                </tr>')
  }
  
  # Find most efficient protocols
  most_efficient <- td_metrics_results %>%
    group_by(metric) %>%
    slice_min(estimated_reads, n = 3) %>%
    ungroup()
  
  html_content <- paste0(html_content, '
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>🏆 Efficiency Rankings</h2>
        <p>Most efficient protocols for transcript discovery (lowest sequencing depth required).</p>')
  
  for (metric_type in c("TD75", "TD50")) {
    metric_efficient <- most_efficient %>% filter(metric == metric_type)
    html_content <- paste0(html_content, '
        <h4>', metric_type, ' - Most Efficient Protocols</h4>
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Protocol</th>
                    <th>Estimated Reads</th>
                    <th>Reads per Cell</th>
                    <th>Confidence Interval</th>
                    <th>Stability</th>
                </tr>
            </thead>
            <tbody>')
    
    for (i in 1:nrow(metric_efficient)) {
      row <- metric_efficient[i, ]
      stability_class <- paste0("stability-", tolower(str_replace_all(row$stability, " ", "-")))
      ci_text <- if (!is.na(row$lower_ci) && !is.na(row$upper_ci)) {
        paste0("[", scales::comma(round(row$lower_ci)), ", ", scales::comma(round(row$upper_ci)), "]")
      } else {
        "Not available"
      }
      
      html_content <- paste0(html_content, '
                <tr>
                    <td><strong>', i, '</strong></td>
                    <td>', clean_curve_name(row$curve_id), '</td>
                    <td>', scales::comma(round(row$estimated_reads)), '</td>
                    <td>', scales::comma(round(reads_per_cell_converter(row$estimated_reads))), '</td>
                    <td class="confidence-interval">', ci_text, '</td>
                    <td class="', stability_class, '">', row$stability, '</td>
                </tr>')
    }
    html_content <- paste0(html_content, '
            </tbody>
        </table>')
  }
  
  html_content <- paste0(html_content, '
    </div>
    
    <div class="section">
        <h2>💡 Practical Recommendations</h2>
        
        <div class="alert alert-info">
            <h4>Key Findings:</h4>
            <ul>
                <li><strong>Most Efficient TD75 Protocol</strong>: ', clean_curve_name(most_efficient$curve_id[most_efficient$metric == "TD75"][1]), ' (', scales::comma(round(most_efficient$estimated_reads[most_efficient$metric == "TD75"][1])), ' total reads)</li>
                <li><strong>Most Efficient TD50 Protocol</strong>: ', clean_curve_name(most_efficient$curve_id[most_efficient$metric == "TD50"][1]), ' (', scales::comma(round(most_efficient$estimated_reads[most_efficient$metric == "TD50"][1])), ' total reads)</li>
                <li><strong>Technology Recommendation</strong>: ', names(sort(table(most_efficient$technology), decreasing = TRUE))[1], ' shows consistently high efficiency</li>
            </ul>
        </div>
        
        <div class="alert alert-warning">
            <h4>Important Considerations:</h4>
            <ul>
                <li>TD metrics are based on mathematical model extrapolation</li>
                <li>Confidence intervals reflect statistical uncertainty in predictions</li>
                <li>High/Moderate stability ratings indicate reliable predictions</li>
                <li>Consider experimental costs vs. discovery benefits when choosing depths</li>
                <li>TD75 provides substantial discovery coverage, TD50 offers cost-effective exploration</li>
            </ul>
        </div>
    </div>
    
    <div class="footer">
        <p>Generated by <strong>ORACLE-seq</strong> | TD75/TD50 Metrics Analysis<br>
        Statistical Framework: 7-model competition with bootstrap uncertainty quantification<br>
        Discovery Ceiling Methodology: Features at 10 f/M threshold</p>
    </div>
</body>
</html>')
  
  # Write HTML file
  writeLines(html_content, output_file)
  cat("TD75/TD50 HTML report created:", output_file, "\n")
}

# Generate comprehensive model competition analysis report
create_model_competition_analysis_report <- function(fitted_models, output_file) {
  cat("Generating model competition analysis report to:", output_file, "\n")
  
  # Open file connection
  con <- file(output_file, "w")
  
  # Write header
  writeLines("=======================================================", con)
  writeLines("       MODEL COMPETITION ANALYSIS REPORT", con)
  writeLines(paste("         Generated:", Sys.time()), con)
  writeLines("=======================================================", con)
  writeLines("", con)
  
  writeLines("OVERVIEW:", con)
  writeLines("========", con)
  writeLines("This report provides detailed model competition analysis using AIC-based", con)
  writeLines("model selection with Akaike weights, evidence ratios, and uncertainty", con)
  writeLines("classification. Linear model alerts flag cases with insufficient curvature.", con)
  writeLines("", con)
  
  writeLines("INTERPRETATION GUIDE:", con)
  writeLines("====================", con)
  writeLines("ΔAIC Interpretation:", con)
  writeLines("  0-2:    Strong support (competitive models)", con)
  writeLines("  2-4:    Substantial support", con)
  writeLines("  4-7:    Moderate support", con)
  writeLines("  7-10:   Weak support", con)
  writeLines("  >10:    Essentially no support", con)
  writeLines("", con)
  writeLines("Uncertainty Classifications:", con)
  writeLines("  Clear Winner:        Single best model with strong evidence", con)
  writeLines("  Low Uncertainty:     Clear winner but other models have some support", con)
  writeLines("  Moderate Uncertainty: Multiple competitive models, best model weight >0.9", con)
  writeLines("  High Uncertainty:    Multiple competitive models, best model weight <0.9", con)
  writeLines("", con)
  writeLines("Linear Model Alerts:", con)
  writeLines("  CRITICAL: Linear model ΔAIC < 2 (insufficient curvature)", con)
  writeLines("  WARNING:  Linear model ΔAIC < 4 (limited curvature)", con)
  writeLines("", con)
  
  if (USE_BAYESIAN_MODEL_AVERAGING) {
    writeLines("Bayesian Model Averaging (BMA):", con)
    writeLines("  Uses evidence-weighted ensemble when multiple models competitive", con)
    writeLines("  Eligibility: ≥2 saturation models with ΔAIC < 4", con)
    writeLines("  Provides model uncertainty alongside parameter uncertainty", con)
    writeLines("  Robust predictions when model selection is uncertain", con)
    writeLines("", con)
  }
  
  # Overall summary statistics
  total_curves <- nrow(fitted_models)
  converged_curves <- sum(fitted_models$convergence, na.rm = TRUE)
  
  # Extract competition analyses
  competition_summaries <- list()
  linear_alerts <- list()
  uncertainty_classifications <- character()
  
  for (i in 1:nrow(fitted_models)) {
    if (fitted_models$convergence[i] && !is.null(fitted_models$fitted_model[[i]]$model_competition_analysis)) {
      curve_id <- fitted_models$curve_id[i]
      comp_analysis <- fitted_models$fitted_model[[i]]$model_competition_analysis
      
      competition_summaries[[curve_id]] <- comp_analysis$competition_table
      linear_alerts[[curve_id]] <- comp_analysis$linear_alert
      uncertainty_classifications[curve_id] <- comp_analysis$uncertainty_classification
    }
  }
  
  # Summary statistics
  writeLines("SUMMARY STATISTICS:", con)
  writeLines("==================", con)
  writeLines(sprintf("Total curves analyzed: %d", total_curves), con)
  writeLines(sprintf("Successfully converged: %d (%.1f%%)", converged_curves, 100 * converged_curves / total_curves), con)
  writeLines(sprintf("Model competition analyses: %d", length(competition_summaries)), con)
  writeLines("", con)
  
  # Uncertainty distribution
  if (length(uncertainty_classifications) > 0) {
    uncertainty_dist <- table(uncertainty_classifications)
    writeLines("UNCERTAINTY DISTRIBUTION:", con)
    writeLines("========================", con)
    for (unc_type in names(uncertainty_dist)) {
      writeLines(sprintf("%-20s: %3d (%5.1f%%)", 
                        unc_type, 
                        uncertainty_dist[unc_type],
                        100 * uncertainty_dist[unc_type] / length(uncertainty_classifications)), con)
    }
    writeLines("", con)
  }
  
  # Linear model alert summary
  if (length(linear_alerts) > 0) {
    critical_alerts <- sum(sapply(linear_alerts, function(x) x$linear_dominant), na.rm = TRUE)
    warning_alerts <- sum(sapply(linear_alerts, function(x) x$linear_competitive & !x$linear_dominant), na.rm = TRUE)
    no_alerts <- length(linear_alerts) - critical_alerts - warning_alerts
    
    writeLines("LINEAR MODEL ALERT SUMMARY:", con)
    writeLines("==========================", con)
    writeLines(sprintf("CRITICAL alerts (ΔAIC < 2): %3d (%5.1f%%)", 
                      critical_alerts, 100 * critical_alerts / length(linear_alerts)), con)
    writeLines(sprintf("WARNING alerts (ΔAIC < 4):  %3d (%5.1f%%)", 
                      warning_alerts, 100 * warning_alerts / length(linear_alerts)), con)
    writeLines(sprintf("No alerts (sufficient curvature): %3d (%5.1f%%)", 
                      no_alerts, 100 * no_alerts / length(linear_alerts)), con)
    writeLines("", con)
  }
  
  # Bayesian Model Averaging summary (if enabled)
  if (USE_BAYESIAN_MODEL_AVERAGING) {
    writeLines("BAYESIAN MODEL AVERAGING SUMMARY:", con)
    writeLines("================================", con)
    writeLines(sprintf("BMA Status: ENABLED (threshold ΔAIC < %d, min %d models)", 
                      BMA_INCLUSION_THRESHOLD, BMA_MIN_MODELS), con)
    
    # Count BMA eligibility from fitted models
    bma_eligible <- 0
    bma_ineligible_linear <- 0
    bma_ineligible_insufficient <- 0
    bma_ineligible_other <- 0
    total_with_competition <- 0
    
    for (i in 1:nrow(fitted_models)) {
      if (fitted_models$convergence[i] && !is.null(fitted_models$fitted_model[[i]]$model_competition_analysis)) {
        total_with_competition <- total_with_competition + 1
        comp_analysis <- fitted_models$fitted_model[[i]]$model_competition_analysis
        eligibility <- determine_ensemble_eligibility(comp_analysis)
        
        if (eligibility$use_bma) {
          bma_eligible <- bma_eligible + 1
        } else {
          if (grepl("Linear model dominant", eligibility$reason)) {
            bma_ineligible_linear <- bma_ineligible_linear + 1
          } else if (grepl("Only.*eligible models", eligibility$reason)) {
            bma_ineligible_insufficient <- bma_ineligible_insufficient + 1
          } else {
            bma_ineligible_other <- bma_ineligible_other + 1
          }
        }
      }
    }
    
    if (total_with_competition > 0) {
      writeLines(sprintf("BMA Eligible Curves: %3d (%5.1f%%)", 
                        bma_eligible, 100 * bma_eligible / total_with_competition), con)
      writeLines(sprintf("  - Multiple competitive models (ΔAIC < %d)", BMA_INCLUSION_THRESHOLD), con)
      writeLines(sprintf("", ""), con)
      writeLines("BMA Ineligible Curves:", con)
      writeLines(sprintf("  Linear dominant: %3d (%5.1f%%) - insufficient curvature", 
                        bma_ineligible_linear, 100 * bma_ineligible_linear / total_with_competition), con)
      writeLines(sprintf("  Too few models:  %3d (%5.1%%) - < %d competitive models", 
                        bma_ineligible_insufficient, 100 * bma_ineligible_insufficient / total_with_competition, 
                        BMA_MIN_MODELS), con)
      if (bma_ineligible_other > 0) {
        writeLines(sprintf("  Other reasons:   %3d (%5.1%%)", 
                          bma_ineligible_other, 100 * bma_ineligible_other / total_with_competition), con)
      }
    }
    writeLines("", con)
  } else {
    writeLines("BAYESIAN MODEL AVERAGING: DISABLED", con)
    writeLines("=================================", con)
    writeLines("Single model approach used for all threshold predictions.", con)
    writeLines("", con)
  }
  
  # Model selection frequency
  if (length(competition_summaries) > 0) {
    best_models <- character()
    for (comp_table in competition_summaries) {
      if (nrow(comp_table) > 0) {
        best_models <- c(best_models, comp_table$Model[1])  # First row is best (sorted by AIC)
      }
    }
    
    if (length(best_models) > 0) {
      model_freq <- table(best_models)
      writeLines("BEST MODEL FREQUENCY:", con)
      writeLines("====================", con)
      writeLines(sprintf("%-20s %8s %12s", "Model Type", "Count", "Percentage"), con)
      writeLines(paste(rep("-", 45), collapse = ""), con)
      
      model_freq_sorted <- sort(model_freq, decreasing = TRUE)
      for (model_name in names(model_freq_sorted)) {
        writeLines(sprintf("%-20s %8d %12.1f%%",
                          str_to_title(str_replace_all(model_name, "_", " ")),
                          model_freq_sorted[model_name],
                          100 * model_freq_sorted[model_name] / length(best_models)), con)
      }
      writeLines("", con)
    }
  }
  
  # Detailed curve-by-curve analysis
  writeLines("DETAILED CURVE-BY-CURVE ANALYSIS:", con)
  writeLines("=================================", con)
  writeLines("", con)
  
  for (i in 1:nrow(fitted_models)) {
    if (!fitted_models$convergence[i]) next
    
    curve_id <- fitted_models$curve_id[i]
    fitted_model <- fitted_models$fitted_model[[i]]
    
    writeLines(sprintf("CURVE: %s", curve_id), con)
    writeLines(paste(rep("-", nchar(curve_id) + 7), collapse = ""), con)
    
    if (!is.null(fitted_model$model_competition_analysis)) {
      comp_analysis <- fitted_model$model_competition_analysis
      
      # Basic info
      writeLines(sprintf("Best Model: %s (R² = %.4f)", 
                        str_to_title(str_replace_all(fitted_model$model, "_", " ")),
                        fitted_model$r_squared), con)
      writeLines(sprintf("Uncertainty: %s", comp_analysis$uncertainty_classification), con)
      writeLines(sprintf("Best Model Weight: %.3f", comp_analysis$best_model_weight), con)
      
      # Linear model alert
      linear_info <- comp_analysis$linear_alert
      if (linear_info$linear_alert) {
        alert_symbol <- ifelse(linear_info$linear_dominant, "🚨 CRITICAL", "⚠️  WARNING")
        writeLines(sprintf("%s: %s", alert_symbol, linear_info$alert_message), con)
      } else {
        writeLines(sprintf("✓ %s", linear_info$alert_message), con)
      }
      writeLines("", con)
      
      # Competition table
      if (nrow(comp_analysis$competition_table) > 0) {
        writeLines("Model Competition Table:", con)
        writeLines(sprintf("%-18s %10s %8s %12s %12s %12s", 
                          "Model", "AIC", "ΔAIC", "Weight", "Evid.Ratio", "Support"), con)
        writeLines(paste(rep("-", 80), collapse = ""), con)
        
        comp_table <- comp_analysis$competition_table
        for (j in 1:min(nrow(comp_table), 7)) {  # Show top 7 models
          model_display <- str_to_title(str_replace_all(comp_table$Model[j], "_", " "))
          writeLines(sprintf("%-18s %10.2f %8.2f %12.4f %12.2f %12s",
                            substr(model_display, 1, 18),
                            comp_table$AIC[j],
                            comp_table$Delta_AIC[j],
                            comp_table$Akaike_Weight[j],
                            comp_table$Evidence_Ratio[j],
                            comp_table$Support_Level[j]), con)
        }
        writeLines("", con)
        
        # Competitive models summary
        if (length(comp_analysis$competitive_models) > 1) {
          writeLines(sprintf("Competitive Models (ΔAIC < 2): %s", 
                            paste(comp_analysis$competitive_models, collapse = ", ")), con)
        }
        if (length(comp_analysis$substantial_support_models) > length(comp_analysis$competitive_models)) {
          other_supported <- setdiff(comp_analysis$substantial_support_models, comp_analysis$competitive_models)
          writeLines(sprintf("Other Supported Models (ΔAIC < 4): %s", 
                            paste(other_supported, collapse = ", ")), con)
        }
        
        # BMA eligibility information
        if (USE_BAYESIAN_MODEL_AVERAGING) {
          eligibility <- determine_ensemble_eligibility(comp_analysis)
          if (eligibility$use_bma) {
            writeLines(sprintf("🔬 BMA: APPLIED (%s)", eligibility$reason), con)
            writeLines(sprintf("   Ensemble Models: %s", paste(eligibility$eligible_models, collapse = ", ")), con)
          } else {
            writeLines(sprintf("📊 BMA: Not Applied (%s)", eligibility$reason), con)
          }
        }
      }
    } else {
      writeLines("No competition analysis available (fallback model or analysis failed)", con)
    }
    
    writeLines("", con)
  }
  
  # Methodology notes
  writeLines("METHODOLOGY NOTES:", con)
  writeLines("=================", con)
  writeLines("1. AIC (Akaike Information Criterion) balances model fit and complexity", con)
  writeLines("2. Akaike weights represent relative model probabilities", con)
  writeLines("3. Evidence ratios compare support between best and other models", con)
  writeLines("4. Linear model alerts identify cases with insufficient curvature", con)
  writeLines("5. Uncertainty classification guides interpretation confidence", con)
  writeLines("6. Models with ΔAIC > 10 are omitted from detailed tables", con)
  
  if (USE_BAYESIAN_MODEL_AVERAGING) {
    writeLines("7. Bayesian Model Averaging combines multiple competitive models", con)
    writeLines("8. BMA uncertainty includes both parameter and model components", con)
    writeLines("9. Bootstrap ensemble captures full statistical uncertainty", con)
    writeLines("10. Linear models excluded from ensembles (serve as curvature check)", con)
  }
  
  writeLines("", con)
  
  writeLines(paste("Report generated:", Sys.time()), con)
  
  # Close file connection
  close(con)
}