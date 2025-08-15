# Model fitting functions for ORACLE-seq rarefaction curve analysis
#
# This module implements the core statistical modeling framework for fitting
# mathematical models to rarefaction curves and quantifying uncertainty.
# 
# ORACLE-seq supports seven different mathematical models to capture the
# diversity of rarefaction curve shapes observed in long-read RNA-seq data:
#
# 1. Linear: y = a*x + b
#    - Simplest model for approximately linear growth phases
#    - Rarely the best fit for mature rarefaction curves
#
# 2. Michaelis-Menten: y = (a*x)/(b + x) 
#    - Classic saturation kinetics model
#    - 'a' = maximum features (asymptote)
#    - 'b' = sequencing depth at half-maximum
#
# 3. Asymptotic Exponential: y = a*(1 - exp(-b*x))
#    - Exponential approach to asymptote
#    - 'a' = maximum features, 'b' = rate parameter
#
# 4. Power Law: y = a*x^b
#    - Scale-free relationship, no built-in saturation
#    - Common in early sequencing phases
#
# 5. Logarithmic: y = a*log(x) + b
#    - Slow, continuous growth without saturation
#    - 'a' = scaling factor, 'b' = baseline
#
# 6. Shifted Logarithmic: y = a*log(x + c) + b
#    - Logarithmic with flexible baseline shift
#    - 'c' parameter handles near-zero x values
#
# 7. Hill (Cooperative Binding): y = (a*x^n)/(b^n + x^n)
#    - Generalization of Michaelis-Menten with cooperativity
#    - 'n' = Hill coefficient (cooperativity parameter)
#    - 'a' = maximum features, 'b' = half-saturation point

#' Generate Starting Values for Nonlinear Model Fitting
#'
#' Creates multiple sets of starting parameter values for robust nonlinear
#' model fitting. Each model type requires different parameter initialization
#' strategies to achieve convergence across diverse rarefaction curve shapes.
#'
#' @param x Numeric vector of sequencing depths (independent variable)
#' @param y Numeric vector of feature counts (dependent variable)
#' @param model_type Character string specifying the model type
#'
#' @return List of starting parameter sets:
#'   - conservative: Lower bounds estimate
#'   - optimistic: Upper bounds estimate  
#'   - data_driven: Empirically derived from data characteristics
#'
#' @details
#' The function uses data-driven heuristics to estimate reasonable starting
#' parameters for each model type:
#' 
#' - For saturation models (Michaelis-Menten, Hill): asymptote estimated
#'   from upper quantiles, half-saturation from empirical midpoint
#' - For exponential models: rate parameter scaled by data range
#' - For power law: exponent initialized around 0.5 (sublinear growth)
#' - For logarithmic: scaling factor from data range over log range
#'
#' Multiple starting value sets improve convergence success rates and
#' help identify global vs. local optima in the parameter space.
#'
#' @examples
#' # Generate starting values for Michaelis-Menten model
#' x <- c(1000, 5000, 10000, 20000, 50000)
#' y <- c(2000, 4500, 6000, 7200, 8500)
#' starts <- generate_starting_values(x, y, "michaelis_menten")
#' # Returns: list(conservative=list(a=6800, b=10000), ...)
generate_starting_values <- function(x, y, model_type) {
  # Calculate data characteristics for parameter estimation
  max_y <- max(y, na.rm = TRUE)
  min_y <- min(y, na.rm = TRUE)
  range_y <- max_y - min_y
  max_x <- max(x, na.rm = TRUE)
  min_x <- min(x, na.rm = TRUE)
  median_x <- median(x, na.rm = TRUE)
  
  starting_values <- switch(model_type,
    linear = {
      # For linear model y = a*x + b
      # Estimate slope from data endpoints and intercept from minimum
      initial_slope <- (max_y - min_y) / (max_x - min_x)
      initial_intercept <- min_y
      list(
        conservative = list(a = initial_slope * 0.8, b = initial_intercept),
        optimistic = list(a = initial_slope * 1.2, b = initial_intercept),
        data_driven = list(a = initial_slope, b = initial_intercept)
      )
    },
    michaelis_menten = {
      # For Michaelis-Menten: y = (a*x)/(b + x)
      # 'a' = asymptote (Vmax), 'b' = half-saturation constant (Km)
      # Find x value where y ≈ max_y/2 for Km estimate
      half_max_x <- x[which.min(abs(y - max_y/2))]
      list(
        conservative = list(a = max_y * 0.8, b = median_x),
        optimistic = list(a = max_y * 1.2, b = max_x * 0.3),
        data_driven = list(a = quantile(y, 0.9), b = half_max_x)
      )
    },
    asymptotic_exp = {
      # For asymptotic exponential: y = a*(1 - exp(-b*x))
      # 'a' = asymptote, 'b' = rate parameter
      # Estimate rate from initial slope and asymptote
      initial_slope <- (y[2] - y[1]) / (x[2] - x[1])
      b_start <- abs(initial_slope / max_y)
      list(
        conservative = list(a = max_y * 0.8, b = b_start * 0.5),
        optimistic = list(a = max_y * 1.2, b = b_start * 2),
        data_driven = list(a = quantile(y, 0.9), b = b_start)
      )
    },
    power_law = {
      # For power law: y = a*x^b
      # 'a' = scaling factor, 'b' = exponent (typically 0.3-0.7 for rarefaction)
      # Most rarefaction curves show sublinear growth (b < 1)
      list(
        conservative = list(a = max_y / (max_x^0.5), b = 0.3),
        optimistic = list(a = max_y / (max_x^0.7), b = 0.7),
        data_driven = list(a = y[1] / x[1], b = 0.5)
      )
    },
    logarithmic = {
      # For logarithmic: y = a*log(x) + b
      # 'a' = scaling factor, 'b' = baseline intercept
      # Estimate scaling from data range over log range
      log_range <- log(max_x) - log(min_x)
      list(
        conservative = list(a = range_y / log_range * 0.8, b = min_y),
        optimistic = list(a = range_y / log_range * 1.2, b = min_y),
        data_driven = list(a = (max_y - min_y) / log_range, b = min_y)
      )
    },
    hill = {
      # For Hill equation: y = (a*x^n)/(b^n + x^n)
      # 'a' = maximum binding, 'b' = half-saturation, 'n' = Hill coefficient
      # Hill coefficient typically 1-3 for biological systems
      half_max_x <- x[which.min(abs(y - max_y/2))]
      list(
        conservative = list(a = max_y * 0.9, b = half_max_x, n = 1),
        optimistic = list(a = max_y * 1.1, b = half_max_x * 0.8, n = 2),
        data_driven = list(a = max_y, b = median_x, n = 1.5)
      )
    }
  )
  
  return(starting_values)
}

#' Get Mathematical Formula for Model Type
#'
#' Returns the nonlinear least squares formula object for each supported
#' model type. These formulas are used by nlsLM() for parameter estimation.
#'
#' @param model_type Character string specifying the model type
#' @return Formula object suitable for nonlinear least squares fitting
#'
#' @details
#' All models are formulated as y ~ f(x, parameters) where:
#' - x = sequencing depth (independent variable)
#' - y = feature count (dependent variable)  
#' - parameters = model-specific coefficients to be estimated
#'
#' The formula objects use R's nonlinear modeling syntax and are designed
#' to work with minpack.lm::nlsLM() for robust parameter estimation.
get_model_formula <- function(model_type) {
  formulas <- list(
    linear = y ~ a * x + b,
    michaelis_menten = y ~ (a * x) / (b + x),
    asymptotic_exp = y ~ a * (1 - exp(-b * x)),
    power_law = y ~ a * x^b,
    logarithmic = y ~ a * log(x) + b,
    hill = y ~ (a * x^n) / (b^n + x^n)
  )
  return(formulas[[model_type]])
}

#' Robust Model Fitting with Multiple Starting Values
#'
#' Fits a specified mathematical model to rarefaction curve data using
#' multiple starting parameter sets to maximize convergence success.
#' Handles both linear and nonlinear models with comprehensive error checking.
#'
#' @param x Numeric vector of sequencing depths
#' @param y Numeric vector of feature counts
#' @param model_type Character string specifying model type
#' @param curve_name Character string for identification (optional)
#'
#' @return List containing:
#'   - model: Model type that was fitted
#'   - convergence: Logical indicating successful convergence
#'   - fit_object: Fitted model object (lm or nls)
#'   - r_squared: Coefficient of determination
#'   - aic: Akaike Information Criterion
#'   - parameters: Parameter estimates with standard errors
#'   - predictions: Fitted values
#'   - residuals: Residual values
#'   - x_values, y_values: Original data
#'   - Individual parameter values (a, b, c, n as applicable)
#'
#' @details
#' The function implements a robust fitting strategy:
#' 
#' 1. **Linear Model**: Uses simple lm() for computational efficiency
#' 2. **Nonlinear Models**: Uses nlsLM() with multiple starting values
#' 3. **Data Validation**: Checks for positive values where required
#' 4. **Model Selection**: Chooses best fit among starting value sets by AIC
#' 5. **Error Handling**: Graceful failure with diagnostic information
#'
#' For models requiring positive values (power_law, logarithmic), the function
#' automatically filters invalid data points and requires at least 80% valid
#' data for fitting.
#'
#' @examples
#' # Fit Michaelis-Menten model to rarefaction data
#' depths <- c(1000, 5000, 10000, 20000, 50000)
#' features <- c(2000, 4500, 6000, 7200, 8500)
#' result <- fit_model_robust(depths, features, "michaelis_menten")
#' if (result$convergence) {
#'   cat("R² =", result$r_squared)
#'   cat("AIC =", result$aic)
#' }
fit_model_robust <- function(x, y, model_type, curve_name = "unknown") {
  # Special handling for linear model - use lm() for efficiency and stability
  if (model_type == "linear") {
    tryCatch({
      # Use ordinary least squares for linear model
      lm_fit <- lm(y ~ x)
      
      # Extract results in consistent format with nonlinear models
      coef_values <- coef(lm_fit)
      names(coef_values) <- c("b", "a")  # Intercept first, then slope
      
      # Calculate predictions and fit statistics
      pred_y <- predict(lm_fit)
      residuals_vals <- residuals(lm_fit)
      r_squared <- summary(lm_fit)$r.squared
      aic_value <- AIC(lm_fit)
      
      result <- list(
        model = model_type,
        curve_name = curve_name,
        fit_object = lm_fit,
        parameters = broom::tidy(lm_fit),
        fit_stats = broom::glance(lm_fit),
        r_squared = r_squared,
        aic = aic_value,
        convergence = TRUE,
        predictions = pred_y,
        residuals = residuals_vals,
        x_values = x,
        y_values = y,
        a = coef_values[["a"]],
        b = coef_values[["b"]]
      )
      
      return(result)
      
    }, error = function(e) {
      return(list(
        model = model_type,
        curve_name = curve_name,
        convergence = FALSE,
        aic = NA,
        error = paste("Linear model failed:", e$message)
      ))
    })
  }
  
  # Nonlinear model fitting with multiple starting values
  starting_sets <- generate_starting_values(x, y, model_type)
  best_fit <- NULL
  best_aic <- Inf
  
  formula <- get_model_formula(model_type)
  
  # Data validation for models requiring positive values
  if (model_type %in% c("power_law", "logarithmic")) {
    valid_indices <- x > 0 & y > 0
    if (sum(valid_indices) < length(x) * 0.8) {
      return(list(
        model = model_type,
        curve_name = curve_name,
        convergence = FALSE,
        aic = NA,
        error = "Too many zero/negative values for logarithmic/power models"
      ))
    }
    x_fit <- x[valid_indices]
    y_fit <- y[valid_indices]
  } else {
    x_fit <- x
    y_fit <- y
  }
  
  # Try multiple starting value sets to find global optimum
  for (start_name in names(starting_sets)) {
    start_params <- starting_sets[[start_name]]
    
    tryCatch({
      # Use Levenberg-Marquardt algorithm for robust nonlinear fitting
      fit <- nlsLM(formula, 
                   start = start_params,
                   data = data.frame(x = x_fit, y = y_fit),
                   control = nls.lm.control(maxiter = 1000))
      
      # Select best fit based on AIC (penalizes overfitting)
      aic_value <- AIC(fit)
      
      if (aic_value < best_aic) {
        best_aic <- aic_value
        best_fit <- fit
      }
    }, error = function(e) {
      # Continue to next starting value set if this one fails
    })
  }
  
  # Check if any starting value set succeeded
  if (is.null(best_fit)) {
    return(list(
      model = model_type,
      curve_name = curve_name,
      convergence = FALSE,
      aic = NA,
      error = "All starting value sets failed to converge"
    ))
  }
  
  # Extract results from best fit
  params <- tidy(best_fit)
  fit_stats <- glance(best_fit)
  coef_values <- coef(best_fit)
  
  # Calculate predictions and goodness of fit
  pred_y <- predict(best_fit, newdata = data.frame(x = x))
  residuals <- y - pred_y
  r_squared <- 1 - sum(residuals^2, na.rm = TRUE) / sum((y - mean(y))^2)
  
  result <- list(
    model = model_type,
    curve_name = curve_name,
    fit_object = best_fit,
    parameters = params,
    fit_stats = fit_stats,
    r_squared = r_squared,
    aic = best_aic,
    convergence = TRUE,
    predictions = pred_y,
    residuals = residuals,
    x_values = x,
    y_values = y
  )
  
  # Add individual parameter values for easy access
  for (param_name in names(coef_values)) {
    result[[param_name]] <- coef_values[[param_name]]
  }
  
  return(result)
}

#' Bootstrap Uncertainty Quantification for Threshold Predictions
#'
#' Performs bootstrap resampling to quantify uncertainty in threshold depth
#' predictions. This is crucial for understanding the reliability of 
#' extrapolated sequencing depth recommendations.
#'
#' @param curve_data Data frame containing rarefaction curve data
#' @param model_type Character string specifying the model type
#' @param threshold Numeric value of the marginal returns threshold (features/million reads)
#' @param n_boot Number of bootstrap iterations (default: 1000)
#'
#' @return List containing:
#'   - depth_ci: 95% confidence interval for threshold depth
#'   - features_ci: 95% confidence interval for features at threshold
#'   - depth_cv: Coefficient of variation for depth estimates
#'   - depth_mean, depth_sd: Summary statistics
#'   - convergence_rate: Proportion of successful bootstrap iterations
#'   - n_valid: Number of valid bootstrap estimates
#'   - valid_depths, valid_features: Raw bootstrap distributions
#'
#' @details
#' The bootstrap procedure:
#' 
#' 1. **Resampling**: Generate n_boot bootstrap samples with replacement
#' 2. **Model Fitting**: Fit the specified model to each bootstrap sample
#' 3. **Threshold Calculation**: Compute threshold depth for each successful fit
#' 4. **Quality Control**: Filter out non-convergent or invalid predictions
#' 5. **Confidence Intervals**: Calculate percentile-based confidence intervals
#'
#' **Interpretation Guidelines:**
#' - **Depth CV < 0.05**: High precision, reliable extrapolation
#' - **Depth CV 0.05-0.15**: Moderate precision, reasonable confidence
#' - **Depth CV 0.15-0.30**: Low precision, interpret with caution
#' - **Depth CV > 0.30**: Very low precision, unreliable extrapolation
#'
#' The function provides progress updates and convergence diagnostics to
#' help assess the quality of uncertainty estimates.
#'
#' @examples
#' # Bootstrap uncertainty for 1 f/M threshold
#' curve_data <- data.frame(
#'   total_reads_unified = c(1000, 5000, 10000, 20000),
#'   featureNum = c(2000, 4500, 6000, 7200),
#'   curve_id = "test_curve"
#' )
#' uncertainty <- bootstrap_threshold_uncertainty(
#'   curve_data, "michaelis_menten", 1.0, n_boot = 500
#' )
#' cat("Depth 95% CI:", uncertainty$depth_ci)
#' cat("Precision (CV):", uncertainty$depth_cv)
bootstrap_threshold_uncertainty <- function(curve_data, model_type, threshold, n_boot = 1000) {
  # Prepare original data for bootstrap resampling
  original_data <- curve_data
  x_orig <- original_data$total_reads_unified
  y_orig <- original_data$featureNum
  n_points <- length(x_orig)
  
  cat("    Running", n_boot, "bootstrap iterations for", threshold, "f/M threshold...\n")
  
  # Initialize storage for bootstrap results
  bootstrap_depths <- numeric(n_boot)
  bootstrap_features <- numeric(n_boot)
  convergence_count <- 0
  
  # Progress tracking points
  progress_points <- seq(100, n_boot, by = 100)
  
  # Bootstrap resampling loop
  for (i in 1:n_boot) {
    # Generate bootstrap sample with replacement
    boot_indices <- sample(n_points, n_points, replace = TRUE)
    x_boot <- x_orig[boot_indices]
    y_boot <- y_orig[boot_indices]
    boot_data <- data.frame(
      total_reads_unified = x_boot,
      featureNum = y_boot,
      curve_id = original_data$curve_id[1]
    )
    
    tryCatch({
      # Fit model to bootstrap sample
      boot_fit <- fit_model_robust(x_boot, y_boot, model_type)
      
      if (boot_fit$convergence) {
        convergence_count <- convergence_count + 1
        
        # Calculate threshold depth using model-specific method
        if (model_type == "hill") {
          depth <- calculate_hill_threshold(boot_fit, coef(boot_fit$fit_object), threshold)
        } else {
          depth <- calculate_threshold_depth(model_type, coef(boot_fit$fit_object), threshold)
        }
        
        # Validate threshold calculation
        if (!is.na(depth) && depth > 0 && is.finite(depth)) {
          features <- predict(boot_fit$fit_object, newdata = data.frame(x = depth))
          features <- as.numeric(features)
          
          # Store valid bootstrap results
          if (is.finite(features) && features > 0) {
            bootstrap_depths[i] <- depth
            bootstrap_features[i] <- features
          }
        }
      }
    }, error = function(e) {
      # Silent failure - bootstrap sample doesn't converge
      # This is expected for some bootstrap samples
    })
    
    # Progress reporting
    if (i %in% progress_points) {
      conv_rate <- round(100 * convergence_count / i, 1)
      cat("      Completed", i, "iterations (", conv_rate, "% convergence rate)\n")
    }
  }
  
  # Filter out failed iterations
  valid_depths <- bootstrap_depths[bootstrap_depths > 0 & is.finite(bootstrap_depths)]
  valid_features <- bootstrap_features[bootstrap_features > 0 & is.finite(bootstrap_features)]
  
  # Quality control: require minimum number of valid bootstrap samples
  if (length(valid_depths) < 10) {
    cat("    WARNING: Only", length(valid_depths), "successful bootstrap iterations\n")
    cat("             Bootstrap uncertainty estimates may be unreliable\n")
    return(list(
      depth_ci = c(NA, NA),
      features_ci = c(NA, NA),
      depth_cv = NA,
      convergence_rate = convergence_count / n_boot,
      n_valid = length(valid_depths)
    ))
  }
  
  # Calculate confidence intervals and summary statistics
  alpha <- 1 - CONFIDENCE_LEVEL  # Default: 95% confidence intervals
  depth_ci <- quantile(valid_depths, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  features_ci <- quantile(valid_features, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  
  # Calculate precision metrics
  depth_mean <- mean(valid_depths, na.rm = TRUE)
  depth_sd <- sd(valid_depths, na.rm = TRUE)
  depth_cv <- depth_sd / depth_mean  # Coefficient of variation
  
  # Report bootstrap results
  cat("    Bootstrap summary:", 
      "mean depth =", scales::comma(round(depth_mean)), 
      "±", scales::comma(round(depth_sd)), 
      "(CV =", sprintf("%.3f", depth_cv), ")\n")
  cat("    95% Confidence interval: [", 
      scales::comma(round(depth_ci[1])), ", ", 
      scales::comma(round(depth_ci[2])), "]\n")
  
  return(list(
    depth_ci = as.numeric(depth_ci),
    features_ci = as.numeric(features_ci),
    depth_cv = depth_cv,
    depth_mean = depth_mean,
    depth_sd = depth_sd,
    convergence_rate = convergence_count / n_boot,
    n_valid = length(valid_depths),
    valid_depths = valid_depths,
    valid_features = valid_features
  ))
}

#' Calculate Akaike Weights for Model Selection
#'
#' Computes Akaike weights from AIC values to quantify the relative
#' support for each model in a candidate set. These weights are used
#' for model selection and Bayesian Model Averaging.
#'
#' @param aic_values Named numeric vector of AIC values
#' @return Named numeric vector of Akaike weights (sum to 1)
#'
#' @details
#' Akaike weights are calculated as:
#' w_i = exp(-0.5 * Δ_i) / Σ(exp(-0.5 * Δ_j))
#' 
#' where Δ_i = AIC_i - AIC_min is the AIC difference from the best model.
#'
#' **Interpretation:**
#' - Weight close to 1: Strong evidence for this model
#' - Weight ~0.5: Moderate evidence, competing models exist
#' - Weight <0.1: Weak evidence, model unlikely to be best
#'
#' The weights represent the probability that each model is the best
#' approximating model given the data and candidate model set.
#'
#' @examples
#' aic_vals <- c(linear = 150.2, michaelis_menten = 142.1, power_law = 148.7)
#' weights <- calculate_akaike_weights(aic_vals)
#' # michaelis_menten will have highest weight (lowest AIC)
calculate_akaike_weights <- function(aic_values) {
  # Remove any NA values and corresponding names
  valid_idx <- !is.na(aic_values)
  valid_aic <- aic_values[valid_idx]
  valid_names <- names(aic_values)[valid_idx]
  
  if (length(valid_aic) == 0) {
    return(setNames(numeric(0), character(0)))
  }
  
  # Calculate delta AIC relative to best model
  min_aic <- min(valid_aic)
  delta_aic <- valid_aic - min_aic
  
  # Calculate Akaike weights using exponential transformation
  # This gives higher weight to models with lower AIC (better fit)
  relative_likelihood <- exp(-0.5 * delta_aic)
  weights <- relative_likelihood / sum(relative_likelihood)
  
  return(setNames(weights, valid_names))
}

#' Classify Model Uncertainty Based on AIC Competition
#'
#' Categorizes the level of model selection uncertainty based on AIC
#' differences and Akaike weights. This helps interpret the reliability
#' of model-based predictions.
#'
#' @param aic_values Named numeric vector of AIC values
#' @param akaike_weights Named numeric vector of Akaike weights
#'
#' @return List containing:
#'   - classification: Uncertainty category (character)
#'   - best_model_weight: Weight of the best model
#'   - competitive_models: Models with strong support (ΔAIC < 2)
#'   - substantial_support_models: Models with some support (ΔAIC < 4)
#'   - moderate_support_models: Models with weak support (ΔAIC < 7)
#'   - delta_aic: AIC differences from best model
#'
#' @details
#' **AIC Difference Interpretation (Burnham & Anderson, 2002):**
#' - ΔAIC < 2: Substantial support, competitive models
#' - ΔAIC 2-4: Some support, consider for averaging
#' - ΔAIC 4-7: Considerably less support
#' - ΔAIC > 7: Essentially no support
#'
#' **Uncertainty Classifications:**
#' - **Clear Winner**: Single model dominates (weight > 0.9, ΔAIC > 7)
#' - **Low Uncertainty**: Best model clear but some alternatives (weight > 0.7)
#' - **Moderate Uncertainty**: Competing models but clear leader (weight > 0.5)
#' - **High Uncertainty**: Multiple competitive models (weight < 0.5)
#'
#' This classification guides interpretation of model predictions and
#' determines whether Bayesian Model Averaging should be used.
#'
#' @examples
#' aic_vals <- c(linear = 150, michaelis_menten = 142, power_law = 143)
#' weights <- calculate_akaike_weights(aic_vals)
#' uncertainty <- classify_model_uncertainty(aic_vals, weights)
#' cat("Uncertainty level:", uncertainty$classification)
classify_model_uncertainty <- function(aic_values, akaike_weights) {
  # Handle case with no valid models
  valid_idx <- !is.na(aic_values)
  valid_aic <- aic_values[valid_idx]
  valid_weights <- akaike_weights[valid_idx]
  
  if (length(valid_aic) == 0) {
    return(list(
      classification = "No Valid Models",
      best_model_weight = 0,
      competitive_models = character(0),
      substantial_support_models = character(0),
      moderate_support_models = character(0),
      delta_aic = numeric(0)
    ))
  }
  
  # Calculate AIC differences from best model
  min_aic <- min(valid_aic)
  delta_aic <- valid_aic - min_aic
  
  # Categorize models by AIC differences (Burnham & Anderson thresholds)
  competitive_models <- names(delta_aic[delta_aic < 2])      # Strong support
  substantial_support <- names(delta_aic[delta_aic < 4])     # Some support
  moderate_support <- names(delta_aic[delta_aic < 7])        # Weak support
  
  best_model_weight <- max(valid_weights)
  
  # Classification logic based on model competition
  if (length(competitive_models) > 1) {
    # Multiple competitive models exist
    if (best_model_weight < 0.5) {
      classification <- "High Uncertainty"
    } else if (best_model_weight < 0.9) {
      classification <- "Moderate Uncertainty"
    } else {
      classification <- "Low Uncertainty"
    }
  } else if (length(substantial_support) > 1) {
    # Some alternative models have support
    classification <- "Low Uncertainty"
  } else {
    # Single model dominates
    classification <- "Clear Winner"
  }
  
  return(list(
    classification = classification,
    best_model_weight = best_model_weight,
    competitive_models = competitive_models,
    substantial_support_models = substantial_support,
    moderate_support_models = moderate_support,
    delta_aic = delta_aic
  ))
}

#' Check for Linear Model Alerts
#'
#' Determines if the linear model being competitive indicates insufficient
#' data for capturing rarefaction curve saturation. This serves as a
#' quality control check for model fitting.
#'
#' @param aic_values Named numeric vector of AIC values
#' @param model_uncertainty_info List from classify_model_uncertainty()
#'
#' @return List containing:
#'   - linear_alert: Logical indicating if linear model is problematically competitive
#'   - linear_competitive: Logical indicating if linear model has substantial support
#'   - interpretation: Character string explaining the implications
#'
#' @details
#' **Linear Model Interpretation:**
#' - **Competitive linear model**: May indicate insufficient sequencing depth
#'   to observe saturation, or genuinely linear discovery phase
#' - **Non-competitive linear model**: Saturation pattern is well-captured
#'   by nonlinear models, appropriate for threshold extrapolation
#'
#' This check helps distinguish between datasets that have reached informative
#' saturation versus those still in linear discovery phases.
check_linear_model_alert <- function(aic_values, model_uncertainty_info) {
  # Check if linear model exists and is competitive
  if (!"linear" %in% names(aic_values) || is.na(aic_values["linear"])) {
    return(list(
      linear_alert = FALSE,
      linear_competitive = FALSE,
      interpretation = "Linear model not fitted or failed to converge"
    ))
  }
  
  # Determine if linear model is competitive
  linear_competitive <- "linear" %in% model_uncertainty_info$substantial_support_models
  
  # Generate alert if linear model is competitive
  if (linear_competitive) {
    return(list(
      linear_alert = TRUE,
      linear_competitive = TRUE,
      interpretation = "Linear model competitive - consider deeper sequencing to observe saturation"
    ))
  } else {
    return(list(
      linear_alert = FALSE,
      linear_competitive = FALSE,
      interpretation = "Nonlinear saturation pattern detected - suitable for threshold extrapolation"
    ))
  }
}

#' Create Comprehensive Model Competition Report
#'
#' Generates a detailed report of model competition, AIC values, and
#' Akaike weights for a given set of fitted models. This is used to
#' determine the best model and assess overall model uncertainty.
#'
#' @param all_fits List of fitted model results (output of fit_model_robust)
#' @param best_model_name Character string of the name of the best model
#'
#' @return List containing:
#'   - competition_table: Data frame of model competition details
#'   - uncertainty_classification: Overall uncertainty level
#'   - best_model_weight: Weight of the best model
#'   - competitive_models: Models with strong support (ΔAIC < 2)
#'   - substantial_support_models: Models with some support (ΔAIC < 4)
#'   - linear_alert: Information about the linear model alert
#'   - summary_stats: Basic statistics about the model set
#'
#' @details
#' The function calculates AIC values, Akaike weights, and interprets
#' the model competition based on AIC differences. It also checks for
#' a linear model alert if the linear model is competitive.
#'
#' @examples
#' # Example usage (assuming MODEL_TYPES and CONFIDENCE_LEVEL are defined)
#' # This function is called internally by fit_best_model_aic
#' # For demonstration, let's create dummy data
#' dummy_fits <- list(
#'   linear = list(aic = 150, convergence = TRUE),
#'   michaelis_menten = list(aic = 142, convergence = TRUE),
#'   power_law = list(aic = 148, convergence = TRUE)
#' )
#' names(dummy_fits) <- c("linear", "michaelis_menten", "power_law")
#' competition_report <- create_model_competition_report(dummy_fits, "michaelis_menten")
#' print(competition_report$competition_table)
create_model_competition_report <- function(all_fits, best_model_name) {
  # Extract AIC values from all fitted models
  aic_values <- sapply(all_fits, function(m) {
    if (m$convergence && !is.null(m$aic)) {
      return(m$aic)
    } else {
      return(NA)
    }
  })
  
  # Calculate Akaike weights
  akaike_weights <- calculate_akaike_weights(aic_values)
  
  # Classify uncertainty
  uncertainty_info <- classify_model_uncertainty(aic_values, akaike_weights)
  
  # Check for linear model alert
  linear_alert_info <- check_linear_model_alert(aic_values, uncertainty_info)
  
  # Create competition table
  valid_models <- names(aic_values)[!is.na(aic_values)]
  
  if (length(valid_models) > 0) {
    competition_table <- data.frame(
      Model = valid_models,
      AIC = aic_values[valid_models],
      Delta_AIC = uncertainty_info$delta_aic[valid_models],
      Akaike_Weight = akaike_weights[valid_models],
      Evidence_Ratio = akaike_weights[best_model_name] / akaike_weights[valid_models],
      stringsAsFactors = FALSE
    )
    
    # Sort by AIC (best first)
    competition_table <- competition_table[order(competition_table$AIC), ]
    
    # Add interpretation column
    competition_table$Support_Level <- sapply(competition_table$Delta_AIC, function(delta) {
      if (delta < 2) return("Strong")
      if (delta < 4) return("Substantial") 
      if (delta < 7) return("Moderate")
      if (delta < 10) return("Weak")
      return("Minimal")
    })
    
  } else {
    competition_table <- data.frame()
  }
  
  return(list(
    competition_table = competition_table,
    uncertainty_classification = uncertainty_info$classification,
    best_model_weight = uncertainty_info$best_model_weight,
    competitive_models = uncertainty_info$competitive_models,
    substantial_support_models = uncertainty_info$substantial_support_models,
    linear_alert = linear_alert_info,
    summary_stats = list(
      n_converged_models = length(valid_models),
      n_competitive_models = length(uncertainty_info$competitive_models),
      best_model_evidence_ratio = max(competition_table$Evidence_Ratio, na.rm = TRUE)
    )
  ))
}

#' Fit All Models and Select Best by AIC with Comprehensive Competition Analysis
#'
#' Iterates through all supported model types and fits them to a given
#' rarefaction curve data subset. It then selects the best model based
#' on AIC and generates a comprehensive report of model competition.
#'
#' @param data_subset Data frame containing rarefaction curve data for a single curve
#'
#' @return List containing:
#'   - model: Name of the best model
#'   - curve_name: Identifier for the curve
#'   - convergence: Logical indicating if best model converged
#'   - aic: Akaike Information Criterion of the best model
#'   - r_squared: R-squared of the best model
#'   - parameters: Parameter estimates of the best model
#'   - predictions: Fitted values of the best model
#'   - residuals: Residuals of the best model
#'   - x_values, y_values: Original data
#'   - model_competition_analysis: Detailed report of model competition
#'
#' @details
#' The function:
#' 
#' 1. **Model Fitting**: Iterates through MODEL_TYPES and fits each model
#'    to the provided data subset using fit_model_robust.
#' 2. **Convergence Check**: Filters for converged models and selects the
#'    one with the lowest AIC.
#' 3. **Fallback**: If no models converge, it falls back to a simple linear
#'    model and reports this as a "fallback".
#' 4. **Competition Analysis**: Generates a detailed report of model
#'    competition using create_model_competition_report.
#' 5. **Legacy Compatibility**: Updates model_uncertainty and aic_differences
#'    fields for backward compatibility.
#'
#' @examples
#' # Example usage (assuming MODEL_TYPES and CONFIDENCE_LEVEL are defined)
#' # This function is called internally by create_ensemble_predictions
#' # For demonstration, let's create dummy data
#' dummy_data <- data.frame(
#'   total_reads_unified = c(1000, 5000, 10000, 20000, 50000),
#'   featureNum = c(2000, 4500, 6000, 7200, 8500),
#'   curve_id = "test_curve"
#' )
#' best_model_result <- fit_best_model_aic(dummy_data)
#' print(best_model_result$model_competition_analysis$competition_table)
fit_best_model_aic <- function(data_subset) {
  curve_name <- unique(data_subset$curve_id)[1]
  x <- data_subset$total_reads_unified
  y <- data_subset$featureNum
  
  cat("Fitting models for:", curve_name, "\n")
  
  if (length(x) == 0 || length(y) == 0) {
    return(list(
      model = "no_data",
      curve_name = curve_name,
      convergence = FALSE,
      aic = NA,
      error = "No data points"
    ))
  }
  
  # Fit all models
  all_fits <- lapply(MODEL_TYPES, function(model_type) {
    fit_model_robust(x, y, model_type, curve_name)
  })
  names(all_fits) <- MODEL_TYPES
  
  # Select best model
  converged_models <- all_fits[sapply(all_fits, function(m) m$convergence)]
  
  if (length(converged_models) == 0) {
    # Fallback to simple linear model
    lm_fit <- lm(y ~ x)
    best_fit <- list(
      model = "linear_fallback",
      curve_name = curve_name,
      fit_object = lm_fit,
      r_squared = summary(lm_fit)$r.squared,
      aic = AIC(lm_fit),
      convergence = TRUE,
      predictions = predict(lm_fit),
      residuals = residuals(lm_fit),
      x_values = x,
      y_values = y,
      model_competition_analysis = list(
        uncertainty_classification = "Fallback Only",
        linear_alert = list(
          linear_alert = TRUE,
          alert_message = "All advanced models failed - using linear fallback"
        )
      )
    )
  } else {
    # Find model with lowest AIC
    aic_values <- sapply(converged_models, function(m) m$aic)
    best_model_name <- names(which.min(aic_values))
    best_fit <- converged_models[[best_model_name]]
    
    # Comprehensive model competition analysis
    competition_analysis <- create_model_competition_report(all_fits, best_model_name)
    
    # Add competition analysis to best fit
    best_fit$model_competition_analysis <- competition_analysis
    
    # Legacy compatibility - update model_uncertainty field
    best_fit$model_uncertainty <- competition_analysis$uncertainty_classification
    best_fit$aic_differences <- competition_analysis$competition_table$Delta_AIC
    names(best_fit$aic_differences) <- competition_analysis$competition_table$Model
    
    # Report results
    cat("  MODEL COMPETITION ANALYSIS:\n")
    cat("  Best model:", best_model_name, "(AIC =", round(best_fit$aic, 2), ")\n")
    cat("  Uncertainty:", competition_analysis$uncertainty_classification, "\n")
    cat("  Best model weight:", round(competition_analysis$best_model_weight, 3), "\n")
    
    if (length(competition_analysis$competitive_models) > 1) {
      cat("  Competitive models (ΔAIC < 2):", paste(competition_analysis$competitive_models, collapse = ", "), "\n")
    }
    
    # Linear model alert
    linear_info <- competition_analysis$linear_alert
    if (linear_info$linear_alert) {
      cat("  ⚠️  LINEAR MODEL ALERT:", linear_info$interpretation, "\n")
    } else {
      cat("  ✓ Linear model check:", linear_info$interpretation, "\n")
    }
    
    # Print competition table for models with substantial support
    substantial_models <- competition_analysis$competition_table[
      competition_analysis$competition_table$Support_Level %in% c("Strong", "Substantial"), 
    ]
    if (nrow(substantial_models) > 0) {
      cat("  Models with substantial support:\n")
      for (i in 1:nrow(substantial_models)) {
        cat("    ", substantial_models$Model[i], 
            ": ΔAIC =", round(substantial_models$Delta_AIC[i], 2),
            ", weight =", round(substantial_models$Akaike_Weight[i], 3),
            ", support =", substantial_models$Support_Level[i], "\n")
      }
    }
  }
  
  best_fit$all_models <- all_fits
  return(best_fit)
}

#' Determine if Curve is Eligible for Bayesian Model Averaging
#'
#' Checks if a given rarefaction curve is suitable for Bayesian Model
#' Averaging (BMA) based on model competition results.
#'
#' @param competition_analysis List from create_model_competition_report()
#'
#' @return List containing:
#'   - use_bma: Logical indicating if BMA is applicable
#'   - eligible_models: Names of competitive models for BMA
#'   - eligible_weights: Akaike weights of competitive models
#'   - eligible_aic: AIC values of competitive models
#'   - eligible_delta_aic: ΔAIC values of competitive models
#'   - reason: Explanation for BMA eligibility
#'
#' @details
#' BMA is only applicable if:
#' 1. The linear model is not dominant (ΔAIC > 2)
#' 2. There are multiple competitive models (weight < 0.5)
#' 3. The best model has substantial support (weight > 0.7)
#'
#' If these conditions are met, it returns the names of eligible models
#' and their weights for BMA.
#'
#' @examples
#' # Example usage (assuming MODEL_TYPES and CONFIDENCE_LEVEL are defined)
#' # This function is called internally by create_ensemble_predictions
#' # For demonstration, let's create dummy data
#' dummy_competition_analysis <- list(
#'   competition_table = data.frame(
#'     Model = c("linear", "michaelis_menten", "power_law"),
#'     AIC = c(150, 142, 148),
#'     Delta_AIC = c(0, 8, 2),
#'     Akaike_Weight = c(0.5, 0.3, 0.2),
#'     stringsAsFactors = FALSE
#'   ),
#'   uncertainty_classification = "Moderate Uncertainty",
#'   best_model_weight = 0.5,
#'   competitive_models = c("michaelis_menten", "power_law"),
#'   substantial_support_models = c("michaelis_menten"),
#'   linear_alert = list(linear_alert = FALSE, linear_competitive = FALSE, interpretation = "Nonlinear saturation pattern detected")
#' )
#' eligibility <- determine_ensemble_eligibility(dummy_competition_analysis)
#' cat("BMA Eligibility:", eligibility$use_bma)
#' cat("Eligible models:", eligibility$eligible_models)
determine_ensemble_eligibility <- function(competition_analysis) {
  # Extract model competition information
  competition_table <- competition_analysis$competition_table
  linear_alert <- competition_analysis$linear_alert
  
  # Exclude linear model from ensemble consideration
  saturation_models <- competition_table[competition_table$Model != "linear", ]
  
  # Check if linear model is dominant (insufficient curvature)
  if (linear_alert$linear_competitive) {
    return(list(
      use_bma = FALSE,
      eligible_models = character(0),
      reason = "Linear model dominant - insufficient curvature"
    ))
  }
  
  # Find models with substantial support (ΔAIC < threshold)
  eligible_models <- saturation_models[saturation_models$Delta_AIC < BMA_INCLUSION_THRESHOLD, ]
  
  # Check minimum model requirement
  if (nrow(eligible_models) < BMA_MIN_MODELS) {
    return(list(
      use_bma = FALSE,
      eligible_models = eligible_models$Model,
      reason = paste("Only", nrow(eligible_models), "eligible models (need", BMA_MIN_MODELS, "minimum)")
    ))
  }
  
  return(list(
    use_bma = TRUE,
    eligible_models = eligible_models$Model,
    eligible_weights = eligible_models$Akaike_Weight,
    eligible_aic = eligible_models$AIC,
    eligible_delta_aic = eligible_models$Delta_AIC,
    reason = paste("Ensemble eligible:", nrow(eligible_models), "competitive models")
  ))
}

#' Bootstrap Ensemble Uncertainty with Full Statistical Rigor
#'
#' Performs bootstrap resampling to quantify uncertainty in Bayesian
#' Model Averaging (BMA) predictions. This is crucial for understanding
#' the reliability of ensemble-based sequencing depth recommendations.
#'
#' @param curve_data Data frame containing rarefaction curve data
#' @param all_fitted_models List of all fitted model results (output of fit_model_robust)
#' @param competitive_model_names Character vector of names of competitive models
#' @param threshold Numeric value of the marginal returns threshold (features/million reads)
#' @param n_boot Number of bootstrap iterations (default: BMA_BOOTSTRAP_SAMPLES)
#'
#' @return List containing:
#'   - ensemble_ci: 95% confidence interval for ensemble prediction
#'   - ensemble_mean: Mean of ensemble predictions
#'   - ensemble_sd: Standard deviation of ensemble predictions
#'   - model_uncertainty: Uncertainty due to model parameter variation
#'   - total_uncertainty: Total uncertainty (bootstrap SD)
#'   - convergence_rate: Proportion of successful bootstrap iterations
#'   - n_valid: Number of valid bootstrap estimates
#'   - weight_means: Mean Akaike weights across bootstrap samples
#'   - weight_uncertainty: Coefficient of variation of Akaike weights
#'   - robustness_score: 1 - (coefficient of variation of ensemble predictions)
#'   - valid_predictions: Raw bootstrap ensemble predictions
#'   - valid_weights: Raw bootstrap Akaike weights
#'
#' @details
#' The bootstrap procedure:
#' 
#' 1. **Resampling**: Generate n_boot bootstrap samples with replacement
#' 2. **Model Fitting**: For each bootstrap sample, fit all competitive models
#' 3. **Threshold Calculation**: Compute threshold depth for each successful fit
#' 4. **Quality Control**: Filter out non-convergent or invalid predictions
#' 5. **Weight Calculation**: Calculate Akaike weights for each bootstrap sample
#' 6. **Ensemble Prediction**: Calculate weighted ensemble prediction
#' 7. **Confidence Intervals**: Calculate percentile-based confidence intervals
#'
#' The function provides progress updates and convergence diagnostics to
#' help assess the quality of uncertainty estimates.
#'
#' @examples
#' # Example usage (assuming MODEL_TYPES and CONFIDENCE_LEVEL are defined)
#' # This function is called internally by create_ensemble_predictions
#' # For demonstration, let's create dummy data
#' dummy_data <- data.frame(
#'   total_reads_unified = c(1000, 5000, 10000, 20000, 50000),
#'   featureNum = c(2000, 4500, 6000, 7200, 8500),
#'   curve_id = "test_curve"
#' )
#' dummy_fits <- list(
#'   linear = list(aic = 150, convergence = TRUE),
#'   michaelis_menten = list(aic = 142, convergence = TRUE),
#'   power_law = list(aic = 148, convergence = TRUE)
#' )
#' names(dummy_fits) <- c("linear", "michaelis_menten", "power_law")
#' ensemble_uncertainty <- bootstrap_ensemble_uncertainty(
#'   dummy_data, dummy_fits, c("michaelis_menten", "power_law"), 1.0, n_boot = 500
#' )
#' cat("Ensemble 95% CI:", ensemble_uncertainty$ensemble_ci)
#' cat("Model Weight CV:", ensemble_uncertainty$weight_uncertainty)
bootstrap_ensemble_uncertainty <- function(curve_data, all_fitted_models, competitive_model_names, threshold, n_boot = BMA_BOOTSTRAP_SAMPLES) {
  # Validate required constants exist
  if (!exists("CONFIDENCE_LEVEL") || is.null(CONFIDENCE_LEVEL)) {
    stop("CONFIDENCE_LEVEL not defined")
  }
  if (!exists("BMA_MIN_MODELS") || is.null(BMA_MIN_MODELS)) {
    stop("BMA_MIN_MODELS not defined")
  }
  
  cat("    DEBUG: Starting bootstrap_ensemble_uncertainty\n")
  cat("      curve_data dims:", dim(curve_data), "\n")
  cat("      competitive_model_names:", competitive_model_names, "\n")
  cat("      threshold:", threshold, "\n")
  cat("      n_boot:", n_boot, "\n")
  
  original_data <- curve_data
  x_orig <- original_data$total_reads_unified
  y_orig <- original_data$featureNum
  n_points <- length(x_orig)
  curve_name <- unique(original_data$curve_id)[1]
  
  cat("    Running ensemble bootstrap (", n_boot, "iterations) for", threshold, "f/M...\n")
  
  # Storage for bootstrap results
  ensemble_predictions <- numeric(n_boot)
  model_weights_matrix <- matrix(NA, nrow = n_boot, ncol = length(competitive_model_names))
  colnames(model_weights_matrix) <- competitive_model_names
  valid_iterations <- 0
  
  # Progress tracking
  progress_points <- seq(100, n_boot, by = 100)
  
  for (i in 1:n_boot) {
    tryCatch({
      # Bootstrap sample
      boot_indices <- sample(n_points, n_points, replace = TRUE)
      x_boot <- x_orig[boot_indices]
      y_boot <- y_orig[boot_indices]
      
      # Fit all competitive models to bootstrap sample
      boot_fits <- list()
      boot_converged <- character()
      
      for (model_type in competitive_model_names) {
        boot_fit <- fit_model_robust(x_boot, y_boot, model_type, curve_name)
        if (boot_fit$convergence) {
          boot_fits[[model_type]] <- boot_fit
          boot_converged <- c(boot_converged, model_type)
        }
      }
      
      # Require at least 2 converged models for meaningful ensemble
      if (length(boot_converged) >= BMA_MIN_MODELS) {
        # Calculate AIC values for converged models
        boot_aic_values <- sapply(boot_fits[boot_converged], function(m) m$aic)
        names(boot_aic_values) <- boot_converged
        
        # Calculate Akaike weights for this bootstrap sample
        boot_akaike_weights <- calculate_akaike_weights(boot_aic_values)
        
        # Calculate threshold predictions for each model
        boot_predictions <- numeric(length(boot_converged))
        names(boot_predictions) <- boot_converged
        
        for (model_name in boot_converged) {
          model_fit <- boot_fits[[model_name]]
          
          # Calculate threshold depth
          if (model_name == "hill") {
            depth <- calculate_hill_threshold(model_fit, coef(model_fit$fit_object), threshold)
          } else {
            depth <- calculate_threshold_depth(model_name, coef(model_fit$fit_object), threshold)
          }
          
          if (!is.na(depth) && depth > 0 && is.finite(depth)) {
            boot_predictions[model_name] <- depth
          } else {
            # Remove this model from bootstrap iteration if threshold calculation fails
            boot_converged <- setdiff(boot_converged, model_name)
          }
        }
        
        # Final check: still have enough models after threshold calculations?
        if (length(boot_converged) >= BMA_MIN_MODELS) {
          # Renormalize weights for remaining models
          boot_akaike_weights <- boot_akaike_weights[boot_converged]
          boot_akaike_weights <- boot_akaike_weights / sum(boot_akaike_weights)
          boot_predictions <- boot_predictions[boot_converged]
          
          # Calculate weighted ensemble prediction
          ensemble_predictions[i] <- sum(boot_akaike_weights * boot_predictions)
          
          # Debug extreme values
          if (is.infinite(ensemble_predictions[i]) || ensemble_predictions[i] > 1e50) {
            cat("    DEBUG: Extreme value detected at iteration", i, "\n")
            cat("      ensemble_predictions[i]:", ensemble_predictions[i], "\n")
            cat("      boot_akaike_weights:", boot_akaike_weights, "\n")
            cat("      boot_predictions:", boot_predictions, "\n")
            cat("      boot_converged models:", boot_converged, "\n")
          }
          
          # Store weights (fill with 0 for missing models)
          for (model_name in competitive_model_names) {
            if (model_name %in% boot_converged) {
              model_weights_matrix[i, model_name] <- boot_akaike_weights[model_name]
            } else {
              model_weights_matrix[i, model_name] <- 0
            }
          }
          
          valid_iterations <- valid_iterations + 1
        }
      }
      
    }, error = function(e) {
      # Silent failure for this bootstrap iteration
    })
    
    # Progress reporting
    if (i %in% progress_points) {
      conv_rate <- round(100 * valid_iterations / i, 1)
      cat("      Completed", i, "iterations (", conv_rate, "% valid)\n")
    }
  }
  
  # Filter out failed iterations
  cat("    DEBUG: Pre-filtering array info:\n")
  cat("      ensemble_predictions length:", length(ensemble_predictions), "\n")
  cat("      ensemble_predictions range:", range(ensemble_predictions, na.rm = TRUE), "\n")
  cat("      model_weights_matrix dims:", dim(model_weights_matrix), "\n")
  cat("      valid_iterations:", valid_iterations, "\n")
  
  # Check for extreme values before filtering
  extreme_values <- sum(is.infinite(ensemble_predictions) | ensemble_predictions > 1e50, na.rm = TRUE)
  if (extreme_values > 0) {
    cat("      WARNING: Found", extreme_values, "extreme/infinite values in ensemble_predictions\n")
  }
  
  valid_predictions <- ensemble_predictions[ensemble_predictions > 0 & is.finite(ensemble_predictions)]
  valid_weights <- model_weights_matrix[complete.cases(model_weights_matrix), , drop = FALSE]
  
  cat("    DEBUG: Post-filtering array info:\n")
  cat("      valid_predictions length:", length(valid_predictions), "\n")
  if (length(valid_predictions) > 0) {
    cat("      valid_predictions range:", range(valid_predictions, na.rm = TRUE), "\n")
  }
  cat("      valid_weights dims:", dim(valid_weights), "\n")
  
  if (length(valid_predictions) < 10) {
    cat("    WARNING: Only", length(valid_predictions), "successful ensemble bootstrap iterations\n")
    return(list(
      ensemble_ci = c(NA, NA),
      ensemble_mean = NA,
      ensemble_sd = NA,
      model_uncertainty = NA,
      total_uncertainty = NA,
      convergence_rate = valid_iterations / n_boot,
      n_valid = length(valid_predictions),
      weight_uncertainty = rep(NA, length(competitive_model_names)),
      robustness_score = NA
    ))
  }
  
  # Calculate ensemble statistics
  alpha <- 1 - CONFIDENCE_LEVEL
  cat("    DEBUG: About to calculate ensemble statistics\n")
  cat("      alpha:", alpha, "\n")
  cat("      CONFIDENCE_LEVEL:", CONFIDENCE_LEVEL, "\n")
  
  tryCatch({
    ensemble_ci <- quantile(valid_predictions, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    cat("      ensemble_ci calculated successfully\n")
  }, error = function(e) {
    cat("      ERROR in quantile calculation:", e$message, "\n")
    stop("Quantile calculation failed: ", e$message)
  })
  
  ensemble_mean <- mean(valid_predictions, na.rm = TRUE)
  ensemble_sd <- sd(valid_predictions, na.rm = TRUE)
  
  cat("    DEBUG: About to calculate model weight statistics\n")
  cat("      valid_weights dims:", dim(valid_weights), "\n")
  cat("      valid_weights colnames:", colnames(valid_weights), "\n")
  
  # Calculate model weight statistics
  tryCatch({
    weight_means <- apply(valid_weights, 2, mean, na.rm = TRUE)
    cat("      weight_means calculated successfully\n")
  }, error = function(e) {
    cat("      ERROR in weight_means calculation:", e$message, "\n")
    stop("Weight means calculation failed: ", e$message)
  })
  
  tryCatch({
    weight_sds <- apply(valid_weights, 2, sd, na.rm = TRUE)
    cat("      weight_sds calculated successfully\n")
  }, error = function(e) {
    cat("      ERROR in weight_sds calculation:", e$message, "\n")
    stop("Weight sds calculation failed: ", e$message)
  })
  
  weight_cvs <- weight_sds / weight_means
  weight_cvs[is.na(weight_cvs) | is.infinite(weight_cvs)] <- 0
  
  # Calculate model uncertainty (between-model variance component)
  # This captures how much predictions vary due to model uncertainty
  cat("    DEBUG: Starting model uncertainty calculation...\n")
  cat("      competitive_model_names:", competitive_model_names, "\n")
  cat("      all_fitted_models names:", names(all_fitted_models), "\n")
  
  tryCatch({
    original_predictions <- sapply(competitive_model_names, function(model_name) {
      cat("      Processing model:", model_name, "\n")
      
      if (model_name %in% names(all_fitted_models) && all_fitted_models[[model_name]]$convergence) {
        cat("        Model exists and converged\n")
        model_fit <- all_fitted_models[[model_name]]
        
        # Safely extract parameters and calculate threshold
        tryCatch({
          # DEBUG: Check model fit structure
          cat("        fit_object class:", class(model_fit$fit_object), "\n")
          
          # DEBUG: Check parameter extraction
          coef_values <- coef(model_fit$fit_object)
          cat("        coef_values:", coef_values, "\n")
          cat("        coef_values length:", length(coef_values), "\n")
          cat("        coef_values names:", names(coef_values), "\n")
          
          # DEBUG: Before threshold calculation
          cat("        About to calculate threshold for", model_name, "\n")
          
          if (model_name == "hill") {
            depth <- calculate_hill_threshold(model_fit, coef_values, threshold)
          } else {
            depth <- calculate_threshold_depth(model_name, coef_values, threshold)
          }
          
          cat("        Threshold calculated:", depth, "\n")
          return(ifelse(is.na(depth) || depth <= 0 || !is.finite(depth), NA, depth))
          
        }, error = function(e) {
          cat("        ERROR: Failed to process model", model_name, "\n")
          cat("        Error message:", e$message, "\n")
          return(NA)
        })
      } else {
        cat("        Model missing or failed to converge\n")
        return(NA)
      }
    })
    
    cat("      ✓ original_predictions calculated successfully\n")
    cat("      original_predictions:", original_predictions, "\n")
    
  }, error = function(e) {
    cat("      ERROR in model uncertainty calculation:", e$message, "\n")
    cat("      Setting original_predictions to all NAs\n")
    original_predictions <<- rep(NA, length(competitive_model_names))
    names(original_predictions) <<- competitive_model_names
  })
  
  cat("    DEBUG: Model uncertainty calculation completed\n")
  
  # Complete the model uncertainty calculation
  original_weights <- weight_means / sum(weight_means)  # Normalize
  valid_original <- !is.na(original_predictions)
  
  if (sum(valid_original) >= 2) {
    weighted_mean_pred <- sum(original_weights[valid_original] * original_predictions[valid_original])
    model_uncertainty <- sqrt(sum(original_weights[valid_original] * (original_predictions[valid_original] - weighted_mean_pred)^2))
  } else {
    model_uncertainty <- NA
  }
  
  # Robustness score: 1 - (coefficient of variation of ensemble predictions)
  robustness_score <- 1 - (ensemble_sd / ensemble_mean)
  robustness_score <- max(0, min(1, robustness_score))  # Bound between 0 and 1
  
  cat("    Ensemble bootstrap results: mean =", round(ensemble_mean), 
      "±", round(ensemble_sd), "(", round(100*ensemble_sd/ensemble_mean, 1), "% CV)\n")
  cat("    Model uncertainty: ±", round(model_uncertainty), "| Robustness:", round(robustness_score, 3), "\n")
  
  return(list(
    ensemble_ci = as.numeric(ensemble_ci),
    ensemble_mean = ensemble_mean,
    ensemble_sd = ensemble_sd,
    model_uncertainty = model_uncertainty,
    total_uncertainty = ensemble_sd,  # Bootstrap SD includes both parameter and model uncertainty
    convergence_rate = valid_iterations / n_boot,
    n_valid = length(valid_predictions),
    weight_means = weight_means,
    weight_uncertainty = weight_cvs,
    robustness_score = robustness_score,
    valid_predictions = valid_predictions,
    valid_weights = valid_weights
  ))
}

#' Create Ensemble Predictions with Full Uncertainty Quantification
#'
#' Combines single model predictions with Bayesian Model Averaging (BMA)
#' to provide a robust and uncertainty-aware prediction for a given threshold.
#'
#' @param fitted_model_result List from fit_best_model_aic()
#' @param thresholds Numeric vector of marginal returns thresholds (features/million reads)
#'
#' @return Data frame with predictions for each threshold, including BMA columns
#'
#' @details
#' The function:
#' 
#' 1. **Single Model Approach**: If BMA is not enabled or not applicable,
#'    it calls find_threshold_depths() to get single model predictions.
#' 2. **BMA Eligibility Check**: Determines if the curve is suitable for BMA
#'    based on model competition results.
#' 3. **BMA Application**: If BMA is applicable, it:
#'    - Fits all competitive models to the original data.
#'    - Runs ensemble bootstrap uncertainty for each threshold.
#'    - Updates the single_model_results data frame with BMA-specific columns.
#' 4. **Legacy Compatibility**: Ensures backward compatibility by
#'    updating model_uncertainty, total_uncertainty, ensemble_models,
#'    ensemble_weights, robustness_score, and bma_reason fields.
#'
#' @examples
#' # Example usage (assuming MODEL_TYPES and CONFIDENCE_LEVEL are defined)
#' # This function is called internally by create_ensemble_predictions
#' # For demonstration, let's create dummy data
#' dummy_data <- data.frame(
#'   total_reads_unified = c(1000, 5000, 10000, 20000, 50000),
#'   featureNum = c(2000, 4500, 6000, 7200, 8500),
#'   curve_id = "test_curve"
#' )
#' dummy_fits <- list(
#'   linear = list(aic = 150, convergence = TRUE),
#'   michaelis_menten = list(aic = 142, convergence = TRUE),
#'   power_law = list(aic = 148, convergence = TRUE)
#' )
#' names(dummy_fits) <- c("linear", "michaelis_menten", "power_law")
#' best_model_result <- fit_best_model_aic(dummy_data)
#' ensemble_predictions <- create_ensemble_predictions(best_model_result, c(1.0, 2.0))
#' print(ensemble_predictions)
create_ensemble_predictions <- function(fitted_model_result, thresholds = MARGINAL_THRESHOLDS) {
  # Check if BMA is enabled
  if (!USE_BAYESIAN_MODEL_AVERAGING) {
    # Fall back to single model approach
    return(find_threshold_depths(fitted_model_result, thresholds))
  }
  
  # Get competition analysis
  if (is.null(fitted_model_result$model_competition_analysis)) {
    cat("  No competition analysis available, using single model approach\n")
    return(find_threshold_depths(fitted_model_result, thresholds))
  }
  
  competition_analysis <- fitted_model_result$model_competition_analysis
  
  # Check ensemble eligibility
  eligibility <- determine_ensemble_eligibility(competition_analysis)
  
  if (!eligibility$use_bma) {
    cat("  BMA not applicable:", eligibility$reason, "\n")
    cat("  Using single model approach\n")
    single_model_results <- find_threshold_depths(fitted_model_result, thresholds)
    
    # Add BMA status columns
    single_model_results$bma_applied <- FALSE
    single_model_results$bma_depth <- NA
    single_model_results$bma_features <- NA
    single_model_results$model_uncertainty <- NA
    single_model_results$total_uncertainty <- NA
    single_model_results$ensemble_models <- NA
    single_model_results$ensemble_weights <- NA
    single_model_results$robustness_score <- NA
    single_model_results$bma_reason <- eligibility$reason
    
    return(single_model_results)
  }
  
  cat("  BMA eligible:", eligibility$reason, "\n")
  cat("  Competitive models:", paste(eligibility$eligible_models, collapse = ", "), "\n")
  
  # Get single model results first
  single_model_results <- find_threshold_depths(fitted_model_result, thresholds)
  
  if (nrow(single_model_results) == 0) {
    cat("  No threshold calculations possible\n")
    return(single_model_results)
  }
  
  # Add BMA columns
  single_model_results$bma_applied <- TRUE
  single_model_results$bma_depth <- NA
  single_model_results$bma_features <- NA
  single_model_results$model_uncertainty <- NA
  single_model_results$total_uncertainty <- NA
  single_model_results$ensemble_models <- paste(eligibility$eligible_models, collapse = "|")
  single_model_results$ensemble_weights <- paste(round(eligibility$eligible_weights, 3), collapse = "|")
  single_model_results$robustness_score <- NA
  single_model_results$bma_reason <- "Ensemble applied"
  
  # Get curve data for bootstrap (need to reconstruct from fitted model)
  curve_data <- data.frame(
    total_reads_unified = fitted_model_result$x_values,
    featureNum = fitted_model_result$y_values,
    curve_id = fitted_model_result$curve_name
  )
  
  # Run ensemble bootstrap for each threshold
  for (i in 1:nrow(single_model_results)) {
    thresh <- single_model_results$threshold[i]
    
    cat("  Running ensemble analysis for", thresh, "f/M threshold...\n")
    
    # Run bootstrap ensemble uncertainty
    tryCatch({
      ensemble_result <- bootstrap_ensemble_uncertainty(
        curve_data = curve_data,
        all_fitted_models = fitted_model_result$all_models,
        competitive_model_names = eligibility$eligible_models,
        threshold = thresh,
        n_boot = BMA_BOOTSTRAP_SAMPLES
      )
    }, error = function(e) {
      cat("ERROR: bootstrap_ensemble_uncertainty failed\n")
      cat("  Curve:", fitted_model_result$curve_name, "\n")
      cat("  Threshold:", thresh, "f/M\n")
      cat("  Competitive models:", paste(eligibility$eligible_models, collapse = ", "), "\n")
      cat("  Error message:", e$message, "\n")
      cat("  Call stack:\n")
      traceback()
      stop("Ensemble bootstrap failed for ", fitted_model_result$curve_name, " at ", thresh, " f/M: ", e$message)
    })
    
    # DEBUG: Check ensemble_result structure
    cat("    DEBUG: ensemble_result structure:\n")
    cat("      class:", class(ensemble_result), "\n")
    cat("      names:", names(ensemble_result), "\n")
    
    # DEBUG: Check each component before accessing
    tryCatch({
      cat("      ensemble_mean exists:", "ensemble_mean" %in% names(ensemble_result), "\n")
      if ("ensemble_mean" %in% names(ensemble_result)) {
        cat("      ensemble_mean value:", ensemble_result$ensemble_mean, "\n")
        cat("      ensemble_mean is.na:", is.na(ensemble_result$ensemble_mean), "\n")
      }
    }, error = function(e) {
      cat("      ERROR accessing ensemble_mean:", e$message, "\n")
    })
    
    # DEBUG: Check array dimensions before assignment
    cat("      single_model_results dims:", dim(single_model_results), "\n")
    cat("      current index i:", i, "\n")
    cat("      nrow(single_model_results):", nrow(single_model_results), "\n")
    
    tryCatch({
      if (!is.na(ensemble_result$ensemble_mean)) {
        cat("    DEBUG: Starting result assignment...\n")
        
        # DEBUG: Check bma_depth column exists and has correct length
        cat("      bma_depth column exists:", "bma_depth" %in% names(single_model_results), "\n")
        cat("      bma_depth length:", length(single_model_results$bma_depth), "\n")
        
        cat("      Assigning ensemble_mean to bma_depth[", i, "]...\n")
        single_model_results$bma_depth[i] <- ensemble_result$ensemble_mean
        cat("      ✓ bma_depth assigned successfully\n")
        
        cat("      Assigning model_uncertainty to index", i, "...\n")
        single_model_results$model_uncertainty[i] <- ensemble_result$model_uncertainty
        cat("      ✓ model_uncertainty assigned successfully\n")
        
        cat("      Assigning total_uncertainty to index", i, "...\n")
        single_model_results$total_uncertainty[i] <- ensemble_result$total_uncertainty
        cat("      ✓ total_uncertainty assigned successfully\n")
        
        cat("      Assigning robustness_score to index", i, "...\n")
        single_model_results$robustness_score[i] <- ensemble_result$robustness_score
        cat("      ✓ robustness_score assigned successfully\n")
        
        # Calculate BMA features at predicted depth
        cat("      About to predict BMA features...\n")
        cat("      fitted_model_result$fit_object class:", class(fitted_model_result$fit_object), "\n")
        cat("      ensemble_mean for prediction:", ensemble_result$ensemble_mean, "\n")
        
        bma_features <- predict(fitted_model_result$fit_object, 
                               newdata = data.frame(x = ensemble_result$ensemble_mean))
        cat("      ✓ predict() completed successfully\n")
        cat("      bma_features value:", bma_features, "\n")
        
        cat("      Assigning bma_features to index", i, "...\n")
        single_model_results$bma_features[i] <- as.numeric(bma_features)
        cat("      ✓ bma_features assigned successfully\n")
        
        # Update ensemble weights with bootstrap-averaged weights
        cat("      Checking weight_means...\n")
        if (!is.null(ensemble_result$weight_means)) {
          cat("      weight_means length:", length(ensemble_result$weight_means), "\n")
          cat("      weight_means values:", ensemble_result$weight_means, "\n")
          
          weight_strings <- paste(round(ensemble_result$weight_means, 3), collapse = "|")
          cat("      weight_strings created:", weight_strings, "\n")
          
          cat("      Assigning ensemble_weights to index", i, "...\n")
          single_model_results$ensemble_weights[i] <- weight_strings
          cat("      ✓ ensemble_weights assigned successfully\n")
        } else {
          cat("      weight_means is NULL, skipping weight assignment\n")
        }
        
        cat("    DEBUG: All assignments completed successfully for index", i, "\n")
      } else {
        cat("    DEBUG: ensemble_mean is NA, skipping assignments\n")
      }
    }, error = function(e) {
      cat("    ERROR during result assignment at index", i, ":\n")
      cat("      Error message:", e$message, "\n")
      cat("      Call stack:\n")
      traceback()
      stop("Result assignment failed at index ", i, ": ", e$message)
    })
  }
  
  cat("  BMA analysis complete\n")
  return(single_model_results)
}