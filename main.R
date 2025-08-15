# Main analysis workflow

# Load required modules for ORACLE-seq analysis
source(file.path(ORACLE_SEQ_DIR, "R", "config.R"))
source(file.path(ORACLE_SEQ_DIR, "R", "utils.R"))
source(file.path(ORACLE_SEQ_DIR, "R", "models.R"))
source(file.path(ORACLE_SEQ_DIR, "R", "analysis.R"))
source(file.path(ORACLE_SEQ_DIR, "R", "plots.R"))

# Main analysis function
run_rarefaction_analysis <- function() {
  # Setup
  setwd(BASE_DIR)
  load_libraries()
  theme_set(theme_minimal())
  set.seed(SEED)
  ensure_directory(PLOTS_DIR)
  
  cat("=== RAREFACTION CURVE ANALYSIS ===\n\n")
  
  # Load data
  cat("Loading data...\n")
  data_list <- load_rarefaction_data()
  rarefaction_data <- data_list$data
  curve_metadata <- data_list$metadata
  
  # Create reads per cell converter
  cat("\nCreating reads per cell converter...\n")
  reads_per_cell_converter <- create_reads_per_cell_converter(rarefaction_data)
  
  # Fit models with bootstrap uncertainty
  cat("\nFitting models with bootstrap...\n")
  fitted_models_with_uncertainty <- rarefaction_data %>%
    group_by(curve_id) %>%
    group_modify(~ {
      current_curve_id <- .y$curve_id
      data_with_id <- .x %>% mutate(curve_id = current_curve_id)
      
      # Fit primary model
      result <- fit_best_model_aic(data_with_id)
      
      # Run bootstrap for all thresholds if model converged
      bootstrap_results <- list()
      if (result$convergence) {
        for (thresh in MARGINAL_THRESHOLDS) {
          cat("\n  Bootstrap analysis for", current_curve_id, "at", thresh, "f/M:\n")
          boot_result <- bootstrap_threshold_uncertainty(
            curve_data = data_with_id,
            model_type = result$model,
            threshold = thresh,
            n_boot = N_BOOTSTRAP
          )
          bootstrap_results[[as.character(thresh)]] <- boot_result
        }
      }
      
      tibble(
        model_type = result$model,
        r_squared = result$r_squared,
        aic = result$aic %||% NA,
        convergence = result$convergence,
        model_uncertainty = result$model_uncertainty %||% "unknown",
        fitted_model = list(result),
        bootstrap_results = list(bootstrap_results)
      )
    })
  
  # Enhanced marginal returns calculation with confidence intervals
  cat("\nCalculating enhanced marginal returns with uncertainty...\n")
  marginal_returns_results <- data.frame()
  
  for (i in 1:nrow(fitted_models_with_uncertainty)) {
    if (!fitted_models_with_uncertainty$fitted_model[[i]]$convergence) next
    
    curve_id <- fitted_models_with_uncertainty$curve_id[i]
    model_result <- fitted_models_with_uncertainty$fitted_model[[i]]
    
    cat("Processing:", curve_id, "\n")
    
    # Get threshold depths with Bayesian Model Averaging if applicable
    thresh_depths <- create_ensemble_predictions(model_result)
    
    if (nrow(thresh_depths) > 0) {
      # Add metadata
      thresh_depths$curve_id <- curve_id
      thresh_depths$technology <- str_extract(curve_id, "^[^_]+")
      thresh_depths$sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
      thresh_depths$feature_type <- str_extract(curve_id, "[^_]+$")
      thresh_depths$model_type <- model_result$model
      thresh_depths$current_max_depth <- max(model_result$x_values)
      thresh_depths$current_max_features <- max(model_result$y_values)
      
      # Add bootstrap statistics if available
      bootstrap_results <- fitted_models_with_uncertainty$bootstrap_results[[i]]
      for (j in 1:nrow(thresh_depths)) {
        thresh <- thresh_depths$threshold[j]
        boot_stats <- bootstrap_results[[as.character(thresh)]]
        
        if (!is.null(boot_stats)) {
          thresh_depths$depth_ci_lower[j] <- boot_stats$depth_ci[1]
          thresh_depths$depth_ci_upper[j] <- boot_stats$depth_ci[2]
          thresh_depths$depth_cv[j] <- boot_stats$depth_cv
          thresh_depths$features_ci_lower[j] <- boot_stats$features_ci[1]
          thresh_depths$features_ci_upper[j] <- boot_stats$features_ci[2]
          thresh_depths$bootstrap_convergence_rate[j] <- boot_stats$convergence_rate
          thresh_depths$extrapolation_factor[j] <- thresh_depths$depth[j] / max(model_result$x_values)
          
          # Categorize model stability
          thresh_depths$model_stability[j] <- case_when(
            thresh_depths$depth_cv[j] < 0.05 ~ "High",
            thresh_depths$depth_cv[j] < 0.15 ~ "Moderate",
            thresh_depths$depth_cv[j] < 0.30 ~ "Low",
            TRUE ~ "Very Low"
          )
        } else {
          # Fill with NAs if bootstrap failed
          thresh_depths$depth_ci_lower[j] <- NA
          thresh_depths$depth_ci_upper[j] <- NA
          thresh_depths$depth_cv[j] <- NA
          thresh_depths$features_ci_lower[j] <- NA
          thresh_depths$features_ci_upper[j] <- NA
          thresh_depths$bootstrap_convergence_rate[j] <- NA
          thresh_depths$extrapolation_factor[j] <- thresh_depths$depth[j] / max(model_result$x_values)
          thresh_depths$model_stability[j] <- "Unknown"
        }
      }
      
      marginal_returns_results <- rbind(marginal_returns_results, thresh_depths)
    }
  }
  
  # Calculate TD75/TD50 metrics
  cat("\nCalculating TD75/TD50 metrics...\n")
  td_metrics_results <- calculate_td_metrics(
    fitted_models_with_uncertainty = fitted_models_with_uncertainty,
    rarefaction_data = rarefaction_data,
    bootstrap_samples = N_BOOTSTRAP
  )

  # Save results with enhanced data
  cat("\nSaving enhanced results...\n")
  saveRDS(fitted_models_with_uncertainty, file.path(STATS_DIR, "fitted_models_with_uncertainty_mvp.rds"))
  write.csv(marginal_returns_results, 
            file.path(STATS_DIR, "enhanced_marginal_returns_mvp.csv"), 
            row.names = FALSE)
  
  # Save TD75/TD50 results
  if (nrow(td_metrics_results) > 0) {
    write.csv(td_metrics_results,
              file.path(STATS_DIR, "td_metrics_results.csv"),
              row.names = FALSE)
    saveRDS(td_metrics_results, file.path(STATS_DIR, "td_metrics_results.rds"))
    cat("TD75/TD50 results saved to:", file.path(STATS_DIR, "td_metrics_results.csv"), "\n")
  }
  
  # Create plots using the enhanced fitted models
  cat("\nCreating comprehensive visualizations...\n")
  
  # Create fitted curves data from enhanced models
  fitted_curves <- create_fitted_curves_data(fitted_models_with_uncertainty, rarefaction_data)
  
  # 1. Reads conversion model plot
  cat("  Creating reads conversion model plot...\n")
  conversion_plot <- create_conversion_model_plot(rarefaction_data, reads_per_cell_converter)
  save_plots(list(conversion_plot), "reads_conversion_model.pdf", width = 10, height = 6)
  
  # 2. Individual model fit plots (two versions)
  cat("  Creating individual model fit plots...\n")
  individual_plots <- create_individual_fit_plots(fitted_models_with_uncertainty, rarefaction_data, fitted_curves)
  save_plots(individual_plots, "individual_model_fits_total_reads.pdf", width = 10, height = 6)
  
  # Create reads per cell version
  individual_plots_rpc <- create_individual_fit_plots_reads_per_cell(fitted_models_with_uncertainty, rarefaction_data, fitted_curves, reads_per_cell_converter)
  save_plots(individual_plots_rpc, "individual_model_fits_reads_per_cell.pdf", width = 10, height = 6)
  
  # 3. Extended saturation curves
  cat("  Creating extended saturation curves...\n")
  extended_plots <- create_extended_saturation_plots(fitted_models_with_uncertainty, rarefaction_data, marginal_returns_results, reads_per_cell_converter)
  save_plots(extended_plots, "extended_saturation_curves_individual.pdf", width = 10, height = 6)
  
  # 4. Combined extended saturation curves
  cat("  Creating combined extended saturation curves...\n")
  combined_extended_plots <- create_combined_extended_saturation_plots(fitted_models_with_uncertainty, rarefaction_data, marginal_returns_results, reads_per_cell_converter)
  save_plots(combined_extended_plots, "extended_saturation_curves_combined.pdf", width = 12, height = 8)
  
  # Create reads per cell version of combined plots
  combined_extended_plots_rpc <- lapply(names(combined_extended_plots), function(plot_name) {
    orig_plot <- combined_extended_plots[[plot_name]]
    # Replace x-axis with reads per cell conversion
    orig_plot + 
      aes(x = reads_per_cell_converter(total_reads_unified)) +
      labs(x = "Estimated Median Reads per Cell")
  })
  names(combined_extended_plots_rpc) <- paste0(names(combined_extended_plots), "_ReadsPerCell")
  save_plots(combined_extended_plots_rpc, "extended_saturation_curves_combined_reads_per_cell.pdf", width = 12, height = 8)
  
  # 5. Features at thresholds bar plots
  cat("  Creating features at thresholds bar plots...\n")
  threshold_plots <- create_threshold_bar_plots(marginal_returns_results)
  save_plots(threshold_plots, "features_at_thresholds_bar_plots.pdf", width = 10, height = 8)
  
  # 6. Reads threshold bar plots (total reads and reads per cell)
  cat("  Creating reads threshold bar plots...\n")
  reads_threshold_plots <- create_reads_threshold_bar_plots(marginal_returns_results, reads_per_cell_converter)
  save_plots(reads_threshold_plots, "reads_at_thresholds_bar_plots.pdf", width = 12, height = 10)
  
  # Create log scale version
  reads_threshold_plots_log <- lapply(names(reads_threshold_plots), function(plot_name) {
    orig_plot <- reads_threshold_plots[[plot_name]]
    orig_plot + scale_x_log10(
      name = paste("Log10", orig_plot$labels$x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  })
  names(reads_threshold_plots_log) <- paste0(names(reads_threshold_plots), "_Log")
  save_plots(reads_threshold_plots_log, "reads_at_thresholds_bar_plots_log_scale.pdf", width = 12, height = 10)
  
  # 7. Marginal returns derivative curves
  cat("  Creating marginal returns derivative curves...\n")
  marginal_returns_plots <- create_marginal_returns_derivative_plots(fitted_models_with_uncertainty, rarefaction_data, marginal_returns_results, reads_per_cell_converter)
  save_plots(marginal_returns_plots, "marginal_returns_derivative_curves.pdf", width = 12, height = 8)
  
  # 8. Bootstrap distribution plots
  cat("  Creating bootstrap distribution plots...\n")
  bootstrap_plots <- create_bootstrap_distribution_plots(fitted_models_with_uncertainty, marginal_returns_results)
  if (length(bootstrap_plots) > 0) {
    save_plots(bootstrap_plots, "bootstrap_uncertainty_distributions.pdf", width = 14, height = 10)
  }
  
  # 8b. TD75/TD50 enhanced rarefaction plots
  if (nrow(td_metrics_results) > 0) {
    cat("  Creating TD75/TD50 enhanced rarefaction plots...\n")
    td_enhanced_plots <- create_td_enhanced_rarefaction_plots(
      fitted_models_with_uncertainty, 
      rarefaction_data, 
      td_metrics_results, 
      reads_per_cell_converter
    )
    if (length(td_enhanced_plots) > 0) {
      save_plots(td_enhanced_plots, "td_enhanced_rarefaction_plots.pdf", width = 12, height = 8)
    }
  }
  
  # 9. Generate model equations report
  cat("  Generating model equations report...\n")
  create_model_equations_report(fitted_models_with_uncertainty, file.path(STATS_DIR, "model_equations_report.txt"))
  
  # 10. Generate summary analysis report
  cat("  Generating summary analysis report...\n")
  create_summary_analysis_report(fitted_models_with_uncertainty, marginal_returns_results, reads_per_cell_converter, file.path(STATS_DIR, "summary_analysis_report.txt"))
  
  # 11. Generate model competition analysis report
  cat("  Generating model competition analysis report...\n")
  create_model_competition_analysis_report(fitted_models_with_uncertainty, file.path(STATS_DIR, "model_competition_analysis_report.txt"))
  
  # 12. Generate TD75/TD50 HTML report
  if (nrow(td_metrics_results) > 0) {
    cat("  Generating TD75/TD50 HTML report...\n")
    ensure_directory(file.path(BASE_DIR, "reports"))
    create_td_metrics_html_report(
      td_metrics_results,
      reads_per_cell_converter, 
      file.path(BASE_DIR, "reports", "td_metrics_analysis.html")
    )
  }
  
  cat("\n=== ENHANCED ANALYSIS COMPLETE ===\n")
  
  return(list(
    fitted_models = fitted_models_with_uncertainty,
    marginal_returns = marginal_returns_results,
    td_metrics = td_metrics_results,
    reads_converter = reads_per_cell_converter
  ))
}

# Helper function to create fitted curves data
create_fitted_curves_data <- function(fitted_models, rarefaction_data) {
  fitted_curve_data <- list()
  
  for (i in 1:nrow(fitted_models)) {
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    
    if (!model_result$convergence) next
    
    orig_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
    x_range <- range(orig_data$total_reads_unified)
    x_smooth <- seq(from = x_range[1], to = x_range[2], length.out = 100)
    
    y_smooth <- predict(model_result$fit_object, newdata = data.frame(x = x_smooth))
    
    curve_data <- tibble(
      curve_id = curve_id,
      total_reads_unified = x_smooth,
      featureNum_fitted = y_smooth
    )
    
    fitted_curve_data[[i]] <- curve_data
  }
  
  return(bind_rows(fitted_curve_data))
}

# Execute the main analysis
cat("Starting ORACLE-seq rarefaction analysis...\n")
results <- run_rarefaction_analysis()
cat("Analysis complete!\n")