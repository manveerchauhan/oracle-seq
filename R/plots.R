# Plotting functions (simplified for MVP)

# Create reads conversion model plot
create_conversion_model_plot <- function(rarefaction_data, reads_per_cell_converter) {
  # Extract the data used for conversion model
  conversion_data <- rarefaction_data %>%
    filter(!is.na(median_reads_per_cell), median_reads_per_cell > 0) %>%
    mutate(
      technology = str_extract(curve_id, "^[^_]+"),
      sample_type = str_extract(curve_id, "(?<=_)[^_]+(?=_)"),
      total_reads_millions = total_reads_unified / 1e6
    ) %>%
    select(curve_id, technology, sample_type, total_reads_unified, total_reads_millions, median_reads_per_cell)
  
  # Get model statistics
  model_r_squared <- attr(reads_per_cell_converter, "r_squared")
  model_intercept <- attr(reads_per_cell_converter, "intercept")
  model_slope <- attr(reads_per_cell_converter, "slope")
  
  # Create fitted line data
  x_range <- range(conversion_data$total_reads_unified)
  x_seq <- seq(x_range[1], x_range[2], length.out = 100)
  fitted_line <- data.frame(
    total_reads_unified = x_seq,
    total_reads_millions = x_seq / 1e6,
    fitted_reads_per_cell = reads_per_cell_converter(x_seq)
  )
  
  # Format equation text
  equation_text <- paste0(
    "y = ", round(model_intercept, 3), 
    " + ", format(model_slope, scientific = TRUE, digits = 3), 
    " × x"
  )
  
  # Get plot colors
  plot_colors <- get_color_map("conversion_plot")
  
  # Calculate n per technology
  tech_summary <- conversion_data %>%
    group_by(technology) %>%
    summarise(n_points = n() / 2, .groups = 'drop')  # Divide by 2 for pseudoreplication
  
  # Create the plot
  p <- ggplot() +
    # Add fitted line
    geom_line(data = fitted_line, 
              aes(x = total_reads_millions, y = fitted_reads_per_cell),
              color = "black", linewidth = 1.2, alpha = 0.8) +
    # Add data points colored by technology
    geom_point(data = conversion_data,
               aes(x = total_reads_millions, y = median_reads_per_cell, 
                   color = technology, shape = sample_type),
               size = 3, alpha = 0.6) +
    # Facet by technology for clear separation
    facet_wrap(~technology, scales = "free") +
    # Color and shape mappings
    scale_color_manual(values = plot_colors, name = "Technology") +
    scale_shape_manual(values = c("SC" = 16, "SN" = 17), 
                      name = "Sample Type",
                      labels = c("SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
    # Add technology-specific model statistics annotation
    geom_label(data = tech_summary,
               aes(x = Inf, y = Inf, label = paste0(
                 equation_text, "\n",
                 "R² = ", round(model_r_squared, 3), "\n",
                 "n = ", n_points, " rarefaction data points"
               )),
               hjust = 1, vjust = 1, size = 4, 
               fill = "white", alpha = 0.8, label.padding = unit(0.5, "lines")) +
    # Axis formatting
    scale_x_continuous(
      name = "Total Sequencing Depth (Million Reads)",
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 6)
    ) +
    scale_y_continuous(
      name = "Measured Median Reads per Cell", 
      labels = scales::comma,
      breaks = scales::pretty_breaks(n = 6)
    ) +
    # Labels and theme
    labs(
      title = "Total Reads to Reads per Cell Conversion Model",
      subtitle = "Linear regression fitted to single-cell and single-nucleus data with measured reads per cell",
      caption = "Used for converting total sequencing depth to estimated per-cell depth across all protocols"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      plot.caption = element_text(size = 9, color = "gray60", hjust = 0),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(title = "Technology", override.aes = list(size = 4)),
      shape = guide_legend(title = "Sample Type", override.aes = list(size = 4))
    )
  
  return(p)
}

# Create individual model fit plots
create_individual_fit_plots <- function(fitted_models, rarefaction_data, fitted_curves) {
  plot_list <- list()
  
  for (i in 1:nrow(fitted_models)) {
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    
    if (!model_result$convergence) next
    
    # Get data
    orig_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
    fitted_data <- fitted_curves %>% filter(curve_id == !!curve_id)
    
    # Extract metadata
    technology <- str_extract(curve_id, "^[^_]+")
    sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
    feature_type <- str_extract(curve_id, "[^_]+$")
    
    # Create plot
    p <- ggplot() +
      geom_line(data = fitted_data, 
                aes(x = total_reads_unified, y = featureNum_fitted),
                color = get_color_map()[technology], 
                linetype = LINETYPE_MAP[sample_type],
                linewidth = 1.5, alpha = 0.8) +
      geom_point(data = orig_data,
                 aes(x = total_reads_unified, y = featureNum),
                 color = get_color_map()[technology], 
                 shape = SHAPE_MAP[sample_type],
                 size = 3, alpha = 0.8) +
      labs(title = paste(technology, sample_type, feature_type),
           subtitle = paste("Model:", model_result$model, "| R² =", round(model_result$r_squared, 3)),
           x = "Total Reads",
           y = "Number of Features Detected") +
      theme_minimal() +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma)
    
    plot_list[[curve_id]] <- p
  }
  
  return(plot_list)
}

# Create individual model fit plots with reads per cell on x-axis
create_individual_fit_plots_reads_per_cell <- function(fitted_models, rarefaction_data, fitted_curves, reads_per_cell_converter) {
  plot_list <- list()
  
  for (i in 1:nrow(fitted_models)) {
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    
    if (!model_result$convergence) next
    
    # Get data and convert to reads per cell
    orig_data <- rarefaction_data %>% 
      filter(curve_id == !!curve_id) %>%
      mutate(reads_per_cell = reads_per_cell_converter(total_reads_unified))
    
    fitted_data <- fitted_curves %>% 
      filter(curve_id == !!curve_id) %>%
      mutate(reads_per_cell = reads_per_cell_converter(total_reads_unified))
    
    # Extract metadata
    technology <- str_extract(curve_id, "^[^_]+")
    sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
    feature_type <- str_extract(curve_id, "[^_]+$")
    
    # Create plot
    p <- ggplot() +
      geom_line(data = fitted_data, 
                aes(x = reads_per_cell, y = featureNum_fitted),
                color = get_color_map()[technology], 
                linetype = LINETYPE_MAP[sample_type],
                linewidth = 1.5, alpha = 0.8) +
      geom_point(data = orig_data,
                 aes(x = reads_per_cell, y = featureNum),
                 color = get_color_map()[technology], 
                 shape = SHAPE_MAP[sample_type],
                 size = 3, alpha = 0.8) +
      labs(title = paste(technology, sample_type, feature_type),
           subtitle = paste("Model:", model_result$model, "| R² =", round(model_result$r_squared, 3)),
           x = "Estimated Median Reads per Cell",
           y = "Number of Features Detected") +
      theme_minimal() +
      scale_x_continuous(labels = scales::comma) +
      scale_y_continuous(labels = scales::comma)
    
    plot_list[[curve_id]] <- p
  }
  
  return(plot_list)
}

# Create threshold bar plots
create_threshold_bar_plots <- function(marginal_returns_results) {
  plot_list <- list()
  feature_types <- unique(marginal_returns_results$feature_type)
  
  for (ft in feature_types) {
    ft_data <- marginal_returns_results %>%
      filter(feature_type == ft) %>%
      mutate(
        protocol = paste(technology, sample_type),
        threshold_label = paste0(threshold, " f/M")
      )
    
    # Complete missing combinations
    complete_grid <- expand.grid(
      protocol = unique(ft_data$protocol),
      threshold = MARGINAL_THRESHOLDS,
      stringsAsFactors = FALSE
    )
    
    ft_data <- complete_grid %>%
      left_join(ft_data, by = c("protocol", "threshold")) %>%
      mutate(
        features = ifelse(is.na(features), 0, features),
        technology = ifelse(is.na(technology), 
                           ifelse(grepl("^PacBio", protocol), "PacBio", "ONT"), 
                           technology),
        threshold_label = paste0(threshold, " f/M")
      )
    
    p <- ggplot(ft_data, aes(x = features, y = reorder_within(protocol, features, threshold_label))) +
      geom_bar(aes(fill = technology), stat = "identity", width = 0.7) +
      facet_wrap(~threshold_label, ncol = 1, scales = "free") +
      scale_y_reordered() +
      scale_fill_manual(values = get_color_map()) +
      scale_x_continuous(labels = scales::comma) +
      labs(title = paste(ft, "at Marginal Return Thresholds"),
           x = "Total Features Detected",
           y = "Protocol") +
      theme_minimal()
    
    plot_list[[ft]] <- p
  }
  
  return(plot_list)
}

# Save plots
save_plots <- function(plot_list, filename, width = 12, height = 8) {
  ensure_directory(PLOTS_DIR)
  pdf_path <- file.path(PLOTS_DIR, filename)
  
  pdf(pdf_path, width = width, height = height)
  for (p in plot_list) {
    print(p)
  }
  dev.off()
  
  cat("Saved plots to:", pdf_path, "\n")
}

# Create extended saturation plots
create_extended_saturation_plots <- function(fitted_models, rarefaction_data, marginal_returns_results, reads_per_cell_converter) {
  plot_list <- list()
  
  for (i in 1:nrow(fitted_models)) {
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    
    if (!model_result$convergence) next
    
    # Get threshold info for this curve
    thresh_info <- marginal_returns_results %>% filter(curve_id == !!curve_id)
    if (nrow(thresh_info) == 0) next
    
    # Get original data
    orig_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
    
    # Extract metadata
    technology <- str_extract(curve_id, "^[^_]+")
    sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
    feature_type <- str_extract(curve_id, "[^_]+$")
    
    # Determine extension range - use dynamic threshold selection
    extend_to_x <- max(orig_data$total_reads_unified, na.rm = TRUE) * 3  # Conservative fallback
    extension_label <- "Extended 3x"
    
    # Try to use configured extension threshold, fall back to dynamic selection
    extension_thresh <- find_extension_threshold(thresh_info, MARGINAL_THRESHOLDS)
    if (nrow(extension_thresh) > 0) {
      extend_to_x <- extension_thresh$depth[1] * 1.3  # Extend 30% past selected threshold
      extension_label <- paste0("Extended to\n", extension_thresh$threshold[1], " f/M")
    } else {
      fold_increase <- round(extend_to_x / max(orig_data$total_reads_unified, na.rm = TRUE), 1)
      extension_label <- paste0("Extended ", fold_increase, "x")
    }
    
    # Create extended x range
    x_extended <- seq(from = min(orig_data$total_reads_unified), to = extend_to_x, length.out = 300)
    
    # Generate extended predictions
    y_extended <- predict(model_result$fit_object, newdata = data.frame(x = x_extended))
    
    # Create extended curve data with reads per cell conversion
    extended_data <- tibble(
      total_reads_unified = x_extended,
      featureNum_fitted = y_extended,
      reads_per_cell = reads_per_cell_converter(x_extended)
    )
    
    # Also convert original data
    orig_data_converted <- orig_data %>%
      mutate(reads_per_cell = reads_per_cell_converter(total_reads_unified))
    
    # Create title
    model_type_display <- case_when(
      model_result$model == "michaelis_menten" ~ "Michaelis-Menten",
      model_result$model == "asymptotic_exp" ~ "Asymptotic Exponential",
      model_result$model == "power_law" ~ "Power Law",
      model_result$model == "logarithmic" ~ "Logarithmic",
      model_result$model == "hill" ~ "Hill Equation",
      TRUE ~ model_result$model
    )
    
    plot_title <- paste0(technology, " ", 
                         case_when(
                           sample_type == "SC" ~ "Single-Cell",
                           sample_type == "SN" ~ "Single-Nucleus",
                           TRUE ~ sample_type
                         ), " ", 
                         feature_type)
    
    # Create subtitle based on extension method used
    if (nrow(extension_thresh) > 0) {
      subtitle <- paste0(model_type_display, " | Extended to ", extension_thresh$threshold[1], " f/M Threshold")
    } else {
      fold_increase <- round(extend_to_x / max(orig_data$total_reads_unified, na.rm = TRUE), 1)
      subtitle <- paste0(model_type_display, " | Extended ", fold_increase, "x (no suitable threshold)")
    }
    
    # Create the plot
    p <- ggplot() +
      # Add extended fitted curve
      geom_line(data = extended_data, 
                aes(x = reads_per_cell, y = featureNum_fitted),
                color = get_color_map()[technology], 
                linetype = LINETYPE_MAP[sample_type],
                linewidth = 1.2, alpha = 0.7) +
      # Add original data points
      geom_point(data = orig_data_converted,
                 aes(x = reads_per_cell, y = featureNum),
                 color = get_color_map()[technology], 
                 shape = SHAPE_MAP[sample_type],
                 size = 2.5, alpha = 0.9) +
      # Add vertical line at current max depth (converted)
      geom_vline(xintercept = reads_per_cell_converter(max_x_orig), linetype = "dashed", alpha = 0.5, color = "gray50") +
      # Add threshold lines if 1 f/M threshold exists (converted)
      {if (nrow(thresh_1fm) > 0) {
        list(
          geom_vline(xintercept = reads_per_cell_converter(thresh_1fm$depth[1]), linetype = "dotted", alpha = 0.7, color = "darkgreen"),
          geom_hline(yintercept = thresh_1fm$features[1], linetype = "dotted", alpha = 0.5, color = "darkgreen")
        )
      }} +
      # Annotate regions (converted positions)
      annotate("text", x = reads_per_cell_converter(max_x_orig) * 0.5, y = max(extended_data$featureNum_fitted) * 0.95, 
               label = "Current\nData", size = 3, alpha = 0.7, hjust = 0.5) +
      annotate("text", x = reads_per_cell_converter(extend_to_x) * 0.8, y = max(extended_data$featureNum_fitted) * 0.95, 
               label = extension_label, size = 3, alpha = 0.7, hjust = 0.5) +
      labs(title = plot_title,
           subtitle = subtitle,
           x = "Estimated Median Reads per Cell",
           y = "Number of Features Detected") +
      theme_minimal() +
      theme(
        text = element_text(size = 10),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray50"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_continuous(
        labels = scales::comma, 
        breaks = scales::pretty_breaks(n = 5)
      ) +
      scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 5))
    
    plot_list[[curve_id]] <- p
  }
  
  return(plot_list)
}

#' Add TD75/TD50 Markers to Rarefaction Plots
#'
#' Adds TD75 and TD50 markers to existing rarefaction curve plots for easy
#' identification of transcript discovery efficiency benchmarks.
#'
#' @param plot_object ggplot object to add markers to
#' @param td_metrics_results Data frame containing TD75/TD50 results for the specific curve
#' @param reads_per_cell_converter Function to convert total reads to reads per cell
#' @param curve_id Character string identifying the current curve
#' @param show_labels Logical indicating whether to show TD75/TD50 labels
#'
#' @return Modified ggplot object with TD75/TD50 markers
#'
#' @details
#' **Marker Design:**
#' - TD75: Solid vertical line in orange with triangle marker
#' - TD50: Dashed vertical line in blue with circle marker  
#' - Optional labels with depth values and confidence intervals
#' - Markers are positioned at estimated sequencing depths
#'
#' **Visual Elements:**
#' - Vertical lines showing required sequencing depths
#' - Point markers at intersection with fitted curve
#' - Optional text annotations with confidence intervals
#' - Color coding: TD75 (orange), TD50 (blue)
#'
#' **Usage Context:**
#' Called by plotting functions to add TD markers to rarefaction plots
#' for practical interpretation of sequencing depth requirements.
#'
#' @examples
#' # Add TD markers to a rarefaction plot
#' plot_with_td <- add_td_markers_to_plot(
#'   base_plot,
#'   td_metrics_results,
#'   reads_per_cell_converter,
#'   "ONT_SC_Genes Discovery"
#' )
add_td_markers_to_plot <- function(plot_object, td_metrics_results, reads_per_cell_converter, curve_id, show_labels = TRUE) {
  # Filter TD results for this specific curve
  curve_td_data <- td_metrics_results %>%
    filter(curve_id == !!curve_id)
  
  if (nrow(curve_td_data) == 0) {
    return(plot_object)  # No TD data for this curve, return unmodified plot
  }
  
  # Define TD marker styling
  td_colors <- c("TD75" = "#FF8C00", "TD50" = "#1E90FF")  # Orange and blue
  td_linetypes <- c("TD75" = "solid", "TD50" = "dashed")
  td_shapes <- c("TD75" = 17, "TD50" = 16)  # Triangle and circle
  
  # Add markers for each TD metric
  for (i in 1:nrow(curve_td_data)) {
    row <- curve_td_data[i, ]
    metric <- row$metric
    depth <- row$estimated_reads
    target_features <- row$target_features
    
    # Convert to reads per cell for consistency
    depth_rpc <- reads_per_cell_converter(depth)
    
    # Create confidence interval text
    ci_text <- if (!is.na(row$lower_ci) && !is.na(row$upper_ci)) {
      paste0("\n[", scales::comma(round(reads_per_cell_converter(row$lower_ci))), 
             "-", scales::comma(round(reads_per_cell_converter(row$upper_ci))), "]")
    } else {
      ""
    }
    
    # Add vertical line at TD depth
    plot_object <- plot_object +
      geom_vline(
        xintercept = depth_rpc,
        color = td_colors[metric],
        linetype = td_linetypes[metric],
        linewidth = 1.2,
        alpha = 0.8
      )
    
    # Add point marker at intersection with curve
    plot_object <- plot_object +
      annotate(
        "point",
        x = depth_rpc,
        y = target_features,
        color = td_colors[metric],
        shape = td_shapes[metric],
        size = 4,
        alpha = 0.9
      )
    
    # Add text label if requested
    if (show_labels) {
      label_text <- paste0(
        metric, ": ", scales::comma(round(depth_rpc)), 
        " reads/cell", ci_text
      )
      
      # Position label at top of plot area
      plot_object <- plot_object +
        annotate(
          "text",
          x = depth_rpc,
          y = Inf,
          label = label_text,
          color = td_colors[metric],
          vjust = if (metric == "TD75") 1.2 else 2.5,  # Stack TD50 below TD75
          hjust = 0.5,
          size = 3,
          fontface = "bold",
          alpha = 0.8
        )
    }
  }
  
  return(plot_object)
}

#' Create TD75/TD50 Enhanced Rarefaction Plots
#'
#' Creates rarefaction curve plots with integrated TD75/TD50 markers for 
#' practical interpretation of sequencing efficiency benchmarks.
#'
#' @param fitted_models Data frame containing fitted model results
#' @param rarefaction_data Original rarefaction curve data
#' @param td_metrics_results Data frame containing TD75/TD50 calculations
#' @param reads_per_cell_converter Function to convert total reads to reads per cell
#'
#' @return List of ggplot objects with TD markers
#'
#' @details
#' **Enhanced Features:**
#' - Standard rarefaction curves with fitted models
#' - TD75/TD50 vertical markers showing required depths
#' - Confidence intervals from bootstrap analysis
#' - Reads per cell conversion for practical interpretation
#' - Protocol comparison across technologies and sample types
#'
#' **Visual Design:**
#' - Consistent color scheme with existing ORACLE-seq plots
#' - Clear TD marker differentiation (TD75 orange, TD50 blue)
#' - Integrated legends and annotations
#' - Publication-quality formatting
#'
#' **Output Format:**
#' Returns named list of plots suitable for save_plots() function.
#' Each plot focuses on a specific protocol with TD markers.
#'
#' @examples
#' # Create TD-enhanced plots
#' td_plots <- create_td_enhanced_rarefaction_plots(
#'   fitted_models_with_uncertainty,
#'   rarefaction_data,
#'   td_metrics_results,
#'   reads_per_cell_converter
#' )
create_td_enhanced_rarefaction_plots <- function(fitted_models, rarefaction_data, td_metrics_results, reads_per_cell_converter) {
  plot_list <- list()
  
  cat("Creating TD75/TD50 enhanced rarefaction plots...\n")
  
  # Process each fitted model
  for (i in 1:nrow(fitted_models)) {
    if (!fitted_models$convergence[i]) next
    
    curve_id <- fitted_models$curve_id[i]
    model_result <- fitted_models$fitted_model[[i]]
    
    cat("  Creating TD-enhanced plot for:", curve_id, "\n")
    
    # Get original data for this curve
    orig_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
    
    # Create extended prediction range
    x_max_orig <- max(orig_data$total_reads_unified)
    x_extend <- x_max_orig * 3  # Extend 3x for TD visualization
    x_seq <- seq(min(orig_data$total_reads_unified), x_extend, length.out = 200)
    
    # Generate fitted curve predictions
    fitted_predictions <- predict(model_result$fit_object, newdata = data.frame(x = x_seq))
    
    # Create fitted curve data with reads per cell conversion
    fitted_data <- tibble(
      total_reads_unified = x_seq,
      featureNum_fitted = fitted_predictions,
      reads_per_cell = reads_per_cell_converter(x_seq),
      curve_type = ifelse(x_seq <= x_max_orig, "observed_range", "extrapolated")
    )
    
    # Convert original data to reads per cell
    orig_data_converted <- orig_data %>%
      mutate(reads_per_cell = reads_per_cell_converter(total_reads_unified))
    
    # Extract metadata for plot labels
    technology <- str_extract(curve_id, "^[^_]+")
    parts <- str_split(curve_id, "_")[[1]]
    sample_type <- if (length(parts) >= 2) parts[2] else "Unknown"
    feature_type <- str_remove(paste(parts[3:length(parts)], collapse = " "), " Discovery")
    
    # Create plot title
    sample_display <- case_when(
      sample_type == "SC" ~ "Single-Cell",
      sample_type == "SN" ~ "Single-Nucleus",
      TRUE ~ sample_type
    )
    plot_title <- paste0(technology, " ", sample_display, " ", feature_type, " Discovery")
    
    # Model info for subtitle
    model_display <- case_when(
      model_result$model == "michaelis_menten" ~ "Michaelis-Menten",
      model_result$model == "asymptotic_exp" ~ "Asymptotic Exponential", 
      model_result$model == "power_law" ~ "Power Law",
      model_result$model == "logarithmic" ~ "Logarithmic",
      model_result$model == "hill" ~ "Hill Equation",
      TRUE ~ str_to_title(str_replace_all(model_result$model, "_", " "))
    )
    
    # Base plot
    p <- ggplot() +
      # Fitted curve (different styles for observed vs extrapolated)
      geom_line(
        data = fitted_data %>% filter(curve_type == "observed_range"),
        aes(x = reads_per_cell, y = featureNum_fitted),
        color = get_color_map()[technology],
        linewidth = 1.5,
        alpha = 0.9
      ) +
      geom_line(
        data = fitted_data %>% filter(curve_type == "extrapolated"),
        aes(x = reads_per_cell, y = featureNum_fitted),
        color = get_color_map()[technology],
        linewidth = 1.2,
        linetype = "dashed",
        alpha = 0.7
      ) +
      # Original data points
      geom_point(
        data = orig_data_converted,
        aes(x = reads_per_cell, y = featureNum),
        color = get_color_map()[technology],
        shape = SHAPE_MAP[sample_type],
        size = 3,
        alpha = 0.9
      ) +
      # Vertical line at current max depth
      geom_vline(
        xintercept = reads_per_cell_converter(x_max_orig),
        linetype = "dotted",
        color = "gray50",
        alpha = 0.6
      ) +
      # Labels and theme
      labs(
        title = plot_title,
        subtitle = paste0(model_display, " Model | R² = ", round(fitted_models$r_squared[i], 3), " | TD75/TD50 Analysis"),
        x = "Estimated Median Reads per Cell",
        y = paste0("Number of ", feature_type, " Detected"),
        caption = "Solid line: observed data range | Dashed line: model extrapolation | TD markers: transcript discovery benchmarks"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray50"),
        plot.caption = element_text(size = 8, color = "gray60", hjust = 0),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      ) +
      scale_x_continuous(
        labels = scales::comma,
        breaks = scales::pretty_breaks(n = 6)
      ) +
      scale_y_continuous(
        labels = scales::comma,
        breaks = scales::pretty_breaks(n = 6)
      )
    
    # Add TD75/TD50 markers
    p <- add_td_markers_to_plot(p, td_metrics_results, reads_per_cell_converter, curve_id, show_labels = TRUE)
    
    # Add annotation for observed vs extrapolated regions
    p <- p +
      annotate(
        "text",
        x = reads_per_cell_converter(x_max_orig * 0.7),
        y = max(fitted_data$featureNum_fitted) * 0.95,
        label = "Observed\nData",
        size = 3,
        alpha = 0.7,
        hjust = 0.5
      ) +
      annotate(
        "text", 
        x = reads_per_cell_converter(x_extend * 0.8),
        y = max(fitted_data$featureNum_fitted) * 0.95,
        label = "Model\nExtrapolation",
        size = 3,
        alpha = 0.7,
        hjust = 0.5
      )
    
    plot_list[[curve_id]] <- p
  }
  
  cat("Created", length(plot_list), "TD-enhanced rarefaction plots\n")
  return(plot_list)
}

# Create combined extended saturation plots
create_combined_extended_saturation_plots <- function(fitted_models, rarefaction_data, marginal_returns_results, reads_per_cell_converter) {
  plot_list <- list()
  
  # Group by feature type  
  feature_types <- c("Genes Discovery", "Isoforms Discovery")
  
  for (ft in feature_types) {
    cat("Processing combined plots for:", ft, "\n")
    
    # Find the best available threshold for extension across all curves for this feature type
    # Use the dynamic threshold selection utility
    all_thresholds_for_ft <- marginal_returns_results %>%
      filter(str_detect(curve_id, paste0(ft, "$")))
    
    if (nrow(all_thresholds_for_ft) == 0) {
      # Fall back to conservative extension
      first_curve_data <- rarefaction_data %>% 
        filter(str_detect(curve_id, paste0(ft, "$"))) %>%
        slice(1)
      max_extend_to <- max(first_curve_data$total_reads_unified, na.rm = TRUE) * 5
      extension_method <- "conservative fallback"
    } else {
      # Find best threshold using dynamic selection
      best_thresh <- find_extension_threshold(all_thresholds_for_ft, MARGINAL_THRESHOLDS)
      if (nrow(best_thresh) > 0) {
        # Use the maximum depth for the selected threshold across all protocols for this feature type
        max_extend_to <- max(best_thresh$depth, na.rm = TRUE) * 1.3  # Extend 30% past furthest threshold
        extension_method <- paste0(best_thresh$threshold[1], " f/M threshold")
      } else {
        # Fallback to median threshold available
        available_thresholds <- unique(all_thresholds_for_ft$threshold)
        if (length(available_thresholds) > 0) {
          median_thresh <- median(available_thresholds)
          median_thresh_data <- all_thresholds_for_ft %>% filter(threshold == median_thresh)
          max_extend_to <- max(median_thresh_data$depth, na.rm = TRUE) * 1.3
          extension_method <- paste0(median_thresh, " f/M threshold (median)")
        } else {
          first_curve_data <- rarefaction_data %>% 
            filter(str_detect(curve_id, paste0(ft, "$"))) %>%
            slice(1)
          max_extend_to <- max(first_curve_data$total_reads_unified, na.rm = TRUE) * 5
          extension_method <- "conservative fallback"
        }
      }
    }
    
    # Get all curves for this feature type
    curves_for_ft <- fitted_models %>% 
      filter(str_detect(curve_id, paste0(ft, "$")), convergence == TRUE)
    
    if (nrow(curves_for_ft) == 0) next

    # Create extended curves for each protocol
    all_extended_data <- data.frame()
    all_original_data <- data.frame()
    
    for (i in 1:nrow(curves_for_ft)) {
      curve_id <- curves_for_ft$curve_id[i]
      model_result <- curves_for_ft$fitted_model[[i]]
      
      # Extract metadata
      technology <- str_extract(curve_id, "^[^_]+")
      sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
      
      # Get original data
      orig_data <- rarefaction_data %>% filter(curve_id == !!curve_id)
      orig_data$technology <- technology
      orig_data$sample_type <- sample_type
      orig_data$protocol <- paste(technology, sample_type)
      
      # Create extended x range
      min_x <- min(orig_data$total_reads_unified, na.rm = TRUE)
      x_extended <- seq(from = min_x, to = max_extend_to, length.out = 500)
      
      # Generate extended predictions
      y_extended <- predict(model_result$fit_object, newdata = data.frame(x = x_extended))
      
      # Create extended curve data
      extended_data <- data.frame(
        curve_id = curve_id,
        total_reads_unified = x_extended,
        featureNum_fitted = y_extended,
        technology = technology,
        sample_type = sample_type,
        protocol = paste(technology, sample_type)
      )
      
      all_extended_data <- rbind(all_extended_data, extended_data)
      all_original_data <- rbind(all_original_data, orig_data)
    }
    
    # Create the combined plot
    feature_name <- str_remove(ft, " Discovery")
    
    # Create subtitle based on extension method used
    if (extension_method == "conservative fallback") {
      max_current_depth <- max(all_original_data$total_reads_unified, na.rm = TRUE)
      fold_increase <- round(max_extend_to / max_current_depth, 1)
      subtitle_text <- paste("Extended", fold_increase, "x current depth (no suitable thresholds found)")
    } else {
      subtitle_text <- paste("Extended to", extension_method, "(furthest curve)")
    }
    
    p <- ggplot() +
      # Add extended fitted curves
      geom_line(data = all_extended_data, 
                aes(x = total_reads_unified, y = featureNum_fitted,
                    color = technology, linetype = sample_type, group = curve_id),
                linewidth = 1.2, alpha = 0.8) +
      # Add original data points
      geom_point(data = all_original_data,
                 aes(x = total_reads_unified, y = featureNum,
                     color = technology, shape = sample_type),
                 size = 2, alpha = 0.7) +
      # Color and line mappings
      scale_color_manual(values = get_color_map(), name = "Technology") +
      scale_linetype_manual(values = LINETYPE_MAP, name = "Sample Type",
                           labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      scale_shape_manual(values = SHAPE_MAP, name = "Sample Type",
                        labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      # Axis formatting
      scale_x_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 6)) +
      scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 6)) +
      # Labels and theme
      labs(title = paste("Extended Saturation Curves:", feature_name),
           subtitle = subtitle_text,
           x = "Total Reads",
           y = paste("Number of", feature_name, "Detected")) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray50"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      ) +
      guides(
        color = guide_legend(title = "Technology", override.aes = list(linewidth = 2)),
        linetype = guide_legend(title = "Sample Type", override.aes = list(linewidth = 2)),
        shape = guide_legend(title = "Sample Type")
      )
    
    # Create unique plot name
    plot_name <- paste0(ft, "_Extended")
    plot_list[[plot_name]] <- p
  }
  
  return(plot_list)
}

# Create reads threshold bar plots
create_reads_threshold_bar_plots <- function(marginal_returns_results, reads_per_cell_converter) {
  plot_list <- list()
  
  # Group by feature type
  feature_types <- unique(marginal_returns_results$feature_type)
  
  for (ft in feature_types) {
    cat("Creating reads threshold bar plots for:", ft, "\n")
    
    # Get data for this feature type
    ft_data <- marginal_returns_results %>%
      filter(feature_type == ft) %>%
      mutate(
        protocol = paste(technology, 
                         case_when(
                           sample_type == "SC" ~ "Single-Cell",
                           sample_type == "SN" ~ "Single-Nucleus",
                           TRUE ~ sample_type
                         )),
        threshold_label = paste0(threshold, " ", str_remove(ft, " Discovery"), " per Million Reads"),
        # Convert depth to millions for better readability
        depth_millions = depth / 1e6
      )
    
    if (nrow(ft_data) == 0) next
    
    # Create complete data frame with all protocol-threshold combinations
    all_protocols <- unique(ft_data$protocol)
    all_thresholds <- unique(ft_data$threshold)
    
    complete_grid <- expand.grid(
      protocol = all_protocols,
      threshold = all_thresholds,
      stringsAsFactors = FALSE
    ) %>%
      left_join(ft_data, 
                by = c("protocol", "threshold")) %>%
      # Fill missing threshold_label and depth_millions
      mutate(
        # Use dynamic threshold labeling instead of hardcoded case_when
        threshold_label = create_dynamic_threshold_labels(ft, threshold, "full")[match(threshold, MARGINAL_THRESHOLDS)],
        # Keep depth as NA for missing combinations
        depth_millions = ifelse(is.na(depth_millions), NA, depth_millions),
        # Extract technology from protocol name for missing entries
        technology = ifelse(is.na(technology), 
                           ifelse(grepl("^PacBio", protocol), "PacBio", "ONT"), 
                           technology),
        # Add model stability columns
        model_stability = ifelse(is.na(model_stability), "Unknown", model_stability),
        model_stability = factor(model_stability, levels = c("High", "Moderate", "Low", "Very Low", "Unknown"))
      )
    
    ft_data <- complete_grid
    
    feature_name <- str_remove(ft, " Discovery")
    
    # Use dynamic factor ordering instead of hardcoded levels
    ft_data$threshold_label <- factor(ft_data$threshold_label, 
                                     levels = create_threshold_factor_levels(ft, MARGINAL_THRESHOLDS, "desc"))
    
    # Remove rows with missing depth
    ft_data <- ft_data %>%
      filter(!is.na(depth_millions))
    
    # Add estimated reads per cell
    ft_data <- ft_data %>%
      mutate(reads_per_cell = reads_per_cell_converter(depth_millions * 1e6))
    
    # Create Total Reads plot (3 thresholds stacked vertically)
    total_reads_plot <- ggplot(ft_data, aes(x = depth_millions, y = reorder_within(protocol, depth_millions, threshold_label))) +
      # Add error bars first (behind bars)
      geom_errorbarh(aes(xmin = depth_ci_lower/1e6, xmax = depth_ci_upper/1e6),
                     height = 0.3, alpha = 0.6, color = "gray30") +
      # Add bars with model stability transparency
      geom_bar(aes(fill = technology, alpha = model_stability), stat = "identity", width = 0.7) +
      facet_wrap(~threshold_label, ncol = 1, scales = "free") +
      scale_y_reordered() +
      scale_fill_manual(values = get_color_map(), name = "Technology") +
      # Add alpha scale for model stability
      scale_alpha_manual(values = c("High" = 1, "Moderate" = 0.8, 
                                   "Low" = 0.6, "Very Low" = 0.4, "Unknown" = 0.7),
                        name = "Model Stability",
                        labels = c("High (CV < 5%)", "Moderate (CV 5-15%)", "Low (CV 15-30%)", 
                                  "Very Low (CV > 30%)", "Unknown")) +
      scale_x_continuous(
        name = "Sequencing Depth Required (Million Reads)",
        labels = scales::comma, 
        breaks = scales::pretty_breaks(n = 8)
      ) +
      labs(
        title = paste("Sequencing Depth Required:", feature_name, "- Total Reads"),
        subtitle = "Error bars show 95% bootstrap confidence intervals | Transparency indicates model stability",
        y = "Protocol"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.25)
      )
    
    # Create Estimated Reads per Cell plot (3 thresholds stacked vertically)
    reads_per_cell_plot <- ggplot(ft_data, aes(x = reads_per_cell, y = reorder_within(protocol, reads_per_cell, threshold_label))) +
      # Add error bars first (behind bars) - convert depth CI to reads per cell
      geom_errorbarh(aes(xmin = reads_per_cell_converter(depth_ci_lower), 
                         xmax = reads_per_cell_converter(depth_ci_upper)),
                     height = 0.3, alpha = 0.6, color = "gray30") +
      # Add bars with model stability transparency
      geom_bar(aes(fill = technology, alpha = model_stability), stat = "identity", width = 0.7) +
      facet_wrap(~threshold_label, ncol = 1, scales = "free") +
      scale_y_reordered() +
      scale_fill_manual(values = get_color_map(), name = "Technology") +
      # Add alpha scale for model stability
      scale_alpha_manual(values = c("High" = 1, "Moderate" = 0.8, 
                                   "Low" = 0.6, "Very Low" = 0.4, "Unknown" = 0.7),
                        name = "Model Stability",
                        labels = c("High (CV < 5%)", "Moderate (CV 5-15%)", "Low (CV 15-30%)", 
                                  "Very Low (CV > 30%)", "Unknown")) +
      scale_x_continuous(
        name = "Estimated Median Reads per Cell",
        labels = scales::comma,
        breaks = scales::pretty_breaks(n = 8)
      ) +
      labs(
        title = paste("Sequencing Depth Required:", feature_name, "- Estimated Reads per Cell"),
        subtitle = "Error bars show 95% bootstrap confidence intervals | Transparency indicates model stability",
        y = "Protocol"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.25)
      )
    
    # Add both plots to the list
    plot_list[[paste0(ft, "_TotalReads")]] <- total_reads_plot
    plot_list[[paste0(ft, "_ReadsPerCell")]] <- reads_per_cell_plot
  }
  
  return(plot_list)
}

# Create marginal returns derivative plots  
create_marginal_returns_derivative_plots <- function(fitted_models, rarefaction_data, marginal_returns_results, reads_per_cell_converter) {
  plot_list <- list()
  
  # Group by feature type
  feature_types <- c("Genes Discovery", "Isoforms Discovery")
  
  for (ft in feature_types) {
    cat("Creating marginal returns derivative plots for:", ft, "\n")
    
    # Get all curves for this feature type
    curves_for_ft <- fitted_models %>% 
      filter(str_detect(curve_id, paste0(ft, "$")), convergence == TRUE)
    
    if (nrow(curves_for_ft) == 0) next
    
    # Collect all derivative data
    all_derivative_data <- data.frame()
    all_threshold_data <- data.frame()
    
    for (i in 1:nrow(curves_for_ft)) {
      curve_id <- curves_for_ft$curve_id[i]
      model_result <- curves_for_ft$fitted_model[[i]]
      
      # Extract metadata
      technology <- str_extract(curve_id, "^[^_]+")
      sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
      
      # Calculate marginal returns using the new function
      derivative_result <- calculate_marginal_returns_simple(model_result)
      
      if (nrow(derivative_result) > 0) {
        derivative_result$curve_id <- curve_id
        derivative_result$technology <- technology
        derivative_result$sample_type <- sample_type
        derivative_result$protocol <- paste(technology, 
                                           case_when(
                                             sample_type == "SC" ~ "Single-Cell",
                                             sample_type == "SN" ~ "Single-Nucleus",
                                             TRUE ~ sample_type
                                           ))
        
        all_derivative_data <- rbind(all_derivative_data, derivative_result)
      }
      
      # Get threshold depths for this curve
      curve_thresholds <- marginal_returns_results %>% filter(curve_id == !!curve_id)
      if (nrow(curve_thresholds) > 0) {
        curve_thresholds$technology <- technology
        curve_thresholds$sample_type <- sample_type
        curve_thresholds$protocol <- paste(technology, 
                                          case_when(
                                            sample_type == "SC" ~ "Single-Cell",
                                            sample_type == "SN" ~ "Single-Nucleus",
                                            TRUE ~ sample_type
                                          ))
        all_threshold_data <- rbind(all_threshold_data, curve_thresholds)
      }
    }
    
    if (nrow(all_derivative_data) == 0) next
    
    feature_name <- str_remove(ft, " Discovery")
    
    # Create Total Reads plot
    total_reads_plot <- ggplot() +
      # Add derivative curves
      geom_line(data = all_derivative_data,
                aes(x = x, y = marginal_returns, 
                    color = technology, linetype = sample_type, group = curve_id),
                linewidth = 1.1, alpha = 0.8) +
      # Add threshold lines if available
      {if (nrow(all_threshold_data) > 0) {
        # Use dynamic threshold annotations instead of hardcoded ones
        threshold_annotations <- create_threshold_annotations(MARGINAL_THRESHOLDS, x_position = Inf, hjust = 1.1)
        c(
          # Add threshold points on curves
          list(geom_point(data = all_threshold_data,
                         aes(x = depth, y = threshold, color = technology, shape = sample_type),
                         size = 3, alpha = 0.8, stroke = 1.2, fill = "white")),
          threshold_annotations
        )
      }} +
      # Color and line mappings
      scale_color_manual(values = get_color_map(), name = "Technology") +
      scale_linetype_manual(values = LINETYPE_MAP, name = "Sample Type",
                           labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      scale_shape_manual(values = SHAPE_MAP, name = "Sample Type",
                        labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      # Axis formatting - log scales
      scale_x_log10(
        name = "Total Reads",
        labels = function(x) ifelse(x >= 1e6, paste0(x/1e6, "M"), paste0(x/1e3, "K")),
        breaks = c(1e4, 1e5, 1e6, 1e7, 1e8, 1e9)
      ) +
      scale_y_log10(
        name = paste("Marginal Returns (New", feature_name, "per Million Reads)"),
        labels = function(y) ifelse(y >= 1, as.character(y), paste0(y*1000, "m")),
        breaks = c(0.01, 0.1, 1, 10, 100, 1000)
      ) +
      # Labels and theme
      labs(
        title = paste("Marginal Returns Analysis:", feature_name, "(Total Reads)"),
        subtitle = "Curves show rate of feature discovery (derivative) | Points mark marginal return thresholds",
        caption = paste("Log-log scale | Horizontal lines show", 
                       paste(MARGINAL_THRESHOLDS, collapse = ", "), 
                       "features per million reads thresholds")
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray50"),
        plot.caption = element_text(size = 9, color = "gray60"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      ) +
      guides(
        color = guide_legend(title = "Technology", override.aes = list(linewidth = 2)),
        linetype = guide_legend(title = "Sample Type", override.aes = list(linewidth = 2)),
        shape = guide_legend(title = "Sample Type")
      )
    
    # Create reads per cell version
    all_derivative_data_rpc <- all_derivative_data %>%
      mutate(reads_per_cell = reads_per_cell_converter(x))
    
    all_threshold_data_rpc <- all_threshold_data %>%
      mutate(reads_per_cell = reads_per_cell_converter(depth))
    
    reads_per_cell_plot <- ggplot() +
      # Add derivative curves (converted to reads per cell)
      geom_line(data = all_derivative_data_rpc,
                aes(x = reads_per_cell, y = marginal_returns, 
                    color = technology, linetype = sample_type, group = curve_id),
                linewidth = 1.1, alpha = 0.8) +
      # Add threshold lines if available
      {if (nrow(all_threshold_data_rpc) > 0) {
        # Use dynamic threshold annotations instead of hardcoded ones
        threshold_annotations <- create_threshold_annotations(MARGINAL_THRESHOLDS, x_position = Inf, hjust = 1.1)
        c(
          # Add threshold points on curves (converted to reads per cell)
          list(geom_point(data = all_threshold_data_rpc,
                         aes(x = reads_per_cell, y = threshold, color = technology, shape = sample_type),
                         size = 3, alpha = 0.8, stroke = 1.2, fill = "white")),
          threshold_annotations
        )
      }} +
      # Color and line mappings
      scale_color_manual(values = get_color_map(), name = "Technology") +
      scale_linetype_manual(values = LINETYPE_MAP, name = "Sample Type",
                           labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      scale_shape_manual(values = SHAPE_MAP, name = "Sample Type",
                        labels = c("Bulk" = "Bulk", "SC" = "Single-Cell", "SN" = "Single-Nucleus")) +
      # Axis formatting - log scales
      scale_x_log10(
        name = "Estimated Median Reads per Cell",
        labels = scales::comma,
        breaks = scales::trans_breaks("log10", function(x) 10^x)
      ) +
      scale_y_log10(
        name = paste("Marginal Returns (New", feature_name, "per Million Reads)"),
        labels = function(y) ifelse(y >= 1, as.character(y), paste0(y*1000, "m")),
        breaks = c(0.01, 0.1, 1, 10, 100, 1000)
      ) +
      # Labels and theme
      labs(
        title = paste("Marginal Returns Analysis:", feature_name, "(Reads per Cell)"),
        subtitle = "Curves show rate of feature discovery (derivative) | Points mark marginal return thresholds",
        caption = "Log-log scale | Horizontal lines show 10, 1, and 0.1 features per million reads thresholds"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray50"),
        plot.caption = element_text(size = 9, color = "gray60"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      ) +
      guides(
        color = guide_legend(title = "Technology", override.aes = list(linewidth = 2)),
        linetype = guide_legend(title = "Sample Type", override.aes = list(linewidth = 2)),
        shape = guide_legend(title = "Sample Type")
      )
    
    # Add both plots to the list
    plot_list[[paste0(ft, "_MarginalReturns_TotalReads")]] <- total_reads_plot
    plot_list[[paste0(ft, "_MarginalReturns_ReadsPerCell")]] <- reads_per_cell_plot
  }
  
  return(plot_list)
}

# Create bootstrap distribution plots
create_bootstrap_distribution_plots <- function(fitted_models, marginal_returns_results) {
  plot_list <- list()
  
  # Get curves with bootstrap results
  curves_with_bootstrap <- fitted_models %>%
    filter(map_lgl(fitted_model, ~ !is.null(.x$bootstrap_results)))
  
  if (nrow(curves_with_bootstrap) == 0) {
    cat("No bootstrap results found for distribution plots\n")
    return(plot_list)
  }
  
  # Group by feature type
  feature_types <- c("Genes Discovery", "Isoforms Discovery")
  
  for (ft in feature_types) {
    cat("Creating bootstrap distribution plots for:", ft, "\n")
    
    # Get curves for this feature type
    ft_curves <- curves_with_bootstrap %>%
      filter(str_detect(curve_id, paste0(ft, "$")))
    
    if (nrow(ft_curves) == 0) next
    
    all_bootstrap_data <- data.frame()
    
    for (i in 1:nrow(ft_curves)) {
      curve_id <- ft_curves$curve_id[i]
      model_result <- ft_curves$fitted_model[[i]]
      bootstrap_results <- model_result$bootstrap_results
      
      if (is.null(bootstrap_results) || length(bootstrap_results) == 0) next
      
      # Extract metadata
      technology <- str_extract(curve_id, "^[^_]+")
      sample_type <- str_extract(curve_id, "(?<=_)[^_]+(?=_)")
      protocol <- paste(technology, 
                       case_when(
                         sample_type == "SC" ~ "Single-Cell",
                         sample_type == "SN" ~ "Single-Nucleus",
                         TRUE ~ sample_type
                       ))
      
      # Convert bootstrap results to data frame format
      for (thresh in names(bootstrap_results)) {
        boot_values <- bootstrap_results[[thresh]]$bootstrap_depths
        if (!is.null(boot_values) && length(boot_values) > 0) {
          bootstrap_df <- data.frame(
            curve_id = curve_id,
            technology = technology,
            sample_type = sample_type,
            protocol = protocol,
            threshold = as.numeric(thresh),
            threshold_label = paste0(thresh, " f/M"),
            bootstrap_depth = boot_values,
            stringsAsFactors = FALSE
          )
          all_bootstrap_data <- rbind(all_bootstrap_data, bootstrap_df)
        }
      }
    }
    
    if (nrow(all_bootstrap_data) == 0) next
    
    # Remove infinite or invalid values
    all_bootstrap_data <- all_bootstrap_data %>%
      filter(is.finite(bootstrap_depth), bootstrap_depth > 0)
    
    if (nrow(all_bootstrap_data) == 0) next
    
    feature_name <- str_remove(ft, " Discovery")
    
    # Get point estimates for vertical lines
    point_estimates <- marginal_returns_results %>%
      filter(str_detect(curve_id, paste0(ft, "$"))) %>%
      mutate(
        protocol = paste(technology, 
                        case_when(
                          sample_type == "SC" ~ "Single-Cell",
                          sample_type == "SN" ~ "Single-Nucleus",
                          TRUE ~ sample_type
                        )),
        threshold_label = paste0(threshold, " f/M")
      )
    
    # Create the bootstrap distribution plot
    p <- ggplot() +
      # Add bootstrap distributions as density plots
      geom_density(data = all_bootstrap_data,
                   aes(x = bootstrap_depth/1e6, fill = technology, color = technology),
                   alpha = 0.3, linewidth = 0.8) +
      # Add vertical lines for point estimates
      geom_vline(data = point_estimates,
                 aes(xintercept = depth/1e6, color = technology),
                 linetype = "dashed", linewidth = 1.2, alpha = 0.8) +
      # Facet by threshold and protocol
      facet_grid(threshold_label ~ protocol, scales = "free", switch = "y") +
      # Color mapping
      scale_fill_manual(values = get_color_map(), name = "Technology") +
      scale_color_manual(values = get_color_map(), name = "Technology") +
      # Axis formatting
      scale_x_continuous(
        name = "Sequencing Depth (Million Reads)",
        labels = scales::comma
      ) +
      scale_y_continuous(
        name = "Bootstrap Density",
        labels = scales::number_format(accuracy = 0.01)
      ) +
      # Labels and theme
      labs(
        title = paste("Bootstrap Uncertainty Distributions:", feature_name),
        subtitle = "Density plots show bootstrap sampling distributions | Dashed lines show point estimates",
        caption = "Based on 1000 bootstrap resamples per threshold"
      ) +
      theme_minimal() +
      theme(
        text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray50"),
        plot.caption = element_text(size = 9, color = "gray60"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      ) +
      guides(
        fill = guide_legend(title = "Technology"),
        color = guide_legend(title = "Technology")
      )
    
    plot_list[[paste0(ft, "_BootstrapDistributions")]] <- p
  }
  
  return(plot_list)
}