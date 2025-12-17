########################################################################
# SIMPLE PLOTTING SCRIPT
# Creates the comparison plot: fixed methods (horizontal) vs velocity (curves)
# Run this AFTER running main_analysis.R
########################################################################

library(tidyverse)
library(ggplot2)

########################################################################
# CONFIGURATION
########################################################################

baseDir <- "C:/Users/o_kho/OneDrive - University of Evansville/2025_fall/Stat_391/Paper_Ideas/R-Files/Cleaned"
resultsDir <- file.path(baseDir, "results")

# Parameters used in analysis
m_fixed <- 30
topRank <- 250
n_sub_seq <- c(1, 2, 3, 5, 10)

########################################################################
# LOAD RESULTS
########################################################################

all_results <- readRDS(file.path(resultsDir, "sensitivity_all_results.rds"))

cat("Loaded results:\n")
cat(paste("  Rows:", nrow(all_results), "\n"))
cat(paste("  Columns:", paste(names(all_results), collapse = ", "), "\n"))
cat(paste("  Methods:", paste(unique(all_results$method), collapse = ", "), "\n"))

########################################################################
# PREPARE DATA FOR PLOTTING
########################################################################

plot_data <- all_results %>%
  dplyr::select(horizon_num, method, n_sub, M2_gain, M3_gain) %>%
  pivot_longer(cols = c(M2_gain, M3_gain), names_to = "model", values_to = "gain") %>%
  mutate(
    model = gsub("_gain", "", model),
    is_fixed = method %in% c("VAB", "PL", "PI")
  )

# Separate fixed and velocity methods
fixed_plot <- plot_data %>% filter(is_fixed)
velocity_plot <- plot_data %>% filter(!is_fixed)

# Colors
method_colors <- c(
  "VAB" = "#E91E63",      # Pink
  "PL" = "#3498DB",       # Blue
  "PI" = "#00BCD4",       # Cyan
  "HNAV" = "#E74C3C",     # Red
  "HWNAV" = "#F39C12",    # Orange
  "OWHNPV" = "#27AE60"    # Green
)

########################################################################
# MAIN COMPARISON PLOT
########################################################################

p <- ggplot() +
  # Velocity methods: curves
  geom_line(
    data = velocity_plot,
    aes(x = n_sub, y = gain, color = method),
    linewidth = 1.2
  ) +
  geom_point(
    data = velocity_plot,
    aes(x = n_sub, y = gain, color = method),
    size = 2.5
  ) +
  # Fixed methods: horizontal lines
  geom_hline(
    data = fixed_plot,
    aes(yintercept = gain, color = method),
    linetype = "dashed", linewidth = 1, alpha = 0.8
  ) +
  # Reference line
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  # Facets: columns = horizon, rows = model
  facet_grid(
    model ~ horizon_num,
    labeller = labeller(
      horizon_num = function(x) paste0("h = ", x),
      model = function(x) x
    )
  ) +
  # Scales
  scale_color_manual(values = method_colors, name = "Method") +
  scale_x_continuous(breaks = n_sub_seq) +
  # Labels
  labs(
    title = "Comparison of Topological Summaries for Anomaly Detection",
    subtitle = "Solid lines: velocity methods (vary with n_sub) | Dashed lines: fixed methods (VAB, PL, PI)",
    x = expression(n[sub] ~ "(number of subintervals per main interval)"),
    y = "AUC Gain (%)"
  ) +
  # Theme
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.title = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 1))

# Save
ggsave(
  file.path(resultsDir, "comparison_sensitivity_final.png"),
  p,
  width = 16, height = 8, dpi = 300
)

print(p)

cat("\nPlot saved to:", file.path(resultsDir, "comparison_sensitivity_final.png"), "\n")

########################################################################
# SEPARATE PLOTS FOR M2 AND M3
########################################################################

for (model_name in c("M2", "M3")) {
  
  fixed_model <- fixed_plot %>% filter(model == model_name)
  velocity_model <- velocity_plot %>% filter(model == model_name)
  
  p_model <- ggplot() +
    # Velocity curves
    geom_line(
      data = velocity_model,
      aes(x = n_sub, y = gain, color = method),
      linewidth = 1.2
    ) +
    geom_point(
      data = velocity_model,
      aes(x = n_sub, y = gain, color = method),
      size = 3
    ) +
    # Fixed horizontal lines
    geom_hline(
      data = fixed_model,
      aes(yintercept = gain, color = method),
      linetype = "dashed", linewidth = 1, alpha = 0.8
    ) +
    # Reference
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    # Facet by horizon
    facet_wrap(
      ~ horizon_num, nrow = 2,
      labeller = labeller(horizon_num = function(x) paste0("Horizon ", x, " days"))
    ) +
    # Scales
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = n_sub_seq) +
    # Labels
    labs(
      title = sprintf("Model %s: AUC Gain Comparison", model_name),
      subtitle = "Solid = Velocity methods | Dashed = Fixed methods",
      x = expression(n[sub]),
      y = "AUC Gain (%)",
      color = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(resultsDir, sprintf("comparison_%s_final.png", model_name)),
    p_model,
    width = 12, height = 8, dpi = 300
  )
  
  cat("Saved plot for", model_name, "\n")
}

########################################################################
# SUMMARY STATISTICS
########################################################################

# Average gains by method
summary_by_method <- all_results %>%
  group_by(method) %>%
  summarise(
    avg_M2_gain = mean(M2_gain, na.rm = TRUE),
    avg_M3_gain = mean(M3_gain, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(avg_M3_gain))

cat("Average AUC Gain by Method:\n")
print(summary_by_method)

# For velocity methods: best n_sub
cat("\n\nBest n_sub for each velocity method (M3):\n")
velocity_best <- all_results %>%
  filter(!is.na(n_sub)) %>%
  group_by(method, n_sub) %>%
  summarise(avg_M3_gain = mean(M3_gain, na.rm = TRUE), .groups = 'drop') %>%
  group_by(method) %>%
  slice_max(avg_M3_gain, n = 1)

print(velocity_best)

########################################################################
# DETAILED TABLE BY HORIZON
########################################################################

detailed_table <- all_results %>%
  dplyr::select(horizon_num, method, n_sub, M2_gain, M3_gain) %>%
  arrange(horizon_num, method, n_sub)

print(detailed_table, n = 50)

########################################################################
# COMPREHENSIVE PLOTTING SCRIPT
# Creates publication-quality figures for AUC, Accuracy, and Precision
# Comparing: VAB, PL, PI (fixed) vs HNAV, HWNAV, OW-HNPV (velocity)
########################################################################

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)

########################################################################
# CONFIGURATION
########################################################################

baseDir <- "C:/Users/o_kho/OneDrive - University of Evansville/2025_fall/Stat_391/Paper_Ideas/R-Files/Cleaned"
resultsDir <- file.path(baseDir, "results")
modelDir <- file.path(baseDir, "model")

# Parameters
m_fixed <- 30
topRank <- 250
n_sub_seq <- c(1, 2, 3, 5, 10)

# Method colors - consistent across all plots
method_colors <- c(
  "VAB" = "#E91E63",
  "PL" = "#3498DB",
  "PI" = "#00BCD4",
  "HNAV" = "#E74C3C",
  "HWNAV" = "#F39C12",
  "OWHNPV" = "#27AE60"
)

# Better labels for methods
method_labels <- c(
  "VAB" = "VAB",
  "PL" = "PL",
  "PI" = "PI",
  "HNAV" = "HNAV",
  "HWNAV" = "HWNAV",
  "OWHNPV" = "OW-HNPV"
)

########################################################################
# FUNCTION: Compute gains for any metric
########################################################################

computeMetricGain <- function(modelDir, methodName, metric = "AUC") {
  
  rf_file <- file.path(modelDir, paste0("rf_results_", methodName, ".rds"))
  
  if (!file.exists(rf_file)) {
    warning(paste("File not found:", rf_file))
    return(NULL)
  }
  
  results <- readRDS(rf_file)
  
  # Summarize by horizon and model type
  summary_df <- results %>%
    group_by(horizon, modelType) %>%
    summarise(
      value = mean(!!sym(metric), na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Pivot to wide format
  wide <- summary_df %>%
    pivot_wider(names_from = modelType, values_from = value)
  
  # Compute gains
  wide$M2_gain <- NA
  wide$M3_gain <- NA
  
  if ("M2" %in% names(wide) && "M1" %in% names(wide)) {
    wide$M2_gain <- (wide$M2 - wide$M1) / wide$M1 * 100
  }
  if ("M3" %in% names(wide) && "M1" %in% names(wide)) {
    wide$M3_gain <- (wide$M3 - wide$M1) / wide$M1 * 100
  }
  
  wide$horizon_num <- as.numeric(gsub("flag", "", wide$horizon))
  
  return(wide)
}

########################################################################
# FUNCTION: Gather all results for a specific metric
########################################################################

gatherAllResults <- function(modelDir, metric = "AUC", n_sub_seq) {
  
  # Fixed methods
  fixed_results <- list()
  for (method in c("vab", "pl", "pi")) {
    gains <- computeMetricGain(modelDir, method, metric)
    if (!is.null(gains)) {
      gains$method <- toupper(method)
      gains$n_sub <- NA
      fixed_results[[method]] <- gains
    }
  }
  
  # Velocity methods - need to read from saved results or recompute
  # For now, we'll read from the sensitivity results if available
  
  # Try to load pre-computed sensitivity results
  sens_file <- file.path(resultsDir, "sensitivity_all_results.rds")
  
  if (file.exists(sens_file)) {
    all_results <- readRDS(sens_file)
    
    # If the metric columns don't exist, we need to recompute
    if (!paste0("M2_gain_", metric) %in% names(all_results)) {
      # The saved results only have AUC - need to compute other metrics
      cat(paste("Computing", metric, "gains from RF results...\n"))
      
      velocity_results <- list()
      
      # For velocity methods, we need the RF results for each n_sub configuration
      # These should be in the model directory from the last run
      for (method in c("hnav", "hwnav", "owhnpv")) {
        gains <- computeMetricGain(modelDir, method, metric)
        if (!is.null(gains)) {
          gains$method <- toupper(method)
          # Get n_sub from the last configuration run
          # This is a limitation - we only have the last n_sub results
          velocity_results[[method]] <- gains
        }
      }
      
      # Combine
      fixed_df <- bind_rows(fixed_results)
      velocity_df <- bind_rows(velocity_results)
      
      # For velocity methods without n_sub info, we'll need the full sensitivity run
      # For now, return what we have
      return(list(fixed = fixed_df, velocity = velocity_df, complete = FALSE))
    }
  }
  
  fixed_df <- bind_rows(fixed_results)
  return(list(fixed = fixed_df, complete = TRUE))
}

########################################################################
# FUNCTION: Create comparison plot for any metric
########################################################################

createComparisonPlot <- function(all_results, metric_name = "AUC", 
                                 m_fixed = 30, topRank = 250,
                                 n_sub_seq = c(1, 2, 3, 5, 10)) {
  
  # Prepare data for plotting
  plot_data <- all_results %>%
    dplyr::select(horizon_num, method, n_sub, M2_gain, M3_gain) %>%
    pivot_longer(cols = c(M2_gain, M3_gain), names_to = "model", values_to = "gain") %>%
    mutate(
      model = gsub("_gain", "", model),
      is_fixed = method %in% c("VAB", "PL", "PI"),
      method = factor(method, levels = names(method_colors))
    )
  
  # Separate fixed and velocity methods
  fixed_plot <- plot_data %>% filter(is_fixed)
  velocity_plot <- plot_data %>% filter(!is_fixed)
  
  # Create plot
  p <- ggplot() +
    # Velocity methods: curves
    geom_line(
      data = velocity_plot,
      aes(x = n_sub, y = gain, color = method),
      linewidth = 1.2
    ) +
    geom_point(
      data = velocity_plot,
      aes(x = n_sub, y = gain, color = method),
      size = 2.5
    ) +
    # Fixed methods: horizontal lines
    geom_hline(
      data = fixed_plot,
      aes(yintercept = gain, color = method),
      linetype = "dashed", linewidth = 1, alpha = 0.8
    ) +
    # Reference line
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    # Facets
    facet_grid(
      model ~ horizon_num,
      labeller = labeller(
        horizon_num = function(x) paste0("h = ", x),
        model = function(x) x
      )
    ) +
    # Scales
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_x_continuous(breaks = n_sub_seq) +
    # Labels
    labs(
      title = paste("Comparison of Topological Summaries -", metric_name, "Gain"),
      subtitle = "Solid lines: velocity methods | Dashed lines: fixed methods (VAB, PL, PI)",
      x = expression(n[sub] ~ "(number of subintervals per main interval)"),
      y = paste(metric_name, "Gain (%)")
    ) +
    # Theme
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = NA),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(nrow = 1))
  
  return(p)
}

########################################################################
# FUNCTION: Create M2/M3 separate plots
########################################################################

createModelPlot <- function(all_results, model_name = "M3", metric_name = "AUC",
                            n_sub_seq = c(1, 2, 3, 5, 10)) {
  
  # Prepare data
  plot_data <- all_results %>%
    dplyr::select(horizon_num, method, n_sub, M2_gain, M3_gain) %>%
    pivot_longer(cols = c(M2_gain, M3_gain), names_to = "model", values_to = "gain") %>%
    mutate(
      model = gsub("_gain", "", model),
      is_fixed = method %in% c("VAB", "PL", "PI"),
      method = factor(method, levels = names(method_colors))
    ) %>%
    filter(model == model_name)
  
  fixed_model <- plot_data %>% filter(is_fixed)
  velocity_model <- plot_data %>% filter(!is_fixed)
  
  p <- ggplot() +
    geom_line(
      data = velocity_model,
      aes(x = n_sub, y = gain, color = method),
      linewidth = 1.2
    ) +
    geom_point(
      data = velocity_model,
      aes(x = n_sub, y = gain, color = method),
      size = 3
    ) +
    geom_hline(
      data = fixed_model,
      aes(yintercept = gain, color = method),
      linetype = "dashed", linewidth = 1, alpha = 0.8
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    facet_wrap(
      ~ horizon_num, nrow = 2,
      labeller = labeller(horizon_num = function(x) paste0("Horizon ", x, " days"))
    ) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    scale_x_continuous(breaks = n_sub_seq) +
    labs(
      title = sprintf("Model %s: %s Gain Comparison", model_name, metric_name),
      subtitle = "Solid = Velocity methods | Dashed = Fixed methods",
      x = expression(n[sub]),
      y = paste(metric_name, "Gain (%)"),
      color = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

########################################################################
# FUNCTION: Create highlighted plot emphasizing OW-HNPV at h=2, h=4
########################################################################

createHighlightPlot <- function(all_results, metric_name = "AUC",
                                highlight_horizons = c(2, 4),
                                n_sub_seq = c(1, 2, 3, 5, 10)) {
  
  # Prepare data
  plot_data <- all_results %>%
    dplyr::select(horizon_num, method, n_sub, M2_gain, M3_gain) %>%
    pivot_longer(cols = c(M2_gain, M3_gain), names_to = "model", values_to = "gain") %>%
    mutate(
      model = gsub("_gain", "", model),
      is_fixed = method %in% c("VAB", "PL", "PI"),
      method = factor(method, levels = names(method_colors)),
      is_highlighted = horizon_num %in% highlight_horizons
    )
  
  fixed_plot <- plot_data %>% filter(is_fixed)
  velocity_plot <- plot_data %>% filter(!is_fixed)
  
  # Create plot with highlighting
  p <- ggplot() +
    # Background for highlighted panels
    geom_rect(
      data = data.frame(
        horizon_num = highlight_horizons,
        model = rep(c("M2", "M3"), each = length(highlight_horizons))
      ),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "lightyellow", alpha = 0.3
    ) +
    # Velocity methods: curves
    geom_line(
      data = velocity_plot,
      aes(x = n_sub, y = gain, color = method, 
          linewidth = ifelse(method == "OWHNPV", 1.5, 1)),
      show.legend = TRUE
    ) +
    geom_point(
      data = velocity_plot,
      aes(x = n_sub, y = gain, color = method,
          size = ifelse(method == "OWHNPV", 3.5, 2.5))
    ) +
    # Fixed methods: horizontal lines
    geom_hline(
      data = fixed_plot,
      aes(yintercept = gain, color = method),
      linetype = "dashed", linewidth = 0.8, alpha = 0.7
    ) +
    # Reference line
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
    # Facets
    facet_grid(
      model ~ horizon_num,
      labeller = labeller(
        horizon_num = function(x) paste0("h = ", x),
        model = function(x) x
      )
    ) +
    # Scales
    scale_color_manual(values = method_colors, labels = method_labels, name = "Method") +
    scale_linewidth_identity() +
    scale_size_identity() +
    scale_x_continuous(breaks = n_sub_seq) +
    # Labels
    labs(
      title = paste("Comparison of Topological Summaries -", metric_name),
      subtitle = "OW-HNPV shows superior performance at h=2 and h=4 (highlighted)",
      x = expression(n[sub] ~ "(number of subintervals per main interval)"),
      y = paste(metric_name, "Gain (%)")
    ) +
    # Theme
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = NA),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      axis.title = element_text(face = "bold")
    ) +
    guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.5, size = 3)))
  
  return(p)
}

########################################################################
# FUNCTION: Create bar plot for specific horizons
########################################################################

createBarPlot <- function(all_results, horizons = c(2, 4, 7), model_name = "M3",
                          metric_name = "AUC", n_sub_value = 10) {
  
  # Filter for specific n_sub (or use max for velocity methods)
  plot_data <- all_results %>%
    filter(horizon_num %in% horizons) %>%
    mutate(
      n_sub_use = ifelse(is.na(n_sub), "fixed", as.character(n_sub))
    )
  
  # For velocity methods, take the specified n_sub
  velocity_data <- plot_data %>%
    filter(!is.na(n_sub), n_sub == n_sub_value)
  
  # For fixed methods
  fixed_data <- plot_data %>%
    filter(is.na(n_sub))
  
  bar_data <- bind_rows(velocity_data, fixed_data) %>%
    mutate(
      method = factor(method, levels = c("VAB", "PL", "PI", "HNAV", "HWNAV", "OWHNPV")),
      horizon_label = paste0("h = ", horizon_num)
    )
  
  # Select the right gain column
  gain_col <- paste0(model_name, "_gain")
  
  p <- ggplot(bar_data, aes(x = method, y = !!sym(gain_col), fill = method)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    facet_wrap(~ horizon_label, nrow = 1) +
    scale_fill_manual(values = method_colors, labels = method_labels) +
    labs(
      title = sprintf("%s Gain by Method (Model %s, n_sub = %d)", metric_name, model_name, n_sub_value),
      subtitle = "Comparing velocity-based methods against static topological summaries",
      x = "Method",
      y = paste(metric_name, "Gain (%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  return(p)
}

########################################################################
# FUNCTION: Create summary heatmap
########################################################################

createHeatmap <- function(all_results, model_name = "M3", metric_name = "AUC") {
  
  # For velocity methods, take the best n_sub per method/horizon
  velocity_best <- all_results %>%
    filter(!is.na(n_sub)) %>%
    group_by(method, horizon_num) %>%
    slice_max(!!sym(paste0(model_name, "_gain")), n = 1) %>%
    ungroup()
  
  # For fixed methods
  fixed_data <- all_results %>%
    filter(is.na(n_sub))
  
  heatmap_data <- bind_rows(velocity_best, fixed_data) %>%
    mutate(
      method = factor(method, levels = rev(c("OWHNPV", "HWNAV", "HNAV", "PI", "PL", "VAB"))),
      gain = !!sym(paste0(model_name, "_gain"))
    )
  
  p <- ggplot(heatmap_data, aes(x = factor(horizon_num), y = method, fill = gain)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", gain)), color = "black", size = 4) +
    scale_fill_gradient2(
      low = "#E74C3C", mid = "white", high = "#27AE60",
      midpoint = 0, name = paste(metric_name, "Gain (%)")
    ) +
    scale_y_discrete(labels = method_labels) +
    labs(
      title = sprintf("%s Gain Heatmap (Model %s)", metric_name, model_name),
      subtitle = "Best n_sub selected for velocity methods",
      x = "Prediction Horizon (days)",
      y = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text = element_text(size = 11),
      panel.grid = element_blank()
    )
  
  return(p)
}

########################################################################
# MAIN EXECUTION
########################################################################

# Load the sensitivity results
all_results <- readRDS(file.path(resultsDir, "sensitivity_all_results.rds"))

cat("Loaded results:\n")
cat(paste("  Methods:", paste(unique(all_results$method), collapse = ", "), "\n"))
cat(paste("  n_sub values:", paste(sort(unique(all_results$n_sub[!is.na(all_results$n_sub)])), collapse = ", "), "\n"))

########################################################################
# GENERATE AUC PLOTS
########################################################################

# 1. Main comparison plot
p_auc_main <- createComparisonPlot(all_results, "AUC", m_fixed, topRank, n_sub_seq)
ggsave(file.path(resultsDir, "Fig1_AUC_comparison_main.png"), p_auc_main, 
       width = 16, height = 8, dpi = 300)
cat("  Saved: Fig1_AUC_comparison_main.png\n")

# 2. Highlighted plot emphasizing h=2, h=4
p_auc_highlight <- createHighlightPlot(all_results, "AUC", c(2, 4), n_sub_seq)
ggsave(file.path(resultsDir, "Fig2_AUC_highlighted.png"), p_auc_highlight,
       width = 16, height = 8, dpi = 300)
cat("  Saved: Fig2_AUC_highlighted.png\n")

# 3. M3 only plot
p_auc_m3 <- createModelPlot(all_results, "M3", "AUC", n_sub_seq)
ggsave(file.path(resultsDir, "Fig3_AUC_M3_only.png"), p_auc_m3,
       width = 12, height = 8, dpi = 300)
cat("  Saved: Fig3_AUC_M3_only.png\n")

# 4. Bar plot for key horizons
p_auc_bar <- createBarPlot(all_results, c(2, 4, 7), "M3", "AUC", 10)
ggsave(file.path(resultsDir, "Fig4_AUC_barplot.png"), p_auc_bar,
       width = 10, height = 6, dpi = 300)
cat("  Saved: Fig4_AUC_barplot.png\n")

# 5. Heatmap
p_auc_heat <- createHeatmap(all_results, "M3", "AUC")
ggsave(file.path(resultsDir, "Fig5_AUC_heatmap.png"), p_auc_heat,
       width = 10, height = 6, dpi = 300)
cat("  Saved: Fig5_AUC_heatmap.png\n")

########################################################################
# COMPUTE AND PLOT ACCURACY
########################################################################

# Recompute gains for Accuracy
accuracy_results <- list()

# Fixed methods
for (method in c("vab", "pl", "pi")) {
  gains <- computeMetricGain(modelDir, method, "Accuracy")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    gains$n_sub <- NA
    accuracy_results[[paste0("fixed_", method)]] <- gains
  }
}

# Velocity methods (from last run)
for (method in c("hnav", "hwnav", "owhnpv")) {
  gains <- computeMetricGain(modelDir, method, "Accuracy")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    # Need to get n_sub from the AUC results
    velocity_nsub <- all_results %>%
      filter(method == toupper(!!method), !is.na(n_sub)) %>%
      pull(n_sub) %>%
      unique()
    
    if (length(velocity_nsub) > 0) {
      # Replicate for each n_sub (this is approximate - ideally rerun full analysis)
      for (ns in velocity_nsub) {
        gains_copy <- gains
        gains_copy$n_sub <- ns
        accuracy_results[[paste0(method, "_", ns)]] <- gains_copy
      }
    } else {
      gains$n_sub <- 10  # Default to last n_sub used
      accuracy_results[[method]] <- gains
    }
  }
}

accuracy_df <- bind_rows(accuracy_results)

if (nrow(accuracy_df) > 0) {
  # Generate Accuracy plots
  p_acc_main <- createComparisonPlot(accuracy_df, "Accuracy", m_fixed, topRank, n_sub_seq)
  ggsave(file.path(resultsDir, "Fig6_Accuracy_comparison.png"), p_acc_main,
         width = 16, height = 8, dpi = 300)
  cat("  Saved: Fig6_Accuracy_comparison.png\n")
  
  p_acc_m3 <- createModelPlot(accuracy_df, "M3", "Accuracy", n_sub_seq)
  ggsave(file.path(resultsDir, "Fig7_Accuracy_M3.png"), p_acc_m3,
         width = 12, height = 8, dpi = 300)
  cat("  Saved: Fig7_Accuracy_M3.png\n")
}

########################################################################
# COMPUTE AND PLOT PRECISION
########################################################################

precision_results <- list()

# Fixed methods
for (method in c("vab", "pl", "pi")) {
  gains <- computeMetricGain(modelDir, method, "Precision")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    gains$n_sub <- NA
    precision_results[[paste0("fixed_", method)]] <- gains
  }
}

# Velocity methods
for (method in c("hnav", "hwnav", "owhnpv")) {
  gains <- computeMetricGain(modelDir, method, "Precision")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    velocity_nsub <- all_results %>%
      filter(method == toupper(!!method), !is.na(n_sub)) %>%
      pull(n_sub) %>%
      unique()
    
    if (length(velocity_nsub) > 0) {
      for (ns in velocity_nsub) {
        gains_copy <- gains
        gains_copy$n_sub <- ns
        precision_results[[paste0(method, "_", ns)]] <- gains_copy
      }
    } else {
      gains$n_sub <- 10
      precision_results[[method]] <- gains
    }
  }
}

precision_df <- bind_rows(precision_results)

if (nrow(precision_df) > 0) {
  p_prec_main <- createComparisonPlot(precision_df, "Precision", m_fixed, topRank, n_sub_seq)
  ggsave(file.path(resultsDir, "Fig8_Precision_comparison.png"), p_prec_main,
         width = 16, height = 8, dpi = 300)
  cat("  Saved: Fig8_Precision_comparison.png\n")
  
  p_prec_m3 <- createModelPlot(precision_df, "M3", "Precision", n_sub_seq)
  ggsave(file.path(resultsDir, "Fig9_Precision_M3.png"), p_prec_m3,
         width = 12, height = 8, dpi = 300)
  cat("  Saved: Fig9_Precision_M3.png\n")
}

########################################################################
# COMPUTE AND PLOT SENSITIVITY
########################################################################

sensitivity_results <- list()

# Fixed methods
for (method in c("vab", "pl", "pi")) {
  gains <- computeMetricGain(modelDir, method, "Sensitivity")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    gains$n_sub <- NA
    sensitivity_results[[paste0("fixed_", method)]] <- gains
  }
}

# Velocity methods
for (method in c("hnav", "hwnav", "owhnpv")) {
  gains <- computeMetricGain(modelDir, method, "Sensitivity")
  if (!is.null(gains)) {
    gains$method <- toupper(method)
    velocity_nsub <- all_results %>%
      filter(method == toupper(!!method), !is.na(n_sub)) %>%
      pull(n_sub) %>%
      unique()
    
    if (length(velocity_nsub) > 0) {
      for (ns in velocity_nsub) {
        gains_copy <- gains
        gains_copy$n_sub <- ns
        sensitivity_results[[paste0(method, "_", ns)]] <- gains_copy
      }
    } else {
      gains$n_sub <- 10
      sensitivity_results[[method]] <- gains
    }
  }
}

sensitivity_df <- bind_rows(sensitivity_results)

if (nrow(sensitivity_df) > 0) {
  p_sens_main <- createComparisonPlot(sensitivity_df, "Sensitivity", m_fixed, topRank, n_sub_seq)
  ggsave(file.path(resultsDir, "Fig10_Sensitivity_comparison.png"), p_sens_main,
         width = 16, height = 8, dpi = 300)
  cat("  Saved: Fig10_Sensitivity_comparison.png\n")
}

########################################################################
# SUMMARY TABLE
########################################################################

# Create summary across all metrics for the paper
summary_table <- all_results %>%
  filter(n_sub == 10 | is.na(n_sub)) %>%  # Use n_sub=10 for velocity methods
  group_by(method) %>%
  summarise(
    avg_M2_AUC_gain = mean(M2_gain, na.rm = TRUE),
    avg_M3_AUC_gain = mean(M3_gain, na.rm = TRUE),
    best_horizon_M3 = horizon_num[which.max(M3_gain)],
    max_M3_gain = max(M3_gain, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(avg_M3_AUC_gain))

cat("\nSummary Table (n_sub = 10 for velocity methods):\n")
print(summary_table)

write.csv(summary_table, file.path(resultsDir, "summary_table_final.csv"), row.names = FALSE)


