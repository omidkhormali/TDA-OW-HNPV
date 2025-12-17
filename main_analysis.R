########################################################################
# MAIN ANALYSIS SCRIPT - SENSITIVITY COMPARISON
# Compares VAB, PL, PI (fixed) vs HNAV, HWNAV, OW-HNPV (varying n_sub)

########################################################################

# Clear workspace
rm(list = ls())
gc()

# Load libraries
library(igraph)
library(TDA)
library(tidyverse)
library(depthTools)
library(randomForest)
library(caret)
library(fda.usc)
library(pROC)
library(lubridate)
library(ggplot2)

options(dplyr.summarise.inform = FALSE)

########################################################################
# CONFIGURATION - CHANGE THIS TO YOUR PATH
########################################################################

baseDir <- "C:/Users/o_kho/OneDrive - University of Evansville/2025_fall/Stat_391/Paper_Ideas/R-Files/Cleaned"

# Parameters
topRank <- 250           # Number of top nodes to keep
threshold <- 0.05        # Price anomaly threshold
m_fixed <- 30            # Fixed number of main intervals
n_sub_seq <- c(1, 2, 3, 5, 10)  # Sequence of n_sub values to test
repNum <- 10             # Number of RF repetitions
filtration <- "sublevel"

# Fixed dimensions for VAB, PL, PI
fixed_dim <- 30

########################################################################
# SETUP DIRECTORIES
########################################################################

folders <- c('data', 'tokenPrice', 'pd', 'graph', 'depth', 'merge', 'model', 'results',
             'vab', 'pl', 'pi', 'hnav', 'hwnav', 'owhnpv')

for (f in folders) {
  dir_path <- file.path(baseDir, f)
  assign(paste0(f, 'Dir'), dir_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
    cat(paste("Created directory:", dir_path, "\n"))
  }
}

cat(paste("Base directory:", baseDir, "\n"))

# Load functions
source(file.path(baseDir, 'functions_clean.R'))

########################################################################
# LOAD DATA
########################################################################

cat("\n========================================\n")
cat("LOADING DATA\n")
cat("========================================\n\n")

allTokens <- readRDS(file.path(dataDir, 'allTokens.rds'))

# Select days with at least 5 tokens
selectedDays <- allTokens %>%
  group_by(time) %>%
  summarise(n = length(unique(name)), .groups = 'drop') %>%
  filter(n >= 5) %>%
  pull(time)

networkDF <- allTokens %>% filter(time %in% selectedDays)
networkDF$time <- ymd(networkDF$time)

# Normalize values
networkDF <- networkDF %>%
  group_by(name) %>%
  mutate(value = log(value) / log(max(value))) %>%
  ungroup() %>%
  dplyr::select(1:4)

cat(paste("Total days:", length(unique(networkDF$time)), "\n"))
cat(paste("Date range:", min(networkDF$time), "to", max(networkDF$time), "\n"))

########################################################################
# STEP 1: COMPUTE GRAPH FEATURES (ONCE)
########################################################################

if (!file.exists(file.path(graphDir, "graphFeature_allTokens.rds"))) {
  featureGraph(networkDF, graphDir)
} else {
  cat("Graph features already exist. Skipping...\n")
}

########################################################################
# STEP 2: COMPUTE PERSISTENCE DIAGRAMS (ONCE)
########################################################################

# Check if PDs already exist
pd_files <- list.files(pdDir, pattern = "^PD_")
if (length(pd_files) == 0) {
  computePD(networkDF, topRank, pdDir, filtration)
} else {
  cat(paste("Found", length(pd_files), "persistence diagrams. Skipping...\n"))
}

########################################################################
# STEP 3: COMPUTE FIXED METHODS (VAB, PL, PI) - ONCE
########################################################################

# VAB
cat("Computing VAB...\n")
computeVAB(networkDF, pdDir, vabDir, num_intervals = fixed_dim, filtration)
computeRollingDepth(vabDir, "vab", depthDir)

# PL
cat("Computing PL...\n")
computePL(networkDF, pdDir, plDir, num_intervals = fixed_dim, filtration)
computeRollingDepth(plDir, "pl", depthDir)

# PI
cat("Computing PI...\n")
computePI(networkDF, pdDir, piDir, res = 6, filtration)  # 6x6 = 36 ~ 30
computeRollingDepth(piDir, "pi", depthDir)

########################################################################
# STEP 4: FIT MODELS FOR FIXED METHODS (ONCE)
########################################################################

# Clear depth files for velocity methods 
old_velocity_depths <- list.files(depthDir, pattern = "^rd_(hnav|hwnav|owhnpv)", full.names = TRUE)
if (length(old_velocity_depths) > 0) {
  file.remove(old_velocity_depths)
}

# Merge data for fixed methods only
dataMerge(graphDir, depthDir, tokenPriceDir, mergeDir, threshold)

# Fit models for fixed methods
fixed_results <- list()

for (method in c("vab", "pl", "pi")) {
  cat(paste("\n--- Fitting RF for", toupper(method), "---\n"))
  fitRFModel(mergeDir, modelDir, method, threshold, repNum)
  summary <- summarizeResults(modelDir, method)
  gains <- computeAUCGain(summary)
  gains$method <- toupper(method)
  gains$n_sub <- NA  # Fixed methods don't depend on n_sub
  fixed_results[[method]] <- gains
}

fixed_df <- bind_rows(fixed_results)
cat("\n\nFixed methods results:\n")
print(fixed_df)

########################################################################
# STEP 5: SENSITIVITY ANALYSIS FOR VELOCITY METHODS
########################################################################

velocity_results <- list()

for (n_sub in n_sub_seq) {
  
  # Clear old velocity depth files
  old_files <- list.files(depthDir, pattern = "^rd_(hnav|hwnav|owhnpv)", full.names = TRUE)
  if (length(old_files) > 0) {
    file.remove(old_files)
    cat(paste("Removed", length(old_files), "old velocity depth files\n"))
  }
  
  # Compute HNAV
  cat("  Computing HNAV...\n")
  computeHNAV(networkDF, pdDir, hnavDir, m = m_fixed, n_sub = n_sub, filtration)
  computeRollingDepth(hnavDir, "hnav", depthDir)
  
  # Compute HWNAV
  cat("  Computing HWNAV...\n")
  computeHWNAV(networkDF, pdDir, hwnavDir, m = m_fixed, n_sub = n_sub, filtration)
  computeRollingDepth(hwnavDir, "hwnav", depthDir)
  
  # Compute OW-HNPV
  cat("  Computing OW-HNPV...\n")
  computeOWHNPV(networkDF, pdDir, owhnpvDir, m = m_fixed, n_sub = n_sub, filtration)
  computeRollingDepth(owhnpvDir, "owhnpv", depthDir)
  
  # Merge data (includes all methods now)
  cat("  Merging data...\n")
  dataMerge(graphDir, depthDir, tokenPriceDir, mergeDir, threshold)
  
  # Fit models for each velocity method
  for (method in c("hnav", "hwnav", "owhnpv")) {
    cat(paste("  Fitting RF for", toupper(method), "...\n"))
    fitRFModel(mergeDir, modelDir, method, threshold, repNum)
    summary <- summarizeResults(modelDir, method)
    gains <- computeAUCGain(summary)
    gains$method <- toupper(method)
    gains$n_sub <- n_sub
    velocity_results[[paste(method, n_sub, sep = "_")]] <- gains
  }
  
  cat(paste("\n  n_sub =", n_sub, "complete!\n"))
}

velocity_df <- bind_rows(velocity_results)

########################################################################
# STEP 6: COMBINE AND SAVE RESULTS
########################################################################

# Combine all results
all_results <- bind_rows(fixed_df, velocity_df)

# Save combined results
saveRDS(all_results, file.path(resultsDir, "sensitivity_all_results.rds"))
write.csv(all_results, file.path(resultsDir, "sensitivity_all_results.csv"), row.names = FALSE)

cat("Results saved!\n")
cat(paste("Total rows:", nrow(all_results), "\n"))
cat(paste("Methods:", paste(unique(all_results$method), collapse = ", "), "\n"))

########################################################################
# STEP 7: CREATE VISUALIZATIONS
########################################################################

# Prepare data for plotting
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

# Define colors
method_colors <- c(
  "VAB" = "#E91E63",
  "PL" = "#3498DB",
  "PI" = "#00BCD4",
  "HNAV" = "#E74C3C",
  "HWNAV" = "#F39C12",
  "OWHNPV" = "#27AE60"
)

########################################################################
# PLOT: COMBINED SENSITIVITY ANALYSIS
########################################################################

p_combined <- ggplot() +
  # Velocity methods as curves
  geom_line(data = velocity_plot,
            aes(x = n_sub, y = gain, color = method),
            linewidth = 1.2) +
  geom_point(data = velocity_plot,
             aes(x = n_sub, y = gain, color = method),
             size = 2.5) +
  # Fixed methods as horizontal lines
  geom_hline(data = fixed_plot,
             aes(yintercept = gain, color = method),
             linetype = "dashed", linewidth = 1, alpha = 0.8) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  # Facets
  facet_grid(model ~ horizon_num,
             labeller = labeller(
               horizon_num = function(x) paste0("h=", x),
               model = function(x) x
             )) +
  # Styling
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = n_sub_seq) +
  labs(
    title = sprintf("Sensitivity of Topological Summaries to n_sub (m=%d, topRank=%d)", 
                    m_fixed, topRank),
    subtitle = "Solid lines: velocity methods (vary with n_sub) | Dashed lines: fixed methods",
    x = expression(n[sub] ~ "(number of subintervals per main interval)"),
    y = "AUC Gain (%)",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  guides(color = guide_legend(nrow = 1))

# Save plot
ggsave(
  file.path(resultsDir, sprintf("sensitivity_combined_m%d_topRank%d.png", m_fixed, topRank)),
  p_combined,
  width = 16, height = 8, dpi = 300
)

cat("Combined plot saved!\n")

########################################################################
# PLOT: SEPARATE PLOTS FOR M2 AND M3
########################################################################

for (model_name in c("M2", "M3")) {
  
  fixed_model <- fixed_plot %>% filter(model == model_name)
  velocity_model <- velocity_plot %>% filter(model == model_name)
  
  p_model <- ggplot() +
    geom_line(data = velocity_model,
              aes(x = n_sub, y = gain, color = method),
              linewidth = 1.2) +
    geom_point(data = velocity_model,
               aes(x = n_sub, y = gain, color = method),
               size = 3) +
    geom_hline(data = fixed_model,
               aes(yintercept = gain, color = method),
               linetype = "dashed", linewidth = 1, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    facet_wrap(~ horizon_num, nrow = 2,
               labeller = labeller(horizon_num = function(x) paste0("Horizon ", x, " days"))) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = n_sub_seq) +
    labs(
      title = sprintf("Model %s: AUC Gain vs n_sub (m=%d, topRank=%d)", 
                      model_name, m_fixed, topRank),
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
    file.path(resultsDir, sprintf("sensitivity_%s_m%d_topRank%d.png", model_name, m_fixed, topRank)),
    p_model,
    width = 12, height = 8, dpi = 300
  )
  
  cat(paste("Plot for", model_name, "saved!\n"))
}

########################################################################
# SUMMARY TABLE
########################################################################

# Average gains across horizons
summary_table <- all_results %>%
  group_by(method, n_sub) %>%
  summarise(
    avg_M2_gain = mean(M2_gain, na.rm = TRUE),
    avg_M3_gain = mean(M3_gain, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(avg_M3_gain))

print(summary_table)

write.csv(summary_table, file.path(resultsDir, "summary_table.csv"), row.names = FALSE)

########################################################################
# BEST CONFIGURATION REPORT
########################################################################

# Best velocity method configuration
best_velocity <- velocity_df %>%
  group_by(method, n_sub) %>%
  summarise(avg_M3_gain = mean(M3_gain, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(avg_M3_gain)) %>%
  head(5)

cat("Top 5 velocity method configurations (by avg M3 gain):\n")
print(best_velocity)

# Best overall
best_overall <- all_results %>%
  group_by(method, n_sub) %>%
  summarise(avg_M3_gain = mean(M3_gain, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(avg_M3_gain)) %>%
  head(1)

cat("\nBest overall configuration:\n")
cat(paste("  Method:", best_overall$method, "\n"))
cat(paste("  n_sub:", ifelse(is.na(best_overall$n_sub), "N/A (fixed)", best_overall$n_sub), "\n"))
cat(paste("  Avg M3 AUC Gain:", round(best_overall$avg_M3_gain, 2), "%\n"))

cat("ANALYSIS COMPLETE!\n")
cat(paste("Results saved to:", resultsDir, "\n"))
