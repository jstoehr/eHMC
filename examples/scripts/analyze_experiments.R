rm(list = ls())

# --------------------------------------------------
# --- LIBRARIES AND OPTIONS
# --------------------------------------------------

library(tidyverse)
library(ggplot2)
library(scales)

repo_root <- normalizePath(".")
fig_dir <- file.path(repo_root, "examples", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# --- SOURCE R FUNCTIONS AND MODELS
# --------------------------------------------------

source(file.path(repo_root, "examples", "R", "io_helpers.R"))
source(file.path(repo_root, "examples", "R", "analysis_helpers.R"))
source(file.path(repo_root, "examples", "R", "plot_helpers.R"))

# --------------------------------------------------
# --- COLLECTING RESULTS
# --------------------------------------------------

config_grid <- build_experiment_grid(
  models = c("banana", "MVNorm", "BLP"),
  seeds = 1:20,
  chains = 50,
  deltas = c(0.651, 0.8),
  n_train = 50000,
  refresh = 0.75,
  rand = c(TRUE, FALSE)
)

summaries <- collect_summaries(config_grid, repo_root)

param_stats <- summaries$param_stats %>%
  add_analysis_labels()

diagnosis <- summaries$diagnosis %>%
  add_analysis_labels()

# --------------------------------------------------
# FIGURES
# --------------------------------------------------

test_list <- list(
  c1 = list(algo_list = c("nuts", "ehmc"), rand = NULL, nrow = 1),
  c2 = list(algo_list = c("nuts", "ehmc", "prhmc"), rand = TRUE, nrow = 2),
  c3 = list(algo_list = c("nuts", "ehmc", "prhmc"), rand = FALSE, nrow = 2)
)

for (test in test_list) {
  prefix <- paste(test$algo_list, collapse = "_")
  
  if (!is.null(test$rand) && "prhmc" %in% test$algo_list) {
    prefix <- paste0(prefix, "_", test$rand)
  }
  
  # --- ESJD
  plot_boxplot(
    diagnosis %>% filter(
      algo %in% test$algo_list, 
      if (is.null(test$rand)) {
        TRUE
      } else {
        algo != "prhmc" | rand_leap == test$rand
      }
    ), 
    metric = "ESJD_norm", 
    legend_nrow = test$nrow
  )
  
  ggsave(
    filename = file.path(
      fig_dir, paste(prefix, "esjd_norm.jpeg", sep = "_")
    ), 
    width = 14, height = 14, unit = "cm", dpi = 300
  )
  
  # --- ESS
  tab_ess <- param_stats %>% filter(
    algo %in% test$algo_list, 
    par != "lp__",
    if (is.null(test$rand)) {
      TRUE
    } else {
      algo != "prhmc" | rand_leap == test$rand
    }
  )
  
  for (m in c("n_eff_norm", "n_eff_var_norm")) {
    plot_boxplot(
      tab_ess, 
      metric = m, 
      legend_nrow = test$nrow
    )
    
    ggsave(
      filename = file.path(
        fig_dir, paste0(prefix, "_", m, ".jpeg")
      ), 
      width = 14, height = 14, unit = "cm", dpi = 300
    )
  }
  
  # --- MIN ESS
  min_ess <- tab_ess %>%
    group_by(model, method_code, delta, seed, x_group, delta_label, algo) %>%
    summarise(
      min_n_eff_norm = min(n_eff_norm, na.rm = TRUE),
      min_n_eff_var_norm = min(n_eff_var_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
  for (m in c("min_n_eff_norm", "min_n_eff_var_norm")) {
    plot_boxplot(
      min_ess,
      metric = m,
      legend_nrow = test$nrow
    )
    ggsave(
      filename = file.path(
        fig_dir, paste0(prefix, "_", m, ".jpeg")
      ),
      width = 14, height = 14, unit = "cm", dpi = 300
    )
  }
}

# --------------------------------------------------
# SUMMARY TABLES
# --------------------------------------------------

# --- TABLES ON ESJD
summary_esjd <- diagnosis %>%
  filter(algo %in% c("nuts", "ehmc")) %>%
  group_by(model, method_code, delta, x_group, delta_label, algo) %>%
  summarise(
    mean_esjd = mean(ESJD_norm, na.rm = TRUE),
    sd_esjd = sd(ESJD_norm, na.rm = TRUE),
    .groups = "drop"
  )

# --- TABLES ON MIN ESS
min_ess <- param_stats %>%
  filter(algo %in% c("nuts", "ehmc"), par != "lp__") %>%
  group_by(model, method_code, delta, seed, x_group, delta_label, algo) %>%
  summarise(
    min_n_eff_norm = min(n_eff_norm, na.rm = TRUE),
    min_n_eff_var_norm = min(n_eff_var_norm, na.rm = TRUE),
    .groups = "drop"
  )

summary_min_ess <- min_ess %>%
  group_by(model, method_code, delta, x_group, delta_label, algo) %>%
  summarise(
    mean_n_eff = mean(min_n_eff_norm, na.rm = TRUE),
    median_n_eff = median(min_n_eff_norm, na.rm = TRUE),
    sd_n_eff = sd(min_n_eff_norm, na.rm = TRUE),
    mean_n_eff_var = mean(min_n_eff_var_norm, na.rm = TRUE),
    median_n_eff_var = median(min_n_eff_var_norm, na.rm = TRUE),
    sd_n_eff_var = sd(min_n_eff_var_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(
    model, 
    method_code, 
    delta, 
    mean_n_eff, sd_n_eff, 
    mean_n_eff_var, sd_n_eff_var
  )

# --------------------------------------------------
# SUBVIEWS
# --------------------------------------------------

subview <- function(df, delta_value, mean_col, sd_col, scale) {
  mean_col <- rlang::as_name(rlang::ensym(mean_col))
  sd_col <- rlang::as_name(rlang::ensym(sd_col))
  
  df %>%
    filter(delta == delta_value) %>%
    transmute(
      model,
      method_code,
      mean = signif(scale * .data[[mean_col]], 3),
      sd = signif(2 * scale * .data[[sd_col]], 2)
    )
}

subview(summary_esjd, 0.651, mean_esjd, sd_esjd, 1e6)
subview(summary_esjd, 0.8, mean_esjd, sd_esjd, 1e6)

subview(summary_min_ess, 0.651, mean_n_eff, sd_n_eff, 1e4)
subview(summary_min_ess, 0.8, mean_n_eff, sd_n_eff, 1e4)

subview(summary_min_ess, 0.651, mean_n_eff_var, sd_n_eff_var, 1e4)
subview(summary_min_ess, 0.8, mean_n_eff_var, sd_n_eff_var, 1e4)
