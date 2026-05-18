build_experiment_grid <- function(
    models = c("banana", "MVNorm", "BLP"),
    seeds = 1:20,
    chains = 50,
    warmup_nuts = 2000,
    warmup_ehmc = c(2000, 5000),
    iter = 2000,
    deltas = c(0.651, 0.8),
    metric = "diag",
    n_train = 50000,
    ess_train = 0.8,
    ess_target = n_train * c(0.05, 0.1),
    resampling = 1,
    epoch = 1,
    refresh = 1,
    rand = FALSE
) {
  config_nuts <- tidyr::expand_grid(
    model = models,
    algo = "nuts",
    seed = seeds,
    chains = chains,
    warmup = warmup_nuts,
    iter = iter,
    delta = deltas,
    metric = "diag",
    n_train = NA_real_,
    ess_train = NA_real_,
    ess_target = NA_real_,
    resampling = NA_real_,
    epoch = NA_real_,
    refresh = NA_real_,
    rand = NA
  )
  
  config_ehmc <- tidyr::expand_grid(
    model = models,
    algo = "ehmc",
    seed = seeds,
    chains = chains,
    warmup = warmup_ehmc,
    iter = iter,
    delta = deltas,
    metric = metric,
    n_train = n_train,
    ess_train = ess_train,
    ess_target = ess_target,
    resampling = resampling,
    epoch = epoch,
    refresh = NA_real_,
    rand = NA
  )
  
  config_prehmc <- tidyr::expand_grid(
    model = models,
    algo = "prhmc",
    seed = seeds,
    chains = chains,
    warmup = warmup_ehmc,
    iter = iter,
    delta = deltas,
    metric = "diag",
    n_train = n_train,
    ess_train = ess_train,
    ess_target = ess_target,
    resampling = resampling,
    epoch = epoch,
    refresh = refresh,
    rand = rand
  )
  
  return(dplyr::bind_rows(config_nuts, config_ehmc, config_prehmc))
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

load_one_summary <- function(opt) {
  model_dir <- file.path(opt$repo_root, "examples", opt$model)
  algo_dir <- file.path(model_dir, opt$algo)
  
  basename <- get_basename(opt)
  algo_name <- get_algoname(opt)
  
  diagnosis_file <- if (opt$algo == "nuts") {
    file.path(
      algo_dir,
      paste0("diagnosis_", algo_name, ".RData")
    )
  } else {
    file.path(
      algo_dir,
      paste0("diagnosis_", algo_name, basename, ".RData")
    )
  }
  
  if (!file.exists(diagnosis_file)) {
    message("Missing summary: ", diagnosis_file)
    return(NULL)
  }
  
  env <- new.env(parent = emptyenv())
  load(diagnosis_file, envir = env)
  
  if (!exists("df_summary", envir = env)) {
    warning("No `df_summary` in ", file)
    return(NULL)
  }
  
  get("df_summary", envir = env)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

collect_summaries <- function(config_grid, repo_root) {
  results <- purrr::pmap(
    config_grid,
    function(...) {
      opt <- list(...)
      opt$repo_root <- repo_root
      load_one_summary(opt)
    }
  )
  
  return(
    list(
      param_stats = results %>%
        purrr::keep(~ !is.null(.x)) %>%
        purrr::map_dfr("param_stats"),
      
      diagnosis = results %>%
        purrr::keep(~ !is.null(.x)) %>%
        purrr::map_dfr("diagnosis")
    )
  )
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

add_analysis_labels <- function(df) {
  df %>%
    dplyr::mutate(
      ess_ratio = dplyr::if_else(
        algo == "nuts",
        NA_real_,
        ess_target / n_train
      ),
      x_group = factor(
        dplyr::if_else(
          algo == "nuts",
          "NUTS",
          sprintf("%.3g", ess_ratio)
        ),
        levels = c("NUTS", "0.05", "0.1")
      ),
      delta_label = paste0("p[0]==", delta),
      method_code = dplyr::case_when(
        algo == "nuts" ~ "NUTS",
        algo == "ehmc" & warmup == 2000 ~ "eHMC-2000",
        algo == "ehmc" & warmup == 5000 ~ "eHMC-5000",
        algo == "prhmc" & warmup == 2000 ~ "prHMC-2000",
        algo == "prhmc" & warmup == 5000 ~ "prHMC-5000",
        TRUE ~ NA_character_
      ),
      method_code = factor(
        method_code,
        levels = c("NUTS", "eHMC-2000", "eHMC-5000", "prHMC-2000", "prHMC-5000")
      )
    ) %>%
    dplyr::filter(!is.na(method_code))
}