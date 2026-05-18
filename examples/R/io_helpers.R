get_basename <- function(opt) {
  return(
    sprintf(
      "seed_%s_train_%s_%s_ess_%s_calib_ehmc_w_%s_r_%s",
      opt$seed,
      opt$n_train,
      100 * opt$ess_train,
      opt$ess_target,
      opt$warmup,
      opt$resampling
    )
  )
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

get_algoname <- function(opt) {
  base <- sprintf(
    "%s_seed_%s_c_%s_w_%s_i_%s_delta_%s_metric_%s",
    opt$algo,
    opt$seed,
    opt$chains,
    opt$warmup,
    opt$iter,
    opt$delta,
    opt$metric
  )
  if (opt$algo != "nuts") {
    base <- sprintf("%s_epoch_%s_", base, opt$epoch)
  }
  if (opt$algo == "prhmc") {
    base <- sprintf("%srefresh_%s_rand_%s_", base, opt$refresh, opt$rand)
  }
  
  return(base)
}

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

load_warmup_sample <- function(opt, model_dir) {
  warmup_dir <- if (opt$resampling == 0) {
    file.path(model_dir, "warmup_wo_resampling")
  } else {
    file.path(model_dir, "warmup_with_resampling")
  }
  
  basename <- get_basename(opt)
  path <- file.path(warmup_dir, paste0(basename, ".parquet"))
  if (!file.exists(path)) {
    stop("Warmup file not found: ", path)
  }
  
  df <- arrow::read_parquet(path)
  
  return(
    list(
      samples = as.matrix(dplyr::select(df, -w)), 
      weights = as.numeric(df$w)
    )
  )
}
