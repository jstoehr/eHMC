# --------------------------------------------------
# --- LIBRARIES AND OPTIONS
# --------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(arrow))
suppressMessages(library(dplyr))
suppressMessages(library(ehmc))
suppressMessages(library(ehmcexamples))

option_list <- list(
  make_option("--repo_root", type = "character", default = ".", 
              help = "Root directory of the repository (default: current directory)"),
  make_option("--model", type = "character",
              help = "Model name (e.g., 'banana', 'BLP', 'MVNorm')"),
  make_option("--algo", type = "character", default = "ehmc",
              help = "Sampling algorithm to use (options: 'ehmc', 'prhmc', 'nuts'; 
              default: 'ehmc')"),
  make_option("--seed", type = "integer", default = 1,
              help = "Random seed for reproducibility (default: 1)"),
  make_option("--chains", type = "integer", default = 50,
              help = "Number of chains to run (default: 50)"),
  make_option("--warmup", type = "integer", default = 2000,
              help = "Number of warmup iterations in nuts / of warmup samples 
              in ehmc/prhmc (default: 1000)"),
  make_option("--iter", type = "integer", default = 2000,
              help = "Number of iterations to run the sampler (it does not include 
              the warmup for nuts. Default: 2000)"),
  make_option("--delta", type = "double", default = 0.651, 
              help = "Target acceptance probability (default: 0.651)"),
  make_option("--metric", type = "character", default = "diag",
              help = "Metric used for the sampler (default: 'diag')"),
  make_option("--n_train", type = "integer", default = 50000,
              help = "Number of draws used for training the proposal used in 
              ehmc/prhmc (default: 50000)"),
  make_option("--ess_train", type = "double", default = 0.8,
              help = "Target ESS for training expressed as a fraction of n_train (default means 0.8 * n_train)"),
  make_option("--ess_target", type = "double", default = NULL,
              help = "Target ESS (default: None)"),
  make_option("--resampling", type = "integer", default = 1,
              help = "Resampling step from proposal (default: 1)"),
  make_option("--epoch", type = "integer", default = 1,
              help = "Number of adaptation epochs (default: 1)"),
  make_option("--refresh", type = "double", default = 1.0,
              help = "Momentum refreshment rate (default: 1.0)"),
  make_option("--rand", type = "integer", default = 1,
              help = "Whether to randomize the number of leapfrog steps when 
              momentum is not refreshed (default: 1)")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$repo_root <- normalizePath(opt$repo_root)
opt$algo <- match.arg(opt$algo, c("ehmc", "prhmc", "nuts"))
opt$metric <- match.arg(opt$metric, c("diag", "dense", "unit"))

if (is.null(opt$model)) stop("--model is required.")

if (is.null(opt$ess_target) && opt$algo != "nuts") {
  opt$ess_target <- floor(0.8 * opt$n_train)
}

opt$rand <- opt$rand > 0

# --------------------------------------------------
# --- SOURCE R FUNCTIONS AND MODELS
# --------------------------------------------------

source(file.path(opt$repo_root, "examples", opt$model, "model.R"))
source(file.path(opt$repo_root, "examples", "R", "io_helpers.R"))
source(file.path(opt$repo_root, "examples", "R", "algo_helpers.R"))
source(file.path(opt$repo_root, "examples", "R", "summary_helpers.R"))

# --------------------------------------------------
# --- WORKING DIRECTORIES AND FILENAMES
# --------------------------------------------------

model_dir <- file.path(opt$repo_root, "examples", opt$model)
pars_dir <- file.path(model_dir, "tuned_pars")
algo_dir <- file.path(model_dir, opt$algo)

dir.create(pars_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(algo_dir, recursive = TRUE, showWarnings = FALSE)

basename <- get_basename(opt)
algo_name <- get_algoname(opt)

rdata_file <- function(root, ...) {
  file.path(root, paste0(..., ".RData"))
}

algo_file <- rdata_file(root = algo_dir, algo_name, basename)

diagnosis_file <- if (opt$algo == "nuts") {
  rdata_file(root = algo_dir, "diagnosis_", algo_name)
} else {
  rdata_file(root = algo_dir, "diagnosis_", algo_name, basename)
}

tuned_pars_file <- rdata_file(
  root = pars_dir,
  basename,
  "_delta_",
  opt$delta,
  "_tuned_pars"
)

# --------------------------------------------------
# --- RUN THE SAMPLER
# --------------------------------------------------

if (!file.exists(algo_file)) {
  
  data <- get_data(data_dir = model_dir)
  fit <- ehmcexamples::get_stanfit(opt$model, data)
  
  target_log_pdf <- function(pars) {
    rstan::log_prob(fit, pars, FALSE)
  }
  
  target_grad_log_pdf <- function(pars) {
    rstan::grad_log_prob(fit, pars, FALSE)
  }
  
  set.seed(opt$seed)
  
  if (opt$algo == "nuts") {
    output <- run_nuts(opt, data)
  } else {
    warmup <- load_warmup_sample(opt, model_dir)
    
    if (file.exists(tuned_pars_file)) {
      load(tuned_pars_file)
    } else {
      tuned_pars <- NULL
    }
    
    output <- run_ehmc(
      opt,
      target_log_pdf,
      target_grad_log_pdf,
      warmup$samples,
      warmup$weights,
      tuned_pars = tuned_pars
    )
    
    if (!file.exists(tuned_pars_file)) {
      tuned_pars <- output$tuned_pars
      save(tuned_pars, file = tuned_pars_file)
    }
  }
  
  save(opt, output, file = algo_file)
  message(opt$algo, " sample done for ", opt$model, " with seed ", opt$seed, " and delta ", opt$delta,
          ". Results saved in ", algo_file)
  
  # --------------------------------------------------
  # --- COMPUTE SUMMARIES AND DIAGNOSTICS
  # --------------------------------------------------
  
  df_summary <- get_summary(opt, output = output, data_dir = model_dir)
  
  save(df_summary, file = diagnosis_file)
  message(opt$algo, " diagnosis done for ", opt$model, " with seed ", opt$seed, " and delta ", opt$delta,
          ". Results saved in ", diagnosis_file)
} else {
  message(algo_file, " already exists.")
}
