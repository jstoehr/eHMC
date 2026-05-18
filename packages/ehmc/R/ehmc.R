#' Self-calibrated randomized Hamiltonian Monte Carlo
#'
#' Runs Hamiltonian Monte Carlo with randomized integration times calibrated
#' from weighted warmup samples.
#'
#' @param iter Integer. Number of sampling iterations per chain.
#' @param log_pdf Function. Log-density of the target distribution.
#' @param grad_log_pdf Function. Gradient of `log_pdf`.
#' @param warmup_samples Matrix. Warmup samples used for calibration. By
#'   default, samples are stored by row.
#' @param sampling_weights Numeric vector. Importance weights associated with
#'   `warmup_samples`. If `NULL`, uniform weights are used.
#' @param chains Integer. Number of chains to run. Default is 1.
#' @param init Numeric vector. Initial value for the chains. If `NULL`, each
#'   chain is initialized by sampling from `warmup_samples` according to
#'   `sampling_weights`.
#' @param control List of control parameters:
#' \describe{
#'   \item{`by_row`}{Logical. Whether `warmup_samples` stores samples by row.
#'   Default is `TRUE`.}
#'   \item{`metric`}{Character. Mass matrix type: `"diag"`, `"dense"`, or
#'   `"unit"`. Default is `"diag"`.}
#'   \item{`inv_M`}{Matrix. Initial inverse mass matrix. Default is the identity
#'   matrix. Ignored when `metric = "unit"`.}
#'   \item{`adapt_M`}{Logical. Whether to adapt the mass matrix. Default is
#'   `TRUE`.}
#'   \item{`reg_param`}{Numeric. Regularization parameter used during mass
#'   matrix adaptation. Default is `1e-1`.}
#'   \item{`stepsize`}{Numeric. Initial leapfrog step size. Default is `1e-1`.}
#'   \item{`adapt_stepsize`}{Logical. Whether to adapt the step size. Default is
#'   `TRUE`.}
#'   \item{`epoch_adapt`}{Integer. Number of adaptation epochs. Default is 1.}
#'   \item{`max_length`}{Numeric. Maximum number of leapfrog steps considered
#'   during adaptation. Default is `2e16`.}
#'   \item{`adapt_delta`}{Numeric. Target acceptance probability for step-size
#'   adaptation. Default is `0.651`. Ignored when `adapt_stepsize = FALSE`.}
#'   \item{`refresh_rate`}{Numeric. Probability of refreshing the momentum at
#'   each iteration. Default is 1.}
#'   \item{`rand_leapfrog`}{Logical. Whether to randomize the number of
#'   leapfrog steps when momentum is not refreshed. Default is `TRUE`.}
#'   \item{`min_eps`}{Numeric. Minimum step size allowed during adaptation.
#'   Default is `1e-12`.}
#'   \item{`tol`}{Numeric. Tolerance used when targeting `adapt_delta`.
#'   Default is `0.05`.}
#'   \item{`max_bisect_iter`}{Integer. Maximum number of bisection iterations
#'   used during step-size adaptation. Default is 50.}
#' }
#' @param tuned_pars Optional list of precomputed tuning parameters. If supplied,
#'   calibration is skipped. Expected entries are `weights__`,
#'   `n_leapfrogs__`, `stepsizes__`, and `inv_M__`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`sim`}{Simulation output:
#'     \describe{
#'       \item{`samples`}{List of matrices, one per chain, containing posterior
#'         draws.}
#'       \item{`lp__`}{List of log-density values, one vector per chain.}
#'       \item{`freq_accept`}{Numeric vector of acceptance frequencies, one per
#'         chain.}
#'       \item{`diagnosis`}{List of data frames, one per chain, containing
#'         sampler diagnostics:
#'         \describe{
#'           \item{`accept_stat__`}{Metropolis acceptance probability for each
#'             transition.}
#'           \item{`stepsize__`}{Leapfrog integrator step size used at the
#'             iteration.}
#'           \item{`n_leapfrog__`}{Number of leapfrog steps performed at the
#'             iteration. With partial momentum refreshment, this may be smaller than the
#'             number of leapfrog steps drawn from the empirical distribution,
#'             because previously computed trajectory points can be reused.}
#'           \item{`index__`}{Index of the atom of the discrete empirical
#'             distribution used to sample the step size and the number of
#'             leapfrog steps.}
#'           \item{`move__`}{Logical indicator of whether a transition
#'             is accepted}
#'           \item{`divergent__`}{Logical indicator of whether a divergent
#'             transition occurred.}
#'           \item{`delta_energy_stat__`}{Change in Hamiltonian energy across
#'             the trajectory.}
#'           \item{`energy__`}{Hamiltonian energy of the proposed state.}
#'         }}
#'     }}
#'   \item{`tuned_pars`}{List of tuning parameters used for sampling.}
#'   \item{`adaptation_settings`}{Resolved control settings used for the run.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example to be added.
#' }
#'
#' @export
ehmc <- function(
    iter,
    log_pdf,
    grad_log_pdf,
    warmup_samples,
    sampling_weights = NULL,
    chains = 1,
    init = NULL,
    control = list(),
    tuned_pars = NULL
) {
  # --- Sample format
  default <- list(by_row = TRUE)
  if (!is.null(control$by_row)) {
    if (!is.logical(control$by_row) || length(control$by_row) != 1L) {
      stop("'control$row' must be a boolean (TRUE or FALSE)")
    } else {
      default$by_row <- control$by_row
    }
  }

  # --- Checking inputs
  # ------ sample
  if (!is.matrix(warmup_samples)) {
    stop("'warmup_samples' must be a matrix.")
  }
  n_sample <- ifelse(default$by_row, nrow(warmup_samples), ncol(warmup_samples))
  dim <- ifelse(default$by_row, ncol(warmup_samples), nrow(warmup_samples))
  if (n_sample == 1) {
    stop("'warmup_samples' must contain more than 1 sample.")
  }

  # --- Formating vector
  if (is.null(sampling_weights)) {
    sampling_weights <- rep(1.0, n_sample)/n_sample
  } else if (!is.vector(sampling_weights)) {
    stop("'sampling_weights' must be a numeric vector.")
  } else if (length(sampling_weights) != n_sample) {
    stop("'sampling_weights' length is incompatible. Expected length", n_sample, " (number of samples) but got length ", length(sampling_weights))
  }

  # --- Checking control inputs
  default <- c(default, list(
    metric = "diag",
    inv_M = diag(1, dim),
    adapt_M = TRUE,
    reg_param = 1e-1,
    stepsize = 1e-1,
    adapt_stepsize = TRUE,
    epoch_adapt = 1,
    max_length = 2e16,
    adapt_delta = 0.651,
    refresh_rate = 1.0,
    rand_leapfrog = TRUE,
    min_eps = 1e-12,
    tol = 0.05,
    max_bisect_iter = 50
  ))
  names_to_replace <- checking_names(names(control), "control", default = names(default))
  checking_inputs(names_to_replace, control, default)
  default[names_to_replace] <- control[names_to_replace]

  # --- Metric type is only possible when the user uses one of the adaptation method
  if (length(default$metric) > 1) {
    stop("'control$metric' should be a single value. Allowed values: 'diag', 'dense', or 'unit'.")
  }
  if (!default$metric %in% c("diag", "dense", "unit")) {
    stop("Invalid value for 'control$metric'. Must be one of: 'diag', 'dense', or 'unit'.")
  }

  if (!is.null(control$inv_M)) {
    m_type <- what_metric_type(default$inv_M)
    if (!is.null(control$metric)) {
      if (control$metric == "unit") {
        default$inv_M <- diag(1, dim)
        warning("'control$inv_M' is ignored when 'control$metric' is set to 'unit'")
      } else if (m_type != control$metric) {
        warning("Mismatch between 'control$inv_M' type (", m_type, ") and 'control$metric' (", control$metric, "). Metric forced to ", m_type)
        default$metric <- m_type
      }
    } else {
      default$metric <- m_type
    }
  }

  # --- Warning for adapt_delta being ignored if adapt_stepsize is FALSE
  if(!default$adapt_stepsize && !is.null(control$adapt_delta)) {
    warning("'control$adapt_delta' is ignored when 'control$adapt_stepsize' is set to FALSE.")
  }
  if (default$adapt_delta > 1) { # --- The positivity check has already been done
    stop("'control$adapt_delta' must be between 0 and 1.")
  }

  # --- Adapting HMC parameters
  start_time <- Sys.time()
  if (is.null(tuned_pars)) {
    hmc_pars <- rcpp_adapt_pars(
      if (default$by_row) t(warmup_samples) else warmup_samples,
      sampling_weights,
      log_pdf,
      grad_log_pdf,
      default
    )
  } else {
    hmc_pars <- tuned_pars
  }
  end_time <- Sys.time()
  time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  minutes <- floor(time_diff / 60)
  seconds <- round(time_diff %% 60)

  message(sprintf("Calibration done in: %d minute(s) %d second(s)", minutes, seconds))

  init_user <- init
  samples <- list()
  lp__ <- list()
  freq_accept <- rep(0, chains)
  diagnosis <- list()

  start_time <- Sys.time()
  for(i in seq_len(chains)) {
    message("SAMPLING FROM MODEL: CHAIN ", i)
    # --- Random initialization of the chain
    if(is.null(init_user)) {
      index <- sample(seq_len(n_sample), 1, prob = sampling_weights)
      init <- if (default$by_row) warmup_samples[index, ] else warmup_samples[, index]
    }

    # --- Running HMC
    ans <- rcpp_run(
      iter,
      init,
      hmc_pars$n_leapfrogs__,
      hmc_pars$stepsizes__,
      hmc_pars$weights__,
      log_pdf,
      grad_log_pdf,
      hmc_pars$inv_M__,
      default$metric,
      default$refresh_rate,
      default$rand_leapfrog
    )

    # --- Formating output
    samples[[i]] <- if (default$by_row) t(ans$samples) else ans$samples
    freq_accept[i] <- ans$freq_accept__
    lp__[[i]] = ans$lp__
    diagnosis[[i]] <- data.frame(
      accept_stat__ = pmin(1.0, exp(ans$log_mh_prob__)),
      stepsize__ = ans$stepsize__,
      n_leapfrog__ = ans$n_leapfrog__,
      index__ = ans$index__ + 1,
      move__ = ans$moves__,
      divergent__ = ans$divergent__,
      delta_energy_stat__ = ans$log_mh_prob__,
      energy__ = ans$energy__
    )
  }
  end_time <- Sys.time()
  time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  minutes <- floor(time_diff / 60)
  seconds <- round(time_diff %% 60)

  message(sprintf("Sampling done in: %d minute(s) %d second(s)", minutes, seconds))

  return(
    list(
      sim = list(
        samples = samples,
        lp__ = lp__,
        freq_accept = freq_accept,
        diagnosis = diagnosis
      ),
      tuned_pars = hmc_pars,
      adaptation_settings = default
    )
  )
}
