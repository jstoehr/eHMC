#' Initialize a Stan fit object for a given example model
#'
#' Creates a minimal `stanfit` object using one fixed-parameter iteration.
#' This is primarily used to access `rstan::log_prob()` and related methods
#' for the selected model without running full posterior sampling.
#'
#' @param model Character string. Name of the Stan model contained in
#'   `get_stanmodels()`.
#' @param data List. Stan-formatted data corresponding to the selected model.
#'
#' @return A `stanfit` object for the specified model.
#'
#' @seealso \code{\link{get_stanmodels}}
#' @export
get_stanfit <- function(model, data) {
  # --- Set the model
  stanmodels <- get_stanmodels()

  stanfit <- rstan::sampling(
    stanmodels[[model]],
    data = data,
    chains = 1,
    iter = 1,
    warmup = 0,
    algorithm = "Fixed_param"
  )
  return(stanfit)
}
