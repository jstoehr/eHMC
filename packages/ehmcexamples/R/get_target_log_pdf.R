#' Create the target log-density function from a Stan model
#'
#' Returns a function evaluating the unnormalized log-density using
#' `rstan::log_prob()`. The function can either use a supplied `stanfit`
#' object, or create one from a model name and Stan data.
#'
#' @param stanfit Optional `stanfit` object. If supplied, `model` and `data`
#'   are ignored.
#' @param model Optional character string. Name of a Stan model contained in
#'   `get_stanmodels()`. Required only when `stanfit` is `NULL`.
#' @param data Optional list. Stan-formatted data. Required only when
#'   `stanfit` is `NULL`.
#'
#' @return A function taking a parameter vector `pars` and returning the
#'   corresponding unnormalized log-density.
#'
#' @seealso \code{\link{get_stanfit}}, \code{\link[rstan]{log_prob}}
#' @export
get_target_log_pdf <- function(stanfit = NULL, model = NULL, data = NULL) {
  if (is.null(stanfit)) {
    if (is.null(model) || is.null(data)) {
      stop("Either provide `stanfit`, or provide both `model` and `data`.")
    }
    stanfit <- get_stanfit(model, data)
  }

  f <- function(pars) {
    rstan::log_prob(stanfit, pars, FALSE)
  }

  return(f)
}
