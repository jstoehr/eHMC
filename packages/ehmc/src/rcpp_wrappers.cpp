// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <variant> // C++17
#include <iostream>
#include <utils.hpp>
#include <algo.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_get_lp__(const arma::mat& thetas_by_col,
                                  const Rcpp::Function& target_log_pdf) {
  // --- Containers
  arma::vec lp__(thetas_by_col.n_cols);
  for (size_t r = 0; r < thetas_by_col.n_cols; r++) {
    lp__.at(r) = Rcpp::as<double>(target_log_pdf(thetas_by_col.col(r)));
  }

  return Rcpp::NumericVector(lp__.begin(), lp__.end());
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List rcpp_get_sampling_weights(const arma::mat& thetas_by_col,
                                     const arma::vec& proposal_log_pdf_val,
                                     const Rcpp::Function& target_log_pdf) {
  // --- Containers
  arma::vec log_weights(thetas_by_col.n_cols);
  double log_weights_max = 0.0;
  for (size_t r = 0; r < thetas_by_col.n_cols; r++) {
    log_weights.at(r) = Rcpp::as<double>(target_log_pdf(thetas_by_col.col(r))) -
      proposal_log_pdf_val.at(r);
    if (log_weights.at(r) > log_weights_max or r == 0) {
      log_weights_max = log_weights.at(r);
    }
  }
  log_weights -= log_weights_max;
  arma::vec weights = arma::exp(log_weights);
  weights /= arma::sum(weights);
  log_weights += log_weights_max;

  return Rcpp::List::create(
    Rcpp::Named("weights", Rcpp::NumericVector(weights.begin(), weights.end())),
    Rcpp::Named("log_weights", Rcpp::NumericVector(log_weights.begin(), log_weights.end()))
  );
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

AdaptControl rcpp_parse_adapt_control(const Rcpp::List& control) {
  AdaptControl out;
  out.metric_str = Rcpp::as<std::string>(control["metric"]);
  out.max_length = Rcpp::as<size_t>(control["max_length"]);
  out.adapt_stepsize = Rcpp::as<bool>(control["adapt_stepsize"]);
  out.stepsize = Rcpp::as<double>(control["stepsize"]);
  out.adapt_m = Rcpp::as<bool>(control["adapt_M"]);
  out.epoch_adapt = Rcpp::as<size_t>(control["epoch_adapt"]);
  out.adapt_delta = Rcpp::as<double>(control["adapt_delta"]);
  out.reg_param = Rcpp::as<double>(control["reg_param"]);
  return out;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List rcpp_adapt_pars(const arma::mat& thetas_by_col,
                           const arma::vec& weights,
                           const Rcpp::Function& target_log_pdf,
                           const Rcpp::Function& target_grad_log_pdf,
                           const Rcpp::List& control) {
  // --- Extract control parameters from the input list
  arma::mat inv_m = control["inv_M"];
  AdaptControl ctrl = rcpp_parse_adapt_control(control);

  // --- Initialize containers for adaptation
  arma::vec weights_ = weights/arma::sum(weights); // To be sure weights are properly normalised
  type_metric metric = string_to_metric(ctrl.metric_str);
  arma::Col<size_t> leaps_vector(2 * weights.n_elem);
  arma::vec stepsizes_vector(2 * weights.n_elem, arma::fill::value(ctrl.stepsize));
  arma::vec mh_vector(2 * weights.n_elem);
  arma::mat inv_m_alt = inv_m;
  arma::mat thetas_by_col_new = thetas_by_col;

  // --- Define the variant to hold the different algo types
  std::variant<algo<unit>, algo<diag>, algo<dense>> algorithm;

  // --- Choose the correct algorithm variant based on the metric
  switch (metric) {
  case 0: algorithm = algo<unit>(target_log_pdf, target_grad_log_pdf, inv_m, ctrl); break;
  case 1: algorithm = algo<diag>(target_log_pdf, target_grad_log_pdf, inv_m, ctrl); break;
  case 2: algorithm = algo<dense>(target_log_pdf, target_grad_log_pdf, inv_m, ctrl); break;
  default: Rcpp::stop("Unknown metric type.");
  }

  // --- Use std::visit to call the appropriate adapt_pars method on the chosen algorithm
  std::visit([&](auto& algo_instance) {
    algo_instance.adapt_pars(thetas_by_col, weights_, leaps_vector, stepsizes_vector,
                             mh_vector, inv_m_alt);
    inv_m = algo_instance.inv_m; // Update inv_m with the result
  }, algorithm);

  // --- Reshape weight to account for forward and backward exploration
  weights_.set_size(2 * weights.n_elem);
  for (arma::uword k = 0; k < weights.n_elem; k++) {
    weights_.at(k) = weights.at(k);
    weights_.at(k + weights.n_elem) = weights.at(k);
  }
  weights_ /= arma::sum(weights_);

  return Rcpp::List::create(
    Rcpp::Named("weights__", Rcpp::NumericVector(weights_.begin(), weights_.end())),
    Rcpp::Named("n_leapfrogs__", Rcpp::NumericVector(leaps_vector.begin(), leaps_vector.end())),
    Rcpp::Named("stepsizes__", Rcpp::NumericVector(stepsizes_vector.begin(), stepsizes_vector.end())),
    Rcpp::Named("accept_stat__", Rcpp::NumericVector(mh_vector.begin(), mh_vector.end())),
    Rcpp::Named("inv_M__", inv_m),
    Rcpp::Named("inv_M_alt__", inv_m_alt)
  );
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List rcpp_run(size_t iter,
                    const arma::vec& start,
                    const arma::Col<size_t>& leaps_vector,
                    const arma::vec& stepsizes_vector,
                    const arma::vec& weights,
                    const Rcpp::Function& target_log_pdf,
                    const Rcpp::Function& target_grad_log_pdf,
                    const arma::mat& inv_m,
                    std::string metric_str,
                    double refresh_rate,
                    bool random) {
  // --- Containers
  arma::vec weights_ = weights/arma::sum(weights); // To be sure weights are properly normalised
  type_metric metric = string_to_metric(metric_str);
  arma::mat thetas_by_col(start.n_rows, iter + 1);
  thetas_by_col.col(0) = start;
  arma::vec lp(iter + 1);
  arma::vec energy(iter);
  arma::vec selected_eps(iter);
  arma::Col<size_t> selected_l(iter);
  arma::uvec indices(iter);
  arma::uvec moves = arma::ones<arma::uvec>(iter);
  arma::uvec divergent = arma::zeros<arma::uvec>(iter);
  arma::vec log_mh_prob(iter);

  // --- Define the variant to hold the different algo types
  std::variant<algo<unit>, algo<diag>, algo<dense>> algorithm;

  // --- Choose the correct algorithm variant based on the metric
  switch (metric) {
  case 0: algorithm = algo<unit>(target_log_pdf, target_grad_log_pdf, inv_m); break;
  case 1: algorithm = algo<diag>(target_log_pdf, target_grad_log_pdf, inv_m); break;
  case 2: algorithm = algo<dense>(target_log_pdf, target_grad_log_pdf, inv_m); break;
  default: Rcpp::stop("Unknown metric type.");
  }

  // --- Use std::visit to call the appropriate run method on the chosen algorithm
  std::visit([&](auto& algo_instance) {
    algo_instance.run(iter, start, leaps_vector, stepsizes_vector, weights_,
                      thetas_by_col, lp, energy, selected_eps, selected_l, indices, moves, log_mh_prob, divergent,
                      refresh_rate, random);
  }, algorithm);

  return Rcpp::List::create(
    Rcpp::Named("samples", thetas_by_col),
    Rcpp::Named("lp__", Rcpp::NumericVector(lp.begin(), lp.end())),
    Rcpp::Named("energy__", Rcpp::NumericVector(energy.begin(), energy.end())),
    Rcpp::Named("freq_accept__", arma::sum(moves)/static_cast<double>(iter)),
    Rcpp::Named("log_mh_prob__", Rcpp::NumericVector(log_mh_prob.begin(), log_mh_prob.end())),
    Rcpp::Named("stepsize__", Rcpp::NumericVector(selected_eps.begin(), selected_eps.end())),
    Rcpp::Named("n_leapfrog__", Rcpp::NumericVector(selected_l.begin(), selected_l.end())),
    Rcpp::Named("index__", Rcpp::NumericVector(indices.begin(), indices.end())),
    Rcpp::Named("moves__", Rcpp::NumericVector(moves.begin(), moves.end())),
    Rcpp::Named("divergent__", Rcpp::NumericVector(divergent.begin(), divergent.end()))
  );
}
