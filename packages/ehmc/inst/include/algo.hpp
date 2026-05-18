#ifndef INST_INCLUDE_ALGO_HPP_
#define INST_INCLUDE_ALGO_HPP_

#include <progress_bar.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

struct AdaptControl {
  std::string metric_str = "diag";
  double stepsize = 0.1;
  size_t max_length = 2e16;
  bool adapt_stepsize = true;
  bool adapt_m = true;
  size_t epoch_adapt = 1;
  double adapt_delta = 0.6;
  double reg_param = 1e-1;
  double min_eps = 1e-12;
  double tol = 0.05;
  size_t max_bisect_iter = 50;
};

template<type_metric T_metric>
class algo {
public:
  // --- Constructors
  algo() = default; // required for defining variants

  algo(
    const Rcpp::Function& target_log_pdf_,
    const Rcpp::Function& target_grad_log_pdf_,
    const arma::mat& inv_m_,
    const AdaptControl& adapt_control_ = AdaptControl()
  ):target_log_pdf(&target_log_pdf_), target_grad_log_pdf(&target_grad_log_pdf_),
  adapt_control(adapt_control_) {
    get_matrix_decomposition<T_metric>(inv_m_, inv_m, chol_m, inv_chol_m, 0.0);
  };

  // --- Destructor
  virtual ~algo(){};

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void leapfrog_step(double eps,
                     arma::vec& theta,
                     arma::vec& v,
                     arma::vec& grad_log_pdf) {
    // --- A single leapfrog step
    v += 0.5 * eps * grad_log_pdf;
    // theta += eps * inv_m * v;
    theta += eps * template_mat_mult<T_metric>(inv_m, v);
    grad_log_pdf = Rcpp::as<arma::vec>((*target_grad_log_pdf)(theta));
    v += 0.5 * eps * grad_log_pdf;
  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void leapfrog_integrator(size_t n_leap,
                           double eps,
                           const arma::vec& grad_log_pdf,
                           arma::vec& theta_new,
                           arma::vec& v_new,
                           arma::vec& grad_log_pdf_new) {
    // --- n_leap leapfrog steps
    // --- Save computation cost by reducing the number of vector summations
    v_new += 0.5 * eps * grad_log_pdf;
    for (size_t l = 0; l < n_leap; l++) {
      // theta_new += eps * inv_m * v_new;
      theta_new += eps * template_mat_mult<T_metric>(inv_m, v_new);
      grad_log_pdf_new = Rcpp::as<arma::vec>((*target_grad_log_pdf)(theta_new));
      v_new += eps * grad_log_pdf_new;
    }
    v_new -= 0.5 * eps * grad_log_pdf_new;
  }

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void longest_batch(
      double weight,
      const arma::vec& theta,
      const arma::vec& v,
      double log_pdf,
      const arma::vec& grad_log_pdf,
      double eps,
      size_t& l,
      double& mh_prob,
      double& sum_w,
      arma::vec& m1_temp,
      arma::mat& m2_temp
  ) {
    // --- Compute the longest excursion on a level set
    // --- Compute the acceptance rate along the level set for a given stepsize
    // --- Recycle visited points to estimate the mass matrix

    // --- Containers
    bool cond = false;
    // ------ Containers for updated position and related log pdf
    arma::vec theta_new = theta;
    double log_pdf_new = log_pdf;
    arma::vec grad_log_pdf_new = grad_log_pdf;
    // ------ Containers for updated momentum
    arma::vec v_new = v;
    // --- Metropolis Hastings acceptance probability
    double prob = 0.0;
    // --- Precomputed elements for updating the mass matrix
    arma::mat theta_2 = theta * theta.t();

    // --- Initialization: compute the potential of the momentum
    arma::vec temp = v;
    inplace_tri_mat_mult<T_metric>(inv_chol_m, temp);
    double pot_v = 0.5 * arma::dot(temp, temp);

    // --- Reset variables
    l = 0;
    mh_prob = 0.0;
    while(!cond and l < adapt_control.max_length) {
      // --- Increase the number of leapfrog steps till "turning back"
      ++l;
      this->leapfrog_step(eps, theta_new, v_new, grad_log_pdf_new);
      cond = arma::dot(theta_new - theta, template_mat_mult<T_metric>(inv_m, v_new)) < 0.0;

      if (adapt_control.adapt_stepsize or adapt_control.adapt_m) {
        // --- Local Metropolis-Hastings acceptance probability
        log_pdf_new = Rcpp::as<double>((*target_log_pdf)(theta_new));
        temp = v_new;
        inplace_tri_mat_mult<T_metric>(inv_chol_m, temp);
        prob = std::min(1.0, exp(log_pdf_new - log_pdf + pot_v - 0.5 * arma::dot(temp, temp)));

        // --- For adaptithe stepsize: average the acceptance rate
        mh_prob += prob;

        // --- Adapt m (importance sampling estimate)
        if (adapt_control.adapt_m) {
          // --- Local update of the first moment
          m1_temp += weight * (prob * theta_new + (1.0 - prob) * theta);
          // --- Local update of the second moment
          m2_temp += weight * (prob * theta_new * theta_new.t() + (1.0 - prob) * theta_2);
          // --- Number of local updates
          sum_w += weight;
        }
      }
    }
    mh_prob /= static_cast<double>(l);
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void adapt_pars_1(
      double weight,
      const arma::vec& theta,
      const arma::vec& v,
      size_t& l,
      double& eps,
      double& accept_stat,
      double& sum_w,
      arma::vec& m1_temp,
      arma::mat& m2_temp
  ) {
    // --- Adaptation method

    // --- Model containers
    double log_pdf = Rcpp::as<double>((*target_log_pdf)(theta));
    arma::vec grad_log_pdf = Rcpp::as<arma::vec>((*target_grad_log_pdf)(theta));

    // --- Metropolis Hastings acceptance probability
    double mh_prob = 0.0;

    // --- Local lambda function
    auto longest_batch_safe = [&]() -> bool {
      const size_t l_buffer = l;
      const double eps_buffer = eps;
      const double sum_w_buffer = sum_w;
      const arma::vec m1_temp_buffer = m1_temp;
      const arma::mat m2_temp_buffer = m2_temp;

      try {
        this->longest_batch(
            weight, theta, v, log_pdf, grad_log_pdf,
            eps, l,
            mh_prob, sum_w, m1_temp, m2_temp
        );
        return true;
      } catch (const std::exception& e) {
        l = l_buffer;
        eps = eps_buffer;
        sum_w = sum_w_buffer;
        m1_temp = m1_temp_buffer;
        m2_temp = m2_temp_buffer;

        Rcpp::warning("longest_batch failed: %s", e.what());
        return false;
      } catch (...) {
        l = l_buffer;
        eps = eps_buffer;
        sum_w = sum_w_buffer;
        m1_temp = m1_temp_buffer;
        m2_temp = m2_temp_buffer;

        Rcpp::warning("longest_batch failed with an unknown error.");
        return false;
      }
    };

    // --- Compute the longest excursion on the level set
    bool success = false;

    while (!success) {
      success = longest_batch_safe();

      if (!success) {
        if (adapt_control.adapt_m) {
          arma::mat temp_inv_m = 0.5 * inv_m;
          get_matrix_decomposition<T_metric>(
            temp_inv_m, inv_m, chol_m, inv_chol_m, adapt_control.reg_param
          );
        }

        if (adapt_control.adapt_stepsize) {
          eps *= 0.5;

          if (eps < adapt_control.min_eps) {
            Rcpp::warning("Unable to build the longest excursion.");
            break;
          }
        } else {
          Rcpp::warning(
            "Unable to build the longest excursion for the user-specified step size."
          );
          break;
        }
      }
    }

    if (!success) {
      // --- We do not have integration time for this point
      l = 0;
      eps = 0.0;
      accept_stat = 0.0;
    } else if (adapt_control.adapt_stepsize) {
      // --- Bisection method to find the stepsize
      accept_stat = mh_prob;

      // --------------------------------------------------
      // --- Case 1: acceptance too low, decrease eps by bisection
      // --------------------------------------------------
      if (mh_prob < adapt_control.adapt_delta - adapt_control.tol and eps > adapt_control.min_eps) {
        double eps_a = 0.0;
        double eps_b = eps;
        for (size_t i = 0; i < adapt_control.max_bisect_iter; ++i) {
          eps = 0.5 * (eps_a + eps_b);
          this->longest_batch(weight, theta, v, log_pdf, grad_log_pdf,
                              eps, l,
                              mh_prob, sum_w, m1_temp, m2_temp);
          accept_stat = mh_prob;
          if (std::abs(mh_prob - adapt_control.adapt_delta) < adapt_control.tol or std::abs(eps_a - eps_b) < 1e-3 * eps_a) {
            break;
          }

          if (mh_prob < adapt_control.adapt_delta - adapt_control.tol) {
            eps_b = eps;
          } else {
            eps_a = eps;
          }
        }
      } else if (mh_prob > adapt_control.adapt_delta + adapt_control.tol) {
        // --------------------------------------------------
        // --- Case 2: acceptance too high, increase eps by bisection
        // --------------------------------------------------
        double eps_a = eps;
        double eps_b = eps * l;

        for (size_t i = 0; i < adapt_control.max_bisect_iter; ++i) {
          eps = 0.5 * (eps_a + eps_b);
          success = longest_batch_safe();

          if (success) {
            accept_stat = mh_prob;
            if (std::abs(mh_prob - adapt_control.adapt_delta) < adapt_control.tol or std::abs(eps_a - eps_b) < 1e-3 * eps_a) {
              break;
            }

            if (mh_prob < adapt_control.adapt_delta - adapt_control.tol) {
              eps_b = eps;
            } else {
              // --- We found a larger value for eps such that
              // --- the acceptance rate is still large enough
              eps_a = eps;
            }
          } else {
            eps_b = eps;
          }
        }
      }
    }
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void adapt_pars(
      const arma::mat& thetas_by_col,
      const arma::vec& weights,
      arma::Col<size_t>& leaps_vector,
      arma::vec& stepsizes_vector,
      arma::vec& mh_vector,
      arma::mat& inv_m_alt
  ) {
    // --- Containers
    arma::vec v = arma::zeros(thetas_by_col.n_rows);
    double sum_w = arma::accu(weights);
    arma::vec m1 = arma::zeros(thetas_by_col.n_rows);
    arma::mat m2 = arma::zeros<arma::mat>(thetas_by_col.n_rows, thetas_by_col.n_rows);
    arma::mat temp_inv_m = arma::zeros<arma::mat>(thetas_by_col.n_rows, thetas_by_col.n_rows);
    arma::mat temp_m = arma::zeros<arma::mat>(thetas_by_col.n_rows, thetas_by_col.n_rows);

    // --- Adapting parameter

    // --------------------------------------------------
    // --- Adapting Mass Matrix
    // --------------------------------------------------
    if (adapt_control.adapt_m) {
      // --- Importance sampling estimates
      for (arma::uword r = 0; r < weights.n_elem; r++) {
        m1 += weights.at(r) * thetas_by_col.col(r);
        m2 += weights.at(r) * thetas_by_col.col(r) * thetas_by_col.col(r).t();
      }
    }

    // --------------------------------------------------
    // --- Adapting n_leapfrog and stepsize
    // --------------------------------------------------
    for (size_t k = 0; k < adapt_control.epoch_adapt; k++) {
      std::cout << " Adaptation no. " << k + 1 << std::endl;
      if (adapt_control.adapt_m) {
        std::cout << " Adapt mass matrix...";
        m1 /= sum_w;
        m2 /= sum_w;
        temp_inv_m = m2 - m1 * m1.t();

        // --- Compute and store decompositions of the mass matrix
        get_matrix_decomposition<T_metric>(temp_inv_m, inv_m, chol_m, inv_chol_m, adapt_control.reg_param);
        std::cout << " DONE!" << std::endl;

        // --- Reset containers
        sum_w = 0.0;
        m1.zeros();
        m2.zeros();
      }
      std::cout << " Adapt leapfrog parameters..." << std::endl;
      double eps = stepsizes_vector.at(0);
      // --- Progress bar
      progress_bar pb(weights.n_elem);
      for (size_t r = 0; r < weights.n_elem; r++) {
        pb.show_progress();
        // --- Sampling momentum
        v.randn();
        inplace_tri_mat_mult<T_metric>(chol_m, v);
        // --- Adaptation in v direction
        stepsizes_vector.at(r) = eps;
        this->adapt_pars_1(weights.at(r), thetas_by_col.col(r), v,
                           leaps_vector.at(r), stepsizes_vector.at(r), mh_vector.at(r),
                           sum_w, m1, m2);
        // --- Adaptation in -v direction
        arma::uword r_ = r + weights.n_elem;
        stepsizes_vector.at(r_) = stepsizes_vector.at(r);
        this->adapt_pars_1(weights.at(r), thetas_by_col.col(r), -1.0 * v,
                           leaps_vector.at(r_), stepsizes_vector.at(r_), mh_vector.at(r_),
                           sum_w, m1, m2);
        eps = (2 * r * eps + stepsizes_vector.at(r) + stepsizes_vector.at(r_))/static_cast<double>(2 * r + 2);
      }
      pb.finish();
    }

    // --------------------------------------------------
    // --- Adapting Mass Matrix
    // --------------------------------------------------
    if (adapt_control.adapt_m) {
      m1 /= sum_w;
      m2 /= sum_w;
      temp_inv_m = m2 - m1 * m1.t();
      arma::mat chol_m_alt = chol_m;
      arma::mat inv_chol_m_alt = inv_chol_m;
      get_matrix_decomposition<T_metric>(temp_inv_m, inv_m_alt, chol_m_alt, inv_chol_m_alt, adapt_control.reg_param);
    }
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void hmc(
      const arma::vec& theta,
      size_t n_leaps,
      double stepsize,
      double& log_pdf,
      arma::vec& grad_log_pdf,
      arma::vec& theta_new,
      double& lp,
      double& energy,
      unsigned& moves,
      double& log_mh_prob,
      unsigned& divergent
  ) {
    // --- Model containers
    double log_pdf_new = 0.0;
    arma::vec grad_log_pdf_new = arma::zeros(grad_log_pdf.n_rows);
    // --- Momentum containers
    arma::vec v(theta.n_rows);
    arma::vec v_new(theta.n_rows);

    // --- Sampling momentum variable
    v = arma::randn(theta.n_rows);
    inplace_tri_mat_mult<T_metric>(chol_m, v);

    // --- Leapfrog integrator
    theta_new = theta;
    v_new = v;

    // --- Energy of the current iteration
    inplace_tri_mat_mult<T_metric>(inv_chol_m, v);
    energy = - log_pdf + 0.5 * arma::dot(v, v);
    try {
      this->leapfrog_integrator(n_leaps, stepsize, grad_log_pdf,
                                theta_new, v_new, grad_log_pdf_new);
      // --- Metropolis Hastings step
      inplace_tri_mat_mult<T_metric>(inv_chol_m, v_new);
      log_pdf_new = Rcpp::as<double>((*target_log_pdf)(theta_new));
      log_mh_prob = energy + log_pdf_new - 0.5 * arma::dot(v_new, v_new);

      if (log(arma::randu()) < log_mh_prob) {
        // --- Accept the move
        log_pdf = log_pdf_new;
        grad_log_pdf = grad_log_pdf_new;
      } else {
        // --- Reject the move
        theta_new = theta;
        moves = 0;
      }
    } catch (...) {
      // --- Reject the move
      theta_new = theta;
      log_mh_prob = -std::numeric_limits<double>::infinity();
      moves = 0;
      divergent = 1;
      Rcpp::warning("Divergent transition.");
    }
    // --- Storing the log_pdf at the current point
    lp = log_pdf;
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void ehmc(
      size_t iter,
      const arma::vec& start,
      const arma::Col<size_t>& leaps_vector,
      const arma::vec& stepsizes_vector,
      const arma::vec& weights,
      arma::mat& thetas_by_col,
      arma::vec& lp,
      arma::vec& energy,
      arma::vec& selected_eps,
      arma::Col<size_t>& selected_l,
      arma::uvec& selected_ind,
      arma::uvec& moves,
      arma::vec& log_mh_prob,
      arma::uvec& divergent
  ) {
    // --- Model containers
    double log_pdf = Rcpp::as<double>((*target_log_pdf)(start));
    lp.at(0) = log_pdf;
    arma::vec grad_log_pdf = Rcpp::as<arma::vec>((*target_grad_log_pdf)(start));
    arma::vec theta;
    arma::vec theta_new;
    // --- Sampling index container
    arma::uword ind;

    std::cout << " eHMC run..." << std::endl;
    progress_bar pb(iter);
    for (size_t k = 0; k < iter; k++) {
      pb.show_progress();

      // --- Sampling integration time
      ind = sample_int_1(weights);
      selected_eps.at(k) = stepsizes_vector.at(ind);
      selected_l.at(k) = leaps_vector.at(ind);
      selected_ind.at(k) = ind;

      // --- HMC steps
      theta = thetas_by_col.col(k);
      this->hmc(theta,
                leaps_vector.at(ind),
                stepsizes_vector.at(ind),
                log_pdf,
                grad_log_pdf,
                theta_new,
                lp.at(k + 1),
                energy.at(k),
                moves.at(k),
                log_mh_prob.at(k),
                divergent.at(k)
      );
      thetas_by_col.col(k + 1) = theta_new;
    }
    std::cout << " DONE!" << std::endl;
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void pr_ehmc(
      size_t iter,
      const arma::vec& start,
      const arma::Col<size_t>& leaps_vector,
      const arma::vec& stepsizes_vector,
      const arma::vec& weights,
      double refresh_rate,
      bool random,
      arma::mat& thetas_by_col,
      arma::vec& lp,
      arma::vec& energy,
      arma::vec& selected_eps,
      arma::Col<size_t>& selected_l,
      arma::uvec& selected_ind,
      arma::uvec& moves,
      arma::vec& log_mh_prob,
      arma::uvec& divergent
  ) {
    // --- Model containers
    double log_pdf = Rcpp::as<double>((*target_log_pdf)(start));
    lp.at(0) = log_pdf;
    double log_pdf_new = 0.0;
    arma::vec grad_log_pdf = Rcpp::as<arma::vec>((*target_grad_log_pdf)(start));
    arma::vec grad_log_pdf_new(grad_log_pdf.n_rows);

    // --- Momentum containers
    arma::vec v_new(start.n_rows);

    // --- Sampling initial momentum
    arma::vec v = arma::randn(start.n_rows);
    inplace_tri_mat_mult<T_metric>(chol_m, v);

    // --- Sampling index container
    arma::uword ind = sample_int_1(weights);

    // --- Containers to handle partial momentum refreshment and recycling
    size_t temp_l = 0;
    arma::uword loc = 0;
    arma::sword lag = 0;
    arma::uword shift = 0;
    int sigma = 1;
    // --- Position containers
    arma::vec theta_new(start.n_rows); // update
    arma::mat temp_t;
    arma::mat mat_t(start.n_rows, 1); // storing intermediates points
    mat_t.col(0) = start;
    // --- Momentum containers
    arma::mat temp_v;
    arma::mat mat_v(start.n_rows, 1); // storing intermediates points
    mat_v.col(0) = v;

    std::cout << " prHMC run..." << std::endl;

    progress_bar pb(iter);
    for (size_t k = 0; k < iter; k++) {
      pb.show_progress();

      if (arma::randu() < refresh_rate or k < 1) {
        // --------------------------------------------------
        // --- The momentum is refreshed
        // --------------------------------------------------

        // --- Sampling number of leapfrog steps
        ind = sample_int_1(weights);
        selected_eps.at(k) = stepsizes_vector.at(ind);
        selected_l.at(k) = leaps_vector.at(ind);
        selected_ind.at(k) = ind;

        // --- Resampling momentum variable
        v = arma::randn(start.n_rows);
        inplace_tri_mat_mult<T_metric>(chol_m, v);

        // --- Initializing container position and momentum
        mat_t.zeros(start.n_rows, leaps_vector.at(ind) + 1);
        mat_v.zeros(start.n_rows, leaps_vector.at(ind) + 1);
        mat_t.col(0) = thetas_by_col.col(k);
        mat_v.col(0) = v;

        // --- Leapfrog integrator
        theta_new = thetas_by_col.col(k);
        v_new = v;

        // --- Energy of the current iteration
        inplace_tri_mat_mult<T_metric>(inv_chol_m, v);
        energy.at(k) = - log_pdf + 0.5 * arma::dot(v, v);

        try {
          for (arma::uword l = 1; l <= leaps_vector.at(ind); l++) {
            this->leapfrog_step(stepsizes_vector.at(ind), theta_new, v_new, grad_log_pdf_new);
            // --- Storing leapfrog intermediate point
            mat_t.col(l) = theta_new;
            mat_v.col(l) = v_new;
          }

          // --- Metropolis Hastings step
          inplace_tri_mat_mult<T_metric>(inv_chol_m, v_new);
          log_pdf_new = Rcpp::as<double>((*target_log_pdf)(theta_new));
          log_mh_prob.at(k) = energy.at(k) + log_pdf_new - 0.5 * arma::dot(v_new, v_new);

          if (log(arma::randu()) < log_mh_prob.at(k)) {
            // --- Accept the  move
            thetas_by_col.col(k + 1) = theta_new;
            log_pdf = log_pdf_new;
            grad_log_pdf = grad_log_pdf_new;

            // --- Update variable for the recycling step
            loc = leaps_vector.at(ind);
            sigma = 1;
          } else {
            // --- Reject the move
            thetas_by_col.col(k + 1) = thetas_by_col.col(k);
            moves.at(k) = 0;
            grad_log_pdf_new = grad_log_pdf;

            // --- Update variable for the recycling step
            loc = 0;
            sigma = -1;
          }
          lp.at(k + 1) = log_pdf;
        } catch (...) {
          // --- Reject the move
          thetas_by_col.col(k + 1) = thetas_by_col.col(k);
          log_mh_prob.at(k) = -std::numeric_limits<double>::infinity();
          moves.at(k) = 0;
          divergent.at(k) = 1;
          grad_log_pdf_new = grad_log_pdf;

          // --- Update variable for the recycling step
          loc = 0;
          sigma = -1;
          Rcpp::warning("Divergent transition occurs.");
        }



      } else {
        // --------------------------------------------------
        // --- The momentum is not refreshed
        // --------------------------------------------------

        bool has_diverged = false;

        if (random) {
          // --- Resample the number of leapfrog steps
          // --- But keep the same step size
          arma::uword ind_ = sample_int_1(weights);
          temp_l = leaps_vector.at(ind_);
          selected_ind.at(k) = ind_;
        } else {
          // --- Fixed L
          temp_l = leaps_vector.at(ind);
          selected_ind.at(k) = ind;
        }
        selected_eps.at(k) = stepsizes_vector.at(ind);

        // --- Use the momentum at the current location
        // --- 0 or the end point if it is the first time you do not refresh
        v = mat_v.col(loc);

        // --- Energy of the current iteration
        inplace_tri_mat_mult<T_metric>(inv_chol_m, v);
        energy.at(k) = - log_pdf + 0.5 * arma::dot(v, v);

        // --- Checking if we need to compute new points
        // --- lag is the number of new points to compute
        if (sigma > 0) {
          lag = loc + temp_l - mat_t.n_cols + 1;
        } else {
          lag = temp_l - loc;
        }


        if (lag > 0) {
          selected_l.at(k) = lag;
          // --- Compute new points along the level set

          // --- Resizing containers for the points traced out by leapfrog
          temp_t.zeros(start.n_rows, lag);
          temp_v.zeros(start.n_rows, lag);

          // --- Computing and adding new points to the leapfrog path
          if (sigma > 0) {
            // --- Starting point is the last point of
            // --- the excursion in the positive direction
            theta_new = mat_t.col(mat_t.n_cols - 1);
            v_new = mat_v.col(mat_v.n_cols - 1);
          } else {
            // --- Starting point is the first point of
            // --- the excursion in the negative direction
            theta_new = mat_t.col(0);
            v_new = -1.0 * mat_v.col(0);
          }

          // --- Add new points to the level set
          try {
            grad_log_pdf_new = Rcpp::as<arma::vec>((*target_grad_log_pdf)(theta_new));
            for (arma::uword l = 0; l < lag; l++) {
              this->leapfrog_step(stepsizes_vector.at(ind), theta_new, v_new, grad_log_pdf_new);
              temp_t.col(l) = theta_new;
              temp_v.col(l) = v_new;
            }

            // --- Concatenation of the new and old points
            if (sigma > 0) {
              shift = 0;
              mat_t = arma::join_rows(mat_t, temp_t);
              mat_v = arma::join_rows(mat_v, temp_v);
            } else {
              shift = lag;

              // --- minus one term to retrace the dynamic in the positive direction
              mat_t = arma::join_rows(temp_t, mat_t);
              mat_v = arma::join_rows(-1.0 * temp_v, mat_v);
            }
          } catch(...) {
            has_diverged = true;
            Rcpp::warning("Divergent transition occurs.");
          }
        } else {
          // --- Go to the existing point
          shift = 0.0;
          theta_new = mat_t.col(loc + sigma * temp_l);
          v_new = mat_v.col(loc + sigma * temp_l);
        }

        if (has_diverged) {
          thetas_by_col.col(k + 1) = thetas_by_col.col(k);
          log_mh_prob.at(k) = -std::numeric_limits<double>::infinity();
          moves.at(k) = 0;
          divergent.at(k) = 1;
          grad_log_pdf_new = grad_log_pdf;

          // --- Update variable for the recycling step
          sigma = -sigma;
          // DO not update loc since we have not added new points
        } else {
          // --- Metropolis Hastings step
          inplace_tri_mat_mult<T_metric>(inv_chol_m, v_new);
          log_pdf_new = Rcpp::as<double>((*target_log_pdf)(theta_new));
          log_mh_prob.at(k) = energy.at(k) + log_pdf_new - 0.5 * arma::dot(v_new, v_new);

          if (log(arma::randu()) < log_mh_prob.at(k)) {
            // --- Accept the  move
            thetas_by_col.col(k + 1) = theta_new;
            log_pdf = log_pdf_new;
            grad_log_pdf = grad_log_pdf_new;

            // --- Update variable for the recycling step
            loc = loc + sigma * temp_l + shift;
          } else {
            // --- Reject the move
            thetas_by_col.col(k + 1) = thetas_by_col.col(k);
            moves.at(k) = 0;
            grad_log_pdf_new = grad_log_pdf;

            // --- Update variable for the recycling step
            loc += shift;
            sigma = -sigma;
          }
        }
        lp.at(k + 1) = log_pdf;
      }
    }
    std::cout << " DONE!" << std::endl;
  };

  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------
  // ------------------------------------------------------

  void run(size_t iter,
           const arma::vec& start,
           const arma::Col<size_t>& leaps_vector,
           const arma::vec& stepsizes_vector,
           const arma::vec& weights,
           arma::mat& thetas_by_col,
           arma::vec& lp,
           arma::vec& energy,
           arma::vec& selected_eps,
           arma::Col<size_t>& selected_l,
           arma::uvec& selected_ind,
           arma::uvec& moves,
           arma::vec& log_mh_prob,
           arma::uvec& divergent,
           double refresh_rate = 1.0,
           bool random = true) {
    // --- Wrapper for the two versions of HMC
    if (refresh_rate < 1.0) {
      this->pr_ehmc(iter, start, leaps_vector, stepsizes_vector, weights, refresh_rate, random,
                    thetas_by_col, lp, energy, selected_eps, selected_l, selected_ind, moves, log_mh_prob, divergent);
    } else {
      this->ehmc(iter, start, leaps_vector, stepsizes_vector, weights,
                 thetas_by_col, lp, energy, selected_eps, selected_l, selected_ind, moves, log_mh_prob, divergent);
    }
  }

public:
  // --- Target distribution
  const Rcpp::Function *target_log_pdf;
  const Rcpp::Function *target_grad_log_pdf;
  // --- Mass Matrix
  arma::mat chol_m;
  arma::mat inv_chol_m;
  arma::mat inv_m;
  // -- Control struct
  AdaptControl adapt_control;
};


#endif /* INST_INCLUDE_ALGO_HPP_ */
