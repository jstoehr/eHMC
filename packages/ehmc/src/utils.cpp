// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <utils.hpp>


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

type_metric string_to_metric(const std::string& metric_str) {
  // --- Convert a string into a type_metric object
  if (metric_str == "unit") {
    return unit;
  } else if (metric_str == "diag") {
    return diag;
  } else if (metric_str == "dense") {
    return dense;
  } else {
    return unknown;
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

//[[Rcpp::export]]
std::string what_metric_type(const arma::mat& m) {
  // --- Determine the type of matrix used as mass matrix
  std::string metric_str;
  if (m.is_diagmat()) {
    if(arma::norm(m.diag() - arma::ones(m.n_rows)) > 100 * arma::datum::eps) {
      metric_str = "diag";
    } else {
      metric_str = "unit";
    }
  } else {
    metric_str = "dense";
  }
  return metric_str;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
inline arma::vec template_mat_mult<unit>(const arma::mat& mat, const arma::vec& x) {
  // --- Multiplication by identity matrix
  return x;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
inline arma::vec template_mat_mult<diag>(const arma::mat& mat, const arma::vec& x) {
  // --- Multiplication by diagonal matrix
  return mat.diag() % x;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void inplace_tri_mat_mult<unit>(const arma::mat& trimat,
                                arma::vec& x) {
  // --- Specialization for identity matrix
  // --- Nothing to do
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void inplace_tri_mat_mult<diag>(const arma::mat& trimat,
                                arma::vec& x) {
  // --- Specialization for diagonal matrix
  // --- Could be replaced by x %= trimat.diag();
  for (auto k = 0; k < trimat.n_rows; k++) {
    x.at(k) *= trimat.at(k, k);
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void sympd_version<unit>(const arma::mat& x,
                         arma::mat& out,
                         double reg_param) {
  // --- Specialization for identity matrix
  // --- Nothing to do
  out = x;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void sympd_version<diag>(const arma::mat& x,
                         arma::mat& out,
                         double reg_param) {
  // --- Regularization for diagonal matrix
  // --- Define small threshold based on precision
  double temp = 100 * arma::datum::eps * arma::norm(x.diag(), 2);

  if (x.diag().min() > temp) {
    // If all diagonal elements are large enough, just return the diagonal matrix
    out = arma::diagmat(x);
  } else {
    // --- For small diagonal elements, regularization is applied
    out = arma::zeros<arma::mat>(x.n_rows, x.n_cols);
    out.diag() = x.diag();
    out.diag().clean(temp); // Remove very small values (below threshold)

    // --- Apply regularization on the diagonal elements
    arma::vec temp_diag = x.diag();
    arma::uvec ind = find(temp_diag > temp); // Find elements above the threshold
    double d = arma::mean(temp_diag.elem(ind)); // Compute the mean of these elements
    out.diag() *= (1.0 - reg_param); // Scale the diagonal
    out.diag() += reg_param * d; // Add the regularized mean
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void sympd_version<dense>(const arma::mat& x,
                          arma::mat& out,
                          double reg_param) {
  // --- Step 1: enforce symmetry
  arma::mat x_sym = 0.5 * (x + x.t());
  // --- Regularization for dense matrices
  // --- Eigen values decomposition
  arma::vec eigval;
  arma::mat eigvec;
  bool success = eig_sym(eigval, eigvec, x_sym);

  if (!success) {
    Rcpp::stop("sympd_version<dense>: eig_sym failed.");
  }
  double scale = arma::norm(x_sym, "fro");
  double eps_floor = 100.0 * arma::datum::eps * scale;

  if (eigval.min() > eps_floor) {
    // If eigenvalues are large enough
    out = x_sym;
  } else {
    // --- Regularization step for small eigenvalues
    double mean_diag = arma::mean(x_sym.diag());
    out = (1.0 - reg_param) * x_sym;
    out.diag() += reg_param * mean_diag;

    eig_sym(eigval, eigvec, out);
    scale = std::max(arma::norm(out, "fro"), 1.0);
    eps_floor = 100.0 * arma::datum::eps * scale;

    for (arma::uword i = 0; i < eigval.n_elem; ++i) {
      if (eigval[i] < eps_floor) {
        eigval[i] = eps_floor;
      }
    }
    out = eigvec * arma::diagmat(eigval) * eigvec.t();
    out = 0.5 * (out + out.t());
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void get_matrix_decomposition<unit>(const arma::mat& temp_inv_m,
                                    arma::mat& inv_m,
                                    arma::mat& chol_m,
                                    arma::mat& inv_chol_m,
                                    double reg_param) {
  // --- Matrix decomposition for the identity matrix
  inv_m = arma::eye(temp_inv_m.n_rows, temp_inv_m.n_cols);
  chol_m = arma::eye(temp_inv_m.n_rows, temp_inv_m.n_cols);
  inv_chol_m = arma::eye(temp_inv_m.n_rows, temp_inv_m.n_cols);
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void get_matrix_decomposition<diag>(const arma::mat& temp_inv_m,
                                    arma::mat& inv_m,
                                    arma::mat& chol_m,
                                    arma::mat& inv_chol_m,
                                    double reg_param) {
  // --- Matrix decomposition for a diagonal matrix
  sympd_version<diag>(temp_inv_m, inv_m, reg_param);
  arma::vec temp = arma::sqrt(inv_m.diag());
  chol_m = arma::diagmat(1.0/temp);
  inv_chol_m = arma::diagmat(temp);
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

template<>
void get_matrix_decomposition<dense>(const arma::mat& temp_inv_m,
                                     arma::mat& inv_m,
                                     arma::mat& chol_m,
                                     arma::mat& inv_chol_m,
                                     double reg_param) {
  // --- Matrix decomposition for a dense matrix
  sympd_version<dense>(temp_inv_m, inv_m, reg_param);
  arma::mat m = arma::inv_sympd(inv_m);
  chol_m = arma::chol(m);
  inv_chol_m = arma::inv(trimatu(chol_m));
}
