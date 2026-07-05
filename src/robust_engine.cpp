#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec robust_mean_cpp(arma::mat X, double threshold) {

   int n = X.n_rows; // number of rows
   int p = X.n_cols; // number of columns

   arma::rowvec sum = arma::zeros<arma::rowvec>(p);
   int valid_count = 0;

   for(int i = 0; i < n; i++) {
     arma::rowvec current_row = X.row(i);

     // If the distance from the center is below the threshold, we keep it
     if(arma::norm(current_row) <= threshold) {
       sum += current_row;
       valid_count++;
     }
   }

   if(valid_count > 0) {
     return sum / valid_count;
   } else {
     return arma::zeros<arma::rowvec>(p);
   }
 }

// Calculates the squared Mahalanobis distance for a matrix of data
//
// [[Rcpp::export]]
arma::vec mahalanobis_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {

  int n = X.n_rows;
  arma::vec dists(n); // Column vector to store the results

  // Calculate the inverse of the covariance matrix
  // We use pinv (pseudo-inverse) for greater numerical stability (Robustness)
  arma::mat Sigma_inv = arma::pinv(Sigma);

  for(int i = 0; i < n; i++) {
    // 1. Calculate the difference between the point and the mean (x_i - mu)
    arma::rowvec diff = X.row(i) - mu;

    // 2. Matrix multiplication: diff * Sigma_inv * diff.transpose
    arma::mat dist_sq = diff * Sigma_inv * diff.t();

    // 3. Extract the resulting single number and save it
    dists(i) = dist_sq(0, 0);
  }

  return dists;
}

// Calculates Huber weights based on squared Mahalanobis distances
//
// [[Rcpp::export]]
arma::vec huber_weights_cpp(arma::vec squared_dists, double chi_sq_cutoff) {

  int n = squared_dists.n_elem;
  arma::vec weights = arma::ones<arma::vec>(n); // Initialize all weights to 1

  for(int i = 0; i < n; i++) {
    // If the squared distance exceeds the cutoff, down-weight the observation
    if(squared_dists(i) > chi_sq_cutoff) {
      weights(i) = std::sqrt(chi_sq_cutoff / squared_dists(i));
    }
  }

  return weights;
}

// Calculates the robust updated mean and covariance for a single cluster (M-step)
//
// [[Rcpp::export]]
Rcpp::List robust_update_cpp(arma::mat X, arma::vec z, arma::vec w) {

  int n = X.n_rows;
  int p = X.n_cols;

  // Element-wise multiplication of probabilities (z) and robust weights (w)
  arma::vec v = z % w;
  double sum_v = arma::sum(v);

  // 1. Calculate the robust weighted mean
  arma::rowvec new_mu = arma::zeros<arma::rowvec>(p);
  for(int i = 0; i < n; i++) {
    new_mu += v(i) * X.row(i);
  }
  new_mu /= sum_v;

  // 2. Calculate the robust weighted covariance
  arma::mat new_sigma = arma::zeros<arma::mat>(p, p);
  for(int i = 0; i < n; i++) {
    arma::rowvec diff = X.row(i) - new_mu;
    new_sigma += v(i) * (diff.t() * diff); // diff.t() * diff is a p x p matrix
  }
  new_sigma /= sum_v;

  // Return both results as an R list
  return Rcpp::List::create(
    Rcpp::Named("mean") = new_mu,
    Rcpp::Named("covariance") = new_sigma
  );
}

// Calculates the Multivariate Normal probability density for a matrix of points
//
// [[Rcpp::export]]
arma::vec dmvnorm_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
  int n = X.n_rows;
  int p = X.n_cols;
  arma::vec densities(n);

  // Inverse of the covariance matrix for the distance calculation
  arma::mat Sigma_inv = arma::pinv(Sigma);

  // Calculate the determinant of the covariance matrix
  // We use log_det for numerical stability, then exponentiate it
  double val;
  double sign;
  arma::log_det(val, sign, Sigma);
  double det_Sigma = std::exp(val);

  // The constant part of the density formula
  double constant = std::pow(2.0 * M_PI, -p / 2.0) * std::pow(det_Sigma, -0.5);

  for(int i = 0; i < n; i++) {
    arma::rowvec diff = X.row(i) - mu;
    arma::mat dist_sq = diff * Sigma_inv * diff.t();

    // Calculate the density for this specific point
    densities(i) = constant * std::exp(-0.5 * dist_sq(0, 0));
  }

  return densities;
}

// NEW FIML FUNCTIONS (FOR HIGH-SPEED HANDLING OF MISSING DATA)
// Calculates the multivariate normal density dynamically ignoring NAs
// [[Rcpp::export]]
arma::vec dmvnorm_fiml_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
  int n = X.n_rows;
  arma::vec densities(n);
  arma::vec mu_vec = mu.t(); // Transpose to facilitate extraction

  for(int i = 0; i < n; i++) {
    arma::rowvec Xi_row = X.row(i);
    arma::uvec obs_idx = arma::find_finite(Xi_row); // Find NON-missing columns
    int p_obs = obs_idx.n_elem;

    if(p_obs == 0) {
      densities(i) = 1e-300;
      continue;
    }

    // Dynamically extracts submatrices omitting NAs for the current row
    arma::vec Xi_obs = Xi_row.elem(obs_idx);
    arma::vec mu_obs = mu_vec.elem(obs_idx);
    arma::mat Sigma_obs = Sigma.submat(obs_idx, obs_idx);

    Sigma_obs.diag() += 1e-6; // Prevents singularities
    arma::mat Sigma_inv = arma::pinv(Sigma_obs);

    double val; double sign;
    arma::log_det(val, sign, Sigma_obs);

    arma::vec diff = Xi_obs - mu_obs;
    double exponent = arma::as_scalar(diff.t() * Sigma_inv * diff);
    double log_density = -0.5 * p_obs * std::log(2.0 * M_PI) - 0.5 * val - 0.5 * exponent;

    double d = std::exp(log_density);
    densities(i) = (d < 1e-300) ? 1e-300 : d;
  }
  return densities;
}

// Calculates the FIML Mahalanobis distance (Conditional to observed variables)
// [[Rcpp::export]]
arma::vec mahalanobis_fiml_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
  int n = X.n_rows;
  arma::vec dists(n, arma::fill::zeros);
  arma::vec mu_vec = mu.t();

  for(int i = 0; i < n; i++) {
    arma::rowvec Xi_row = X.row(i);
    arma::uvec obs_idx = arma::find_finite(Xi_row);
    int p_obs = obs_idx.n_elem;

    if(p_obs == 0) continue;

    arma::vec Xi_obs = Xi_row.elem(obs_idx);
    arma::vec mu_obs = mu_vec.elem(obs_idx);
    arma::mat Sigma_obs = Sigma.submat(obs_idx, obs_idx);

    Sigma_obs.diag() += 1e-6;
    arma::mat Sigma_inv = arma::pinv(Sigma_obs);

    arma::vec diff = Xi_obs - mu_obs;
    dists(i) = arma::as_scalar(diff.t() * Sigma_inv * diff);
  }
  return dists;
}

// Calculates Huber weights adjusting degrees of freedom for each case (FIML)
// [[Rcpp::export]]
arma::vec huber_weights_fiml_cpp(arma::mat X, arma::vec squared_dists, double alpha = 0.05) {
  int n = X.n_rows;
  arma::vec w(n, arma::fill::ones);

  for(int i = 0; i < n; i++) {
    arma::uvec obs_idx = arma::find_finite(X.row(i));
    int p_obs = obs_idx.n_elem;
    if(p_obs == 0) continue;

    // Dynamic Chi-Square calculated natively in C++
    double cutoff = R::qchisq(1.0 - alpha, p_obs, 1, 0);
    if(squared_dists(i) > cutoff) {
      w(i) = std::sqrt(cutoff / squared_dists(i));
    }
  }
  return w;
}

// Updates mean and covariance with Pairwise approach (Maximum C++ speed)
// [[Rcpp::export]]
Rcpp::List robust_update_fiml_cpp(arma::mat X, arma::vec z, arma::vec w) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::rowvec new_mu(p, arma::fill::zeros);
  arma::mat new_sigma(p, p, arma::fill::zeros);

  // Pairwise Mean Calculation
  for(int j = 0; j < p; j++) {
    double sum_val = 0;
    double sum_w = 0;
    for(int i = 0; i < n; i++) {
      if(std::isfinite(X(i, j))) {
        double weight = z(i) * w(i);
        sum_val += weight * X(i, j);
        sum_w += weight;
      }
    }
    new_mu(j) = (sum_w > 1e-6) ? (sum_val / sum_w) : 0;
  }

  // Pairwise Covariance Calculation
  for(int j = 0; j < p; j++) {
    for(int k = 0; k <= j; k++) {
      double sum_val = 0;
      double sum_w = 0;
      for(int i = 0; i < n; i++) {
        if(std::isfinite(X(i, j)) && std::isfinite(X(i, k))) {
          double weight = z(i) * w(i);
          sum_val += weight * (X(i, j) - new_mu(j)) * (X(i, k) - new_mu(k));
          sum_w += weight;
        }
      }
      if(sum_w > 1e-6) {
        new_sigma(j, k) = sum_val / sum_w;
        new_sigma(k, j) = new_sigma(j, k); // Symmetrization
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("mean") = new_mu,
    Rcpp::Named("covariance") = new_sigma
  );
}
