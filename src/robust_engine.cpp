#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Internal helper (not exported to R): numerically safe symmetric inverse.
// Tries the fast Cholesky-based inv_sympd() first (appropriate for the
// well-conditioned, positive-definite matrices this engine works with) and
// only falls back to the much slower SVD-based pinv() if that fails (e.g.
// due to near-singularity). This avoids paying the pinv() cost on every
// single EM iteration / profile / random start, which otherwise dominates
// runtime for larger datasets.
// ---------------------------------------------------------------------------
static arma::mat safe_inv_sympd(const arma::mat& Sigma) {
  arma::mat Sigma_inv;
  bool ok = arma::inv_sympd(Sigma_inv, Sigma);
  if (!ok) {
    Sigma_inv = arma::pinv(Sigma);
  }
  return Sigma_inv;
}

// Trimmed Centroid Based on Distance to the Coordinate-Wise Median
 //
 // Low-level utility that computes a simple robust location estimate: an
 // observation is included in the average only if its Euclidean distance to
 // the coordinate-wise median of the data is below \code{threshold}. This is
 // a one-step, easy-to-reason-about trimmed mean used internally as a cheap
 // fallback initializer; it is not part of the main robust EM/MCMC pipeline.
 //
 // @param X A numeric matrix.
 // @param threshold Maximum allowed Euclidean distance to the coordinate-wise
 //   median for an observation to be included in the average.
 // @return A numeric row vector with the trimmed centroid. Returns a vector
 //   of zeros (with a warning) if no observation falls within
 //   \code{threshold} of the median.
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::rowvec robust_mean_cpp(arma::mat X, double threshold) {
   int n = X.n_rows;
   int p = X.n_cols;
   
   // Coordinate-wise median as the reference center for trimming.
   arma::rowvec center(p);
   for (int j = 0; j < p; j++) {
     center(j) = arma::median(X.col(j));
   }
   
   arma::rowvec sum = arma::zeros<arma::rowvec>(p);
   int valid_count = 0;
   
   for (int i = 0; i < n; i++) {
     arma::rowvec current_row = X.row(i);
     double dist_to_center = arma::norm(current_row - center);
     if (dist_to_center <= threshold) {
       sum += current_row;
       valid_count++;
     }
   }
   
   if (valid_count > 0) {
     return sum / valid_count;
   } else {
     Rcpp::warning("robust_mean_cpp: no observation within `threshold` of the median; returning a zero vector.");
     return arma::zeros<arma::rowvec>(p);
   }
 }

// Squared Mahalanobis Distances (Complete Data)
 //
 // @param X A numeric matrix of complete (no missing values) observations.
 // @param mu A numeric row vector, the location parameter.
 // @param Sigma A numeric positive-(semi)definite covariance matrix.
 // @return A numeric vector of squared Mahalanobis distances, one per row of \code{X}.
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::vec mahalanobis_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
   int n = X.n_rows;
   arma::vec dists(n);
   arma::mat Sigma_inv = safe_inv_sympd(Sigma);
   
   for(int i = 0; i < n; i++) {
     arma::rowvec diff = X.row(i) - mu;
     arma::mat dist_sq = diff * Sigma_inv * diff.t();
     dists(i) = dist_sq(0, 0);
   }
   
   return dists;
 }

// Huber Weights from Squared Mahalanobis Distances (Complete Data)
 //
 // Observations whose squared Mahalanobis distance exceeds the chi-squared
 // cutoff are down-weighted proportionally to \code{sqrt(cutoff / dist)};
 // all other observations receive full weight 1.
 //
 // @param squared_dists Numeric vector of squared Mahalanobis distances.
 // @param chi_sq_cutoff The chi-squared quantile used as the outlier threshold.
 // @return A numeric vector of weights in (0, 1].
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::vec huber_weights_cpp(arma::vec squared_dists, double chi_sq_cutoff) {
   int n = squared_dists.n_elem;
   arma::vec weights = arma::ones<arma::vec>(n);
   
   for(int i = 0; i < n; i++) {
     if(squared_dists(i) > chi_sq_cutoff) {
       weights(i) = std::sqrt(chi_sq_cutoff / squared_dists(i));
     }
   }
   
   return weights;
 }

// Weighted (Robust, LASSO-Penalized) M-Step Update (Complete Data)
 //
 // Computes the Huber- and posterior-probability-weighted mean and
 // covariance matrix for a single profile, with optional soft-thresholding
 // (LASSO-type shrinkage) applied to the mean vector.
 //
 // @param X A numeric matrix of complete (no missing values) observations.
 // @param z Posterior membership probabilities for this profile (length \code{nrow(X)}).
 // @param w Huber robustness weights (length \code{nrow(X)}), as returned by \code{huber_weights_cpp()}.
 // @param lambda Non-negative soft-thresholding penalty. Shrinks each mean
 //   component toward zero; use only on centered/scaled data (see
 //   \code{\link{robust_lpa}} for details), otherwise shrinkage toward the
 //   raw-scale origin is not meaningful.
 // @return A list with \code{mean} (numeric row vector) and \code{covariance} (matrix).
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 Rcpp::List robust_update_cpp(arma::mat X, arma::vec z, arma::vec w, double lambda) {
   int n = X.n_rows;
   int p = X.n_cols;
   
   arma::vec v = z % w;
   double sum_v = arma::sum(v);
   
   arma::rowvec new_mu = arma::zeros<arma::rowvec>(p);
   for(int i = 0; i < n; i++) {
     new_mu += v(i) * X.row(i);
   }
   new_mu /= sum_v;
   
   if (lambda > 0.0) {
     for (int j = 0; j < p; ++j) {
       if (new_mu(j) > lambda) {
         new_mu(j) -= lambda;
       } else if (new_mu(j) < -lambda) {
         new_mu(j) += lambda;
       } else {
         new_mu(j) = 0.0;
       }
     }
   }
   
   arma::mat new_sigma = arma::zeros<arma::mat>(p, p);
   for(int i = 0; i < n; i++) {
     arma::rowvec diff = X.row(i) - new_mu;
     new_sigma += v(i) * (diff.t() * diff);
   }
   new_sigma /= sum_v;
   
   return Rcpp::List::create(
     Rcpp::Named("mean") = new_mu,
     Rcpp::Named("covariance") = new_sigma
   );
 }

// Multivariate Normal Density (Complete Data)
 //
 // @param X A numeric matrix of complete (no missing values) observations.
 // @param mu A numeric row vector, the mean.
 // @param Sigma A numeric positive-definite covariance matrix.
 // @return A numeric vector of density values, one per row of \code{X}.
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::vec dmvnorm_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
   int n = X.n_rows;
   int p = X.n_cols;
   arma::vec densities(n);
   arma::mat Sigma_inv = safe_inv_sympd(Sigma);
   
   // log_det() on a matrix that is not (numerically) positive definite can
   // return a complex/garbage exponent; Sigma is expected to already be
   // regularized to PD upstream (see force_pd() in robust_lpa.R), but we
   // guard the determinant sign here too so a bad Sigma yields a density of
   // 0 rather than propagating NaNs into the log-likelihood.
   double val;
   double sign;
   arma::log_det(val, sign, Sigma);
   if (sign <= 0) {
     return arma::zeros<arma::vec>(n);
   }
   double det_Sigma = std::exp(val);
   double constant = std::pow(2.0 * M_PI, -p / 2.0) * std::pow(det_Sigma, -0.5);
   
   for(int i = 0; i < n; i++) {
     arma::rowvec diff = X.row(i) - mu;
     arma::mat dist_sq = diff * Sigma_inv * diff.t();
     densities(i) = constant * std::exp(-0.5 * dist_sq(0, 0));
   }
   
   return densities;
 }

// Multivariate Normal Density Under Missingness (FIML)
 //
 // For each row, evaluates the multivariate normal density on the observed
 // subset of variables only (Full Information Maximum Likelihood), skipping
 // entries that are \code{NA}.
 //
 // @param X A numeric matrix, possibly containing \code{NA} values.
 // @param mu A numeric row vector, the mean (full dimension).
 // @param Sigma A numeric positive-definite covariance matrix (full dimension).
 // @return A numeric vector of density values, one per row of \code{X}. Rows
 //   with no observed variables get a density floored at \code{1e-300}.
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::vec dmvnorm_fiml_cpp(arma::mat X, arma::rowvec mu, arma::mat Sigma) {
   int n = X.n_rows;
   arma::vec densities(n);
   arma::vec mu_vec = mu.t();
   
   for(int i = 0; i < n; i++) {
     arma::rowvec Xi_row = X.row(i);
     arma::uvec obs_idx = arma::find_finite(Xi_row);
     int p_obs = obs_idx.n_elem;
     
     if(p_obs == 0) {
       densities(i) = 1e-300;
       continue;
     }
     
     arma::vec Xi_obs = Xi_row.elem(obs_idx);
     arma::vec mu_obs = mu_vec.elem(obs_idx);
     arma::mat Sigma_obs = Sigma.submat(obs_idx, obs_idx);
     
     Sigma_obs.diag() += 1e-6;
     arma::mat Sigma_inv = safe_inv_sympd(Sigma_obs);
     
     double val; double sign;
     arma::log_det(val, sign, Sigma_obs);
     if (sign <= 0) {
       densities(i) = 1e-300;
       continue;
     }
     
     arma::vec diff = Xi_obs - mu_obs;
     double exponent = arma::as_scalar(diff.t() * Sigma_inv * diff);
     double log_density = -0.5 * p_obs * std::log(2.0 * M_PI) - 0.5 * val - 0.5 * exponent;
     
     double d = std::exp(log_density);
     densities(i) = (d < 1e-300) ? 1e-300 : d;
   }
   return densities;
 }

// Squared Mahalanobis Distances Under Missingness (FIML)
 //
 // For each row, computes the squared Mahalanobis distance using only the
 // observed subset of variables.
 //
 // @param X A numeric matrix, possibly containing \code{NA} values.
 // @param mu A numeric row vector, the location parameter (full dimension).
 // @param Sigma A numeric positive-definite covariance matrix (full dimension).
 // @return A numeric vector of squared Mahalanobis distances, one per row of
 //   \code{X}. Rows with no observed variables get distance 0.
 // @keywords internal
 // @noRd
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
     arma::mat Sigma_inv = safe_inv_sympd(Sigma_obs);
     
     arma::vec diff = Xi_obs - mu_obs;
     dists(i) = arma::as_scalar(diff.t() * Sigma_inv * diff);
   }
   return dists;
 }

// Huber Weights from Squared Mahalanobis Distances Under Missingness (FIML)
 //
 // Like \code{huber_weights_cpp()}, but the chi-squared cutoff is
 // computed per-row using the row-specific number of observed variables
 // (\code{p_obs}), since the null distribution of the squared Mahalanobis
 // distance depends on the observed dimensionality.
 //
 // @param X A numeric matrix, possibly containing \code{NA} values (only
 //   used to determine, per row, which variables are observed).
 // @param squared_dists Numeric vector of squared Mahalanobis distances, as
 //   returned by \code{mahalanobis_fiml_cpp()}.
 // @param alpha Significance level used for the chi-squared cutoff (default \code{0.05}).
 // @return A numeric vector of weights in (0, 1].
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 arma::vec huber_weights_fiml_cpp(arma::mat X, arma::vec squared_dists, double alpha = 0.05) {
   int n = X.n_rows;
   arma::vec w(n, arma::fill::ones);
   
   for(int i = 0; i < n; i++) {
     arma::uvec obs_idx = arma::find_finite(X.row(i));
     int p_obs = obs_idx.n_elem;
     if(p_obs == 0) continue;
     
     double cutoff = R::qchisq(1.0 - alpha, p_obs, 1, 0);
     if(squared_dists(i) > cutoff) {
       w(i) = std::sqrt(cutoff / squared_dists(i));
     }
   }
   return w;
 }

// Weighted (Robust, LASSO-Penalized) M-Step Update Under Missingness (FIML)
 //
 // Like \code{robust_update_cpp()}, but the mean and covariance are
 // accumulated only from observed entries for each variable (and each pair
 // of variables, for the covariance), i.e. pairwise-available-case FIML
 // estimation, combined with posterior probability weights \code{z} and
 // Huber robustness weights \code{w}.
 //
 // @param X A numeric matrix, possibly containing \code{NA} values.
 // @param z Posterior membership probabilities for this profile (length \code{nrow(X)}).
 // @param w Huber robustness weights (length \code{nrow(X)}), as returned by \code{huber_weights_fiml_cpp()}.
 // @param lambda Non-negative soft-thresholding penalty applied to the mean vector (see \code{robust_update_cpp()}).
 // @return A list with \code{mean} (numeric row vector) and \code{covariance} (matrix).
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 Rcpp::List robust_update_fiml_cpp(arma::mat X, arma::vec z, arma::vec w, double lambda) {
   int n = X.n_rows;
   int p = X.n_cols;
   
   arma::rowvec new_mu(p, arma::fill::zeros);
   arma::mat new_sigma(p, p, arma::fill::zeros);
   
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
   
   if (lambda > 0.0) {
     for (int j = 0; j < p; ++j) {
       if (new_mu(j) > lambda) {
         new_mu(j) -= lambda;
       } else if (new_mu(j) < -lambda) {
         new_mu(j) += lambda;
       } else {
         new_mu(j) = 0.0;
       }
     }
   }
   
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
         new_sigma(k, j) = new_sigma(j, k);
       }
     }
   }
   
   return Rcpp::List::create(
     Rcpp::Named("mean") = new_mu,
     Rcpp::Named("covariance") = new_sigma
   );
 }

// Internal MCMC helpers (not exported to R): categorical sampling,
// inverse-Gaussian, inverse-Wishart (via Bartlett decomposition), and
// inverse-Gamma draws, plus a numerically stable log-sum-exp used for the
// allocation step. None of these are user-facing and none need Roxygen
// documentation.
int sample_class(arma::vec probs) {
  double u = R::runif(0.0, 1.0);
  double cumsum = 0.0;
  for(arma::uword i = 0; i < probs.n_elem; i++) {
    cumsum += probs(i);
    if(u <= cumsum) return i;
  }
  return probs.n_elem - 1;
}

double rinvgauss(double mu, double lambda) {
  double nu = R::rnorm(0.0, 1.0);
  double y = nu * nu;
  double x = mu + (mu * mu * y) / (2.0 * lambda) - (mu / (2.0 * lambda)) * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  double u = R::runif(0.0, 1.0);
  if (u <= mu / (mu + x)) {
    return x;
  } else {
    return (mu * mu) / x;
  }
}

// NOTE: `df` is a double (not int) because, under robust = TRUE, the
// degrees of freedom used here is an *effective* (Huber-weighted) sample
// size -- see the "eff_counts" accumulation in robust_mcmc_cpp() below --
// which is generally non-integer. R::rchisq() accepts a real-valued degrees
// of freedom, so this is well-defined for any df > 0.
arma::mat rinvwishart(double df, arma::mat S) {
  int p = S.n_rows;
  arma::mat Z(p, p, arma::fill::randn);
  arma::mat L_S = arma::chol(S, "lower");
  arma::mat L_Z = arma::trimatl(Z);

  for(int i = 0; i < p; i++) {
    L_Z(i, i) = std::sqrt(R::rchisq(df - i));
  }

  arma::mat X = L_S * arma::inv(L_Z.t());
  return X * X.t();
}

double rinvgamma_cpp(double shape, double rate) {
  return 1.0 / R::rgamma(shape, 1.0 / rate);
}

double log_sum_exp_cpp(arma::vec x) {
  double max_x = arma::max(x);
  return max_x + std::log(arma::sum(arma::exp(x - max_x)));
}

// Gibbs Sampler for Bayesian (Laplace-Shrinkage) Latent Profile Analysis
 //
 // Runs a single MCMC chain for a Gaussian mixture model with a Bayesian
 // Lasso (Laplace) prior on the profile means, supporting missing data via
 // row-wise available-case handling in the allocation step and in the
 // sufficient statistics. Called internally by \code{\link{robust_lpa}} once
 // per chain when \code{engine = "MCMC"} (see \code{n_chains} there); not
 // intended to be called directly by end users.
 //
 // When \code{robust = TRUE}, every sweep computes a Huber weight for each
 // observation from its squared Mahalanobis distance to its *currently
 // assigned* profile's previous-sweep mean/covariance (using the same
 // \code{1 - alpha} chi-squared cutoff as the EM engine's M-step; see
 // huber_weights_fiml_cpp() above). These weights down-weight each
 // observation's contribution to the sufficient statistics (weighted mean
 // and scatter matrix) used to draw that sweep's new mean/covariance, and
 // the resulting *effective* (weighted) per-profile sample size replaces the
 // raw allocation count in the conjugate updates' precision/degrees-of-
 // freedom terms. This mirrors, in a Gibbs-sampling setting, exactly what
 // robust_m_step()/robust_update_cpp() already do for the EM engine's
 // M-step: the allocation step itself (the Bayesian analogue of the E-step)
 // is not down-weighted, and neither is the mixing-proportion update
 // (pi_g), which still uses the raw allocation counts -- only the
 // sufficient statistics feeding the mean/covariance draws are robustified.
 // When \code{robust = FALSE}, all weights are fixed at 1 and this reduces
 // exactly to the original (non-robust) sampler.
 //
 // @param X A numeric matrix, possibly containing \code{NA} values.
 // @param G Integer, the number of latent profiles.
 // @param model Integer from 1 to 6, the variance-covariance parameterization
 //   (see \code{\link{robust_lpa}} for the full description of each model).
 // @param mcmc_iter Integer, the total number of MCMC iterations (single chain).
 // @param prior_laplace Positive numeric, the Laplace (Bayesian Lasso) prior scale for the means.
 // @param robust Logical, defaults to \code{true}. If \code{true}, Huber
 //   down-weighting of outlying observations is applied as described above.
 //   If \code{false}, classical (non-robust) Gibbs updates are used.
 // @param alpha Significance level for the Huber down-weighting threshold,
 //   defaults to \code{0.05}. Ignored when \code{robust = false}.
 // @return A list with \code{mu_chain} (list of length \code{mcmc_iter}, each
 //   a list of \code{G} mean vectors), \code{sigma_chain} (list of length
 //   \code{mcmc_iter}, each a list of \code{G} covariance matrices), and
 //   \code{pi_chain} (a \code{mcmc_iter x G} matrix of mixing proportions).
 // @keywords internal
 // @noRd
 // [[Rcpp::export]]
 Rcpp::List robust_mcmc_cpp(arma::mat X, int G, int model, int mcmc_iter, double prior_laplace,
                             bool robust = true, double alpha = 0.05) {
   int n = X.n_rows;
   int p = X.n_cols;

   arma::vec z_alloc(n, arma::fill::zeros);
   arma::vec pi_g(G, arma::fill::value(1.0 / G));
   // Per-observation Huber robustness weight, recomputed every sweep from
   // the previous sweep's parameters (see huber_weights_fiml_cpp() for the
   // complete-data-agnostic, per-row chi-squared cutoff this mirrors).
   // Initialized to 1 (no down-weighting) before the first allocation.
   arma::vec obs_weight(n, arma::fill::ones);

   Rcpp::List mu(G);
   Rcpp::List sigma(G);
   Rcpp::List tau(G);

   for(int g = 0; g < G; g++) {
     mu[g] = arma::zeros<arma::rowvec>(p);
     sigma[g] = arma::eye<arma::mat>(p, p);
     tau[g] = arma::ones<arma::vec>(p);
   }

   Rcpp::List mu_chain(mcmc_iter);
   Rcpp::List sigma_chain(mcmc_iter);
   arma::mat pi_chain(mcmc_iter, G);

   for(int iter = 0; iter < mcmc_iter; iter++) {

     for(int i = 0; i < n; i++) {
       arma::vec log_probs(G, arma::fill::zeros);
       arma::vec dist_sq_by_g(G, arma::fill::value(arma::datum::inf));
       arma::rowvec Xi_row = X.row(i);
       arma::uvec obs_idx = arma::find_finite(Xi_row);
       int p_obs = obs_idx.n_elem;

       if(p_obs > 0) {
         arma::vec Xi_obs = Xi_row.elem(obs_idx);

         for(int g = 0; g < G; g++) {
           arma::rowvec current_mu_full = mu[g];
           arma::mat current_sigma_full = sigma[g];

           arma::vec mu_obs = current_mu_full.elem(obs_idx);
           arma::mat Sigma_obs = current_sigma_full.submat(obs_idx, obs_idx);
           Sigma_obs.diag() += 1e-6;

           arma::mat Sigma_inv = safe_inv_sympd(Sigma_obs);
           double val;
           double sign;
           arma::log_det(val, sign, Sigma_obs);
           if (sign <= 0) {
             log_probs(g) = -arma::datum::inf;
             continue;
           }

           arma::vec diff = Xi_obs - mu_obs;
           double dist_sq = arma::as_scalar(diff.t() * Sigma_inv * diff);
           dist_sq_by_g(g) = dist_sq;

           log_probs(g) = std::log(pi_g(g)) - 0.5 * dist_sq - 0.5 * val;
         }

         arma::vec probs = arma::exp(log_probs - log_sum_exp_cpp(log_probs));
         int g_star = sample_class(probs);
         z_alloc(i) = g_star;

         if (robust) {
           double cutoff = R::qchisq(1.0 - alpha, p_obs, 1, 0);
           double dsq = dist_sq_by_g(g_star);
           obs_weight(i) = (dsq > cutoff) ? std::sqrt(cutoff / dsq) : 1.0;
         } else {
           obs_weight(i) = 1.0;
         }
       } else {
         // No observed variables for this row: it cannot inform any
         // profile's sufficient statistics either way (every X(i,j) term
         // below is skipped by the std::isfinite() guards regardless of
         // group), but it must still be *allocated* somewhere every sweep.
         // Its allocation posterior given zero data is exactly the current
         // prior pi_g (Bayes' rule with an uninformative likelihood), so
         // sample from that directly. (A previous version of this sampler
         // left z_alloc(i) frozen at its initial value -- always profile 1
         // -- forever for such rows, which silently inflated that
         // profile's allocation count/pi_g whenever a fully-missing row
         // was present.)
         z_alloc(i) = sample_class(pi_g);
         obs_weight(i) = 1.0;
       }
     }

     Rcpp::List raw_S(G);
     arma::vec n_counts(G, arma::fill::zeros);
     // Effective (Huber-weighted) per-profile sample size: the sum of
     // obs_weight() over observations currently allocated to profile g. With
     // robust = FALSE (all weights 1) this equals n_counts(g) exactly, so
     // every formula below reduces to the original (non-robust) sampler.
     arma::vec eff_counts(G, arma::fill::zeros);

     for(int g = 0; g < G; g++) {
       arma::uvec indices = arma::find(z_alloc == g);
       int n_g = indices.n_elem;
       n_counts(g) = n_g;

       arma::mat S_g(p, p, arma::fill::zeros);
       double eff_n_g = 0.0;

       if(n_g > 0) {
         arma::rowvec x_bar(p, arma::fill::zeros);
         arma::vec valid_weight(p, arma::fill::zeros);

         for(int i = 0; i < n_g; i++) {
           int global_idx = indices(i);
           double w_i = obs_weight(global_idx);
           for(int j = 0; j < p; j++) {
             if(std::isfinite(X(global_idx, j))) {
               x_bar(j) += w_i * X(global_idx, j);
               valid_weight(j) += w_i;
             }
           }
         }

         for(int j = 0; j < p; j++) {
           if(valid_weight(j) > 1e-8) x_bar(j) /= valid_weight(j);
         }

         for(int i = 0; i < n_g; i++) {
           eff_n_g += obs_weight(indices(i));
         }

         arma::vec current_tau = tau[g];
         arma::mat current_sigma = sigma[g];
         arma::mat Sigma_inv = safe_inv_sympd(current_sigma);

         // D_tau is diagonal by construction (independent Laplace scale
         // mixture per coordinate), so its inverse is just the elementwise
         // reciprocal -- no need for a general (and much slower) matrix
         // inverse here.
         arma::vec D_tau_inv_diag = 1.0 / arma::clamp(current_tau, 1e-8, arma::datum::inf);
         arma::mat V_post = safe_inv_sympd(eff_n_g * Sigma_inv + arma::diagmat(D_tau_inv_diag));
         arma::rowvec mu_post = x_bar * (eff_n_g * Sigma_inv) * V_post;

         arma::mat chol_V = arma::chol(V_post, "lower");
         arma::vec Z_norm(p, arma::fill::randn);
         arma::rowvec new_mu = mu_post + (chol_V * Z_norm).t();
         mu[g] = new_mu;

         for(int j = 0; j < p; j++) {
           // Always resample tau_j (the Bayesian Lasso auxiliary prior
           // variance for mu_j) from its full conditional, every sweep --
           // required for a valid Gibbs sampler. A previous version of this
           // sampler skipped this update whenever |mu_j| <= 1e-6 (to avoid
           // dividing by an exactly-zero mu_j), which silently froze tau_j
           // at its previous value precisely in the sparse regime the
           // Laplace prior exists to explore, biasing the chain away from
           // its intended stationary distribution. Flooring |mu_j| at 1e-6
           // avoids the division by zero without skipping the update: as
           // mu_j -> 0 the true full conditional's mean parameter mu' ->
           // Inf, so this floor only matters in that (measure-zero) limit.
           double mu_j_abs = std::max(std::abs(new_mu(j)), 1e-6);
           double mu_prime = prior_laplace / mu_j_abs;
           current_tau(j) = 1.0 / rinvgauss(mu_prime, prior_laplace * prior_laplace);
         }
         tau[g] = current_tau;

         for(int i = 0; i < n_g; i++) {
           int global_idx = indices(i);
           double w_i = obs_weight(global_idx);
           for(int j = 0; j < p; j++) {
             for(int k = 0; k <= j; k++) {
               if(std::isfinite(X(global_idx, j)) && std::isfinite(X(global_idx, k))) {
                 double diff_j = X(global_idx, j) - new_mu(j);
                 double diff_k = X(global_idx, k) - new_mu(k);
                 S_g(j, k) += w_i * diff_j * diff_k;
                 if(j != k) S_g(k, j) = S_g(j, k);
               }
             }
           }
         }
       }
       eff_counts(g) = eff_n_g;
       raw_S[g] = S_g;
     }

     // Draw the mixing proportions pi_g from their Dirichlet(n_1+1, ..., n_G+1)
     // full conditional -- the conjugate posterior for a Dirichlet(1,...,1)
     // prior given multinomial allocation counts n_counts -- via G independent
     // Gamma(shape = n_g + 1, rate = 1) draws normalized to sum to 1 (the
     // standard construction of a Dirichlet variate). The mixing-proportion
     // update intentionally still uses the raw allocation count (not
     // eff_counts), mirroring the EM engine, where pi_g is updated from the
     // (non-robustified) posterior responsibilities z.
     //
     // A previous version of this sampler instead *fixed* pi_g at this
     // Dirichlet's posterior mean, (n_g + 1) / (n + G), every sweep -- a valid
     // point estimate but not a valid Gibbs draw. Whenever the hard allocation
     // counts happened to repeat from one sweep to the next (which they
     // typically do once the chain has settled into a stable classification,
     // e.g. a well-separated fit with entropy near 1), pi_g was silently
     // *exactly constant* across those sweeps: its true posterior uncertainty
     // was never sampled at all, and in the extreme this produced a fully
     // degenerate (zero-variance) pi_chain -- visible downstream as
     // Gelman-Rubin R-hat = NaN and effective sample size = 0 for pi in
     // $mcmc_diagnostics, even when every other parameter's chain looked
     // healthy.
     {
       arma::vec gamma_draws(G);
       for(int g = 0; g < G; g++) {
         gamma_draws(g) = R::rgamma(n_counts(g) + 1.0, 1.0);
       }
       pi_g = gamma_draws / arma::sum(gamma_draws);
     }

     if(model == 1) {
       double pooled_diag = 0;
       double total_eff_n = arma::sum(eff_counts);
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         pooled_diag += arma::accu(S.diag());
       }
       double var_val = rinvgamma_cpp(0.5 * total_eff_n * p + 1.0, 0.5 * pooled_diag + 1.0);
       for(int g = 0; g < G; g++) {
         sigma[g] = arma::eye<arma::mat>(p, p) * var_val;
       }

     } else if(model == 2) {
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         arma::mat new_sig(p, p, arma::fill::zeros);
         for(int j = 0; j < p; j++) {
           new_sig(j, j) = rinvgamma_cpp(0.5 * eff_counts(g) + 1.0, 0.5 * S(j, j) + 1.0);
         }
         sigma[g] = new_sig;
       }

     } else if(model == 3) {
       arma::mat S_pool(p, p, arma::fill::zeros);
       double total_eff_n = arma::sum(eff_counts);
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         S_pool += S;
       }
       arma::mat pooled_sig = rinvwishart(total_eff_n + p + 1.0, S_pool + arma::eye<arma::mat>(p, p));
       for(int g = 0; g < G; g++) {
         sigma[g] = pooled_sig;
       }

     } else if(model == 4) {
       arma::mat S_pool(p, p, arma::fill::zeros);
       double total_eff_n = arma::sum(eff_counts);
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         S_pool += S;
       }
       arma::mat pooled_sig = rinvwishart(total_eff_n + p + 1.0, S_pool + arma::eye<arma::mat>(p, p));

       arma::vec pooled_sd = arma::sqrt(pooled_sig.diag());
       arma::mat R_pool = pooled_sig;
       for(int j = 0; j < p; j++) {
         for(int k = 0; k < p; k++) {
           R_pool(j, k) /= (pooled_sd(j) * pooled_sd(k));
         }
       }

       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         arma::mat new_sig = R_pool;
         arma::vec new_sd(p);

         for(int j = 0; j < p; j++) {
           new_sd(j) = std::sqrt(rinvgamma_cpp(0.5 * eff_counts(g) + 1.0, 0.5 * S(j, j) + 1.0));
         }

         for(int j = 0; j < p; j++) {
           for(int k = 0; k < p; k++) {
             new_sig(j, k) *= (new_sd(j) * new_sd(k));
           }
         }
         sigma[g] = new_sig;
       }

     } else if(model == 5) {
       double pooled_diag = 0;
       double total_eff_n = arma::sum(eff_counts);
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         pooled_diag += arma::accu(S.diag());
       }

       double var_val = rinvgamma_cpp(0.5 * total_eff_n * p + 1.0, 0.5 * pooled_diag + 1.0);
       double sd_val = std::sqrt(var_val);

       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         arma::mat temp_sig = rinvwishart(eff_counts(g) + p + 1.0, S + arma::eye<arma::mat>(p, p));

         arma::vec temp_sd = arma::sqrt(temp_sig.diag());
         arma::mat R_g = temp_sig;

         for(int j = 0; j < p; j++) {
           for(int k = 0; k < p; k++) {
             R_g(j, k) /= (temp_sd(j) * temp_sd(k));
             R_g(j, k) *= (sd_val * sd_val);
           }
         }
         sigma[g] = R_g;
       }

     } else if(model == 6) {
       for(int g = 0; g < G; g++) {
         arma::mat S = raw_S[g];
         sigma[g] = rinvwishart(eff_counts(g) + p + 1.0, S + arma::eye<arma::mat>(p, p));
       }
     }

     pi_chain.row(iter) = pi_g.t();
     mu_chain[iter] = Rcpp::clone(mu);
     sigma_chain[iter] = Rcpp::clone(sigma);
   }

   return Rcpp::List::create(
     Rcpp::Named("mu_chain") = mu_chain,
     Rcpp::Named("sigma_chain") = sigma_chain,
     Rcpp::Named("pi_chain") = pi_chain
   );
 }
