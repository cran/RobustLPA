#ifndef RLPA_UTILS_H
#define RLPA_UTILS_H

#include <RcppArmadillo.h>

// Shared numerical and random-variate helpers for the LPA and growth-model
// engines (not exported to R).

const double LOG_2PI = std::log(2.0 * M_PI);

// Lower-triangular Cholesky factor L of a (symmetrized) covariance matrix,
// S = L L'. If S is not numerically positive definite, an escalating ridge
// (relative to the average variance) is added until the factorization
// succeeds; as a last resort a diagonal factor is returned, so callers never
// receive an invalid factor. Covariances reaching this engine have normally
// already been regularized upstream (.force_pd() in R/robust_lpa.R), so the
// fallbacks are a numerical safety net only.
inline arma::mat chol_lower_safe(const arma::mat& S) {
  arma::mat Ssym = 0.5 * (S + S.t());
  arma::mat L;
  if (arma::chol(L, Ssym, "lower")) return L;
  double scale = std::max(arma::mean(arma::abs(Ssym.diag())), 1e-8);
  double jit = 1e-10 * scale;
  for (int k = 0; k < 14; ++k) {
    arma::mat Sj = Ssym;
    Sj.diag() += jit;
    if (arma::chol(L, Sj, "lower")) return L;
    jit *= 10.0;
  }
  return arma::diagmat(arma::sqrt(arma::clamp(Ssym.diag(), 1e-8, arma::datum::inf)));
}

// Symmetric positive-definite inverse via the safe Cholesky factor above.
inline arma::mat inv_spd_safe(const arma::mat& S) {
  arma::mat L = chol_lower_safe(S);
  arma::mat Linv = arma::inv(arma::trimatl(L));
  return Linv.t() * Linv;
}

// Log-density constant for a p_o-dimensional Gaussian (dist = 0) or
// multivariate t with `nu` degrees of freedom (dist = 1), given log|Sigma_oo|.
inline double log_const(int dist, int p_o, double logdet, double nu) {
  if (dist == 1) {
    return std::lgamma(0.5 * (nu + p_o)) - std::lgamma(0.5 * nu) -
      0.5 * p_o * std::log(nu * M_PI) - 0.5 * logdet;
  }
  return -0.5 * (p_o * LOG_2PI + logdet);
}

inline double log_kernel(int dist, int p_o, double maha, double nu) {
  if (dist == 1) return -0.5 * (nu + p_o) * std::log1p(maha / nu);
  return -0.5 * maha;
}

inline int sample_class(const arma::vec& probs) {
  double u = R::runif(0.0, 1.0);
  double cumsum = 0.0;
  for (arma::uword i = 0; i < probs.n_elem; i++) {
    cumsum += probs(i);
    if (u <= cumsum) return i;
  }
  return probs.n_elem - 1;
}

inline double rinvgauss(double mu, double lambda) {
  double nu = R::rnorm(0.0, 1.0);
  double y = nu * nu;
  double x = mu + (mu * mu * y) / (2.0 * lambda) -
    (mu / (2.0 * lambda)) * std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y);
  double u = R::runif(0.0, 1.0);
  return (u <= mu / (mu + x)) ? x : (mu * mu) / x;
}

// Inverse-Wishart draw via the Bartlett decomposition. `df` is real-valued
// (it is an effective, Huber-weighted sample size under the heuristic
// robust sampler); R::rchisq() accepts non-integer degrees of freedom.
inline arma::mat rinvwishart(double df, const arma::mat& S) {
  int p = S.n_rows;
  arma::mat Z(p, p, arma::fill::randn);
  arma::mat L_S = chol_lower_safe(S);
  arma::mat L_Z = arma::trimatl(Z);
  for (int i = 0; i < p; i++) L_Z(i, i) = std::sqrt(R::rchisq(df - i));
  arma::mat Xm = L_S * arma::inv(arma::trimatu(L_Z.t()));
  return Xm * Xm.t();
}

inline double rinvgamma_cpp(double shape, double rate) {
  return 1.0 / R::rgamma(shape, 1.0 / rate);
}

inline double log_sum_exp_cpp(const arma::vec& x) {
  double max_x = x.max();
  if (!std::isfinite(max_x)) return max_x;
  return max_x + std::log(arma::sum(arma::exp(x - max_x)));
}

inline arma::mat cov_to_cor(const arma::mat& S) {
  arma::vec sd = arma::sqrt(arma::clamp(S.diag(), 1e-12, arma::datum::inf));
  return S / (sd * sd.t());
}

// Log-likelihood contribution of nu for the observed-data multivariate-t
// densities, given each row's squared Mahalanobis distance to its allocated
// profile and number of observed variables (terms not involving nu, such as
// log|Sigma|, cancel in Metropolis ratios and are omitted).
inline double log_marg_nu(double nu, const arma::vec& d, const arma::vec& pobs) {
  double out = 0.0;
  for (arma::uword i = 0; i < d.n_elem; ++i) {
    const double q = pobs(i);
    if (q <= 0) continue;
    out += std::lgamma(0.5 * (nu + q)) - std::lgamma(0.5 * nu) - 0.5 * q * std::log(nu) -
      0.5 * (nu + q) * std::log1p(d(i) / nu);
  }
  return out;
}

// Univariate slice sampler with stepping out and shrinkage (Neal, 2003).
template <typename F>
inline double slice_sample(double x0, F logf, double w = 1.0, int m = 50) {
  double fx0 = logf(x0);
  double y = fx0 - R::rexp(1.0);
  double L = x0 - w * R::runif(0.0, 1.0);
  double Rr = L + w;
  int J = static_cast<int>(std::floor(m * R::runif(0.0, 1.0)));
  int K = m - 1 - J;
  while (J > 0 && logf(L) > y) { L -= w; J--; }
  while (K > 0 && logf(Rr) > y) { Rr += w; K--; }
  for (int it = 0; it < 200; ++it) {
    double x1 = L + R::runif(0.0, 1.0) * (Rr - L);
    if (logf(x1) > y) return x1;
    if (x1 < x0) L = x1; else Rr = x1;
  }
  return x0;
}

#endif
