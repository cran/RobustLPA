#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "rlpa_utils.h"

using namespace Rcpp;

// ===========================================================================
// Robust growth mixture model (latent class growth analysis) engine
//
// Long-format data: observation j of person i has outcome index outc(j)
// (0-based, one of K outcomes), time time(j) and value y(j); observations
// are sorted by person, and person i owns rows starts(i) .. starts(i+1)-1.
// For outcome k the fixed-effect design row is the polynomial basis
// (1, t, ..., t^degree) placed in block k of a K * (degree + 1) vector, and
// the random-effect design row is (1, t, ..., t^(nrand - 1)) placed in block
// k of a K * nrand vector. Within profile g:
//   y_i = X_i beta_g + Z_i b_i + e_i,  b_i ~ N(0, D_g / u_i),
//   e_i ~ N(0, R_ig / u_i),  R_ig = diag(sigma2_g[outcome]),
// with u_i = 1 (Gaussian) or u_i ~ Gamma(nu / 2, nu / 2) (multivariate t,
// Pinheiro, Liu & Wu, 2001), so that marginally
// y_i ~ N or t_nu(X_i beta_g, V_ig), V_ig = Z_i D_g Z_i' + R_ig.
// All person-level computations use the Woodbury identity with D = L L':
//   M = I + L' Z' R^{-1} Z L,  log|V| = log|R| + log|M|,
//   (y - Xb)' V^{-1} (y - Xb) = r'R^{-1}r - s'M^{-1}s,  s = L'Z'R^{-1} r,
//   E[b | y] = L M^{-1} s,  Var[b | y] (per unit scale) = L M^{-1} L',
// which only requires factorizing matrices of the random-effect dimension
// and remains valid when D is (near) singular.
// ===========================================================================

namespace {

struct PersonWork {
  arma::vec r;       // residual y - X beta
  arma::vec rinv;    // 1 / sigma2 of each row
  arma::mat Zt;      // Z L  (n_i x Q)
  arma::mat Lm;      // chol of M = I + Zt' R^{-1} Zt (lower)
  arma::vec s;       // Zt' R^{-1} r
  double logdetR;
};

inline double poly_term(double t, int power) {
  double out = 1.0;
  for (int l = 0; l < power; ++l) out *= t;
  return out;
}

// Fill the residuals and Woodbury pieces of person i for one profile.
inline void person_work(PersonWork& pw, int s0, int s1, const arma::vec& y, const arma::vec& time,
                        const arma::ivec& outc, const arma::vec& beta, const arma::mat& L,
                        const arma::vec& sig2, int degree, int nrand, int Q) {
  const int ni = s1 - s0;
  const int nb = degree + 1;
  pw.r.set_size(ni);
  pw.rinv.set_size(ni);
  pw.logdetR = 0.0;
  if (Q > 0) pw.Zt.zeros(ni, Q);
  for (int j = 0; j < ni; ++j) {
    const int row = s0 + j;
    const int k = outc(row);
    const double t = time(row);
    double xb = 0.0;
    for (int l = 0; l < nb; ++l) xb += beta(k * nb + l) * poly_term(t, l);
    pw.r(j) = y(row) - xb;
    pw.rinv(j) = 1.0 / sig2(k);
    pw.logdetR += std::log(sig2(k));
    if (Q > 0) {
      for (int l = 0; l < nrand; ++l) pw.Zt.row(j) += poly_term(t, l) * L.row(k * nrand + l);
    }
  }
  if (Q > 0) {
    arma::mat ZtR = pw.Zt.each_col() % pw.rinv;
    arma::mat M = pw.Zt.t() * ZtR;
    M.diag() += 1.0;
    pw.Lm = chol_lower_safe(M);
    pw.s = ZtR.t() * pw.r;
  }
}

// Squared Mahalanobis distance and log|V| from the Woodbury pieces.
inline void person_maha(const PersonWork& pw, int Q, double& maha, double& logdet) {
  maha = arma::dot(pw.r % pw.rinv, pw.r);
  logdet = pw.logdetR;
  if (Q > 0) {
    arma::vec v = arma::solve(arma::trimatl(pw.Lm), pw.s);
    maha -= arma::dot(v, v);
    logdet += 2.0 * arma::sum(arma::log(pw.Lm.diag()));
  }
  if (maha < 0.0) maha = 0.0;
}

}  // namespace

// Person-Level Log-Density of a Growth-Mixture Profile
//
// For every person, the exact marginal log-density of all of his/her
// observed values under one profile (Gaussian, dist = 0, or multivariate
// t, dist = 1), and the squared Mahalanobis distance of the observed
// vector to the profile's mean trajectory. Persons observed on different
// occasions and outcomes are handled exactly (full-information likelihood).
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List gmm_class_estep_cpp(const arma::vec& y, const arma::vec& time, const arma::ivec& outc,
                               const arma::ivec& starts, const arma::vec& beta, const arma::mat& D,
                               const arma::vec& sig2, int K, int degree, int nrand, int dist, double nu) {
  const int N = starts.n_elem - 1;
  const int Q = K * nrand;
  arma::mat L;
  if (Q > 0) L = chol_lower_safe(D);
  arma::vec logdens(N, arma::fill::zeros), maha(N, arma::fill::zeros), nobs(N, arma::fill::zeros);
  arma::vec logdetv(N, arma::fill::zeros);
  PersonWork pw;
  for (int i = 0; i < N; ++i) {
    const int s0 = starts(i), s1 = starts(i + 1);
    const int ni = s1 - s0;
    if (ni == 0) continue;
    person_work(pw, s0, s1, y, time, outc, beta, L, sig2, degree, nrand, Q);
    double d, logdet;
    person_maha(pw, Q, d, logdet);
    maha(i) = d;
    nobs(i) = ni;
    logdetv(i) = logdet;
    logdens(i) = log_const(dist, ni, logdet, nu) + log_kernel(dist, ni, d, nu);
  }
  return Rcpp::List::create(Rcpp::Named("logdens") = logdens, Rcpp::Named("maha") = maha,
                            Rcpp::Named("nobs") = nobs, Rcpp::Named("logdet") = logdetv);
}

// Sufficient Statistics of the ECM Step for One Growth-Mixture Profile
//
// Given posterior probabilities z, robustness weights w (all 1, Huber
// weights, or the E-step expectations u of the t latent scales) and the
// current parameters, returns:
//   A, c   : weighted normal equations of the fixed effects,
//            A = sum_i z_i w_i X_i' R^{-1} X_i,
//            c = sum_i z_i w_i X_i' R^{-1} (y_i - Z_i bhat_i);
//   Dnum, Dden : sum_i [z_i w_i bhat_i bhat_i' + c_i C_i] and sum_i c_i;
//   bhat   : E[b_i | y_i] for every person (N x Q);
//   zcz    : z_ij' C_i z_ij for every observation (conditional variance of
//            the random part of each observation),
// with C_i = Var(b_i | y_i) per unit scale and c_i = z_i (t model: exact ECM)
// or z_i w_i (Gaussian / Huber-weighted).
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List gmm_class_suffstats_cpp(const arma::vec& y, const arma::vec& time, const arma::ivec& outc,
                                   const arma::ivec& starts, const arma::vec& beta, const arma::mat& D,
                                   const arma::vec& sig2, int K, int degree, int nrand,
                                   const arma::vec& z, const arma::vec& w, bool t_scatter) {
  const int N = starts.n_elem - 1;
  const int Q = K * nrand;
  const int nb = degree + 1;
  const int P = K * nb;
  arma::mat L;
  if (Q > 0) L = chol_lower_safe(D);

  arma::mat A(P, P, arma::fill::zeros);
  arma::vec cvec(P, arma::fill::zeros);
  arma::mat Dnum(std::max(Q, 1), std::max(Q, 1), arma::fill::zeros);
  double Dden = 0.0;
  arma::mat bhat(N, std::max(Q, 1), arma::fill::zeros);
  arma::vec zcz(y.n_elem, arma::fill::zeros);
  PersonWork pw;

  for (int i = 0; i < N; ++i) {
    const int s0 = starts(i), s1 = starts(i + 1);
    const int ni = s1 - s0;
    if (ni == 0) continue;
    const double a = z(i) * w(i);
    const double cw = t_scatter ? z(i) : a;
    person_work(pw, s0, s1, y, time, outc, beta, L, sig2, degree, nrand, Q);

    arma::vec bh;
    arma::mat B;
    if (Q > 0) {
      arma::vec v = arma::solve(arma::trimatl(pw.Lm), pw.s);
      arma::vec Minv_s = arma::solve(arma::trimatu(pw.Lm.t()), v);
      bh = L * Minv_s;
      arma::mat Lm_inv = arma::inv(arma::trimatl(pw.Lm));
      B = L * Lm_inv.t();                 // C_i = B B'
      bhat.row(i) = bh.t();
      if (cw > 0.0 || a > 0.0) {
        Dnum += a * (bh * bh.t()) + cw * (B * B.t());
        Dden += cw;
      }
    }

    for (int j = 0; j < ni; ++j) {
      const int row = s0 + j;
      const int k = outc(row);
      const double t = time(row);
      double zb = 0.0;
      if (Q > 0) {
        arma::rowvec zrow(Q, arma::fill::zeros);
        for (int l = 0; l < nrand; ++l) zrow(k * nrand + l) = poly_term(t, l);
        zb = arma::dot(zrow, bh);
        arma::rowvec zB = zrow * B;
        zcz(row) = arma::dot(zB, zB);
      }
      if (a <= 0.0) continue;
      const double ri = pw.rinv(j);
      const double ystar = y(row) - zb;
      for (int l1 = 0; l1 < nb; ++l1) {
        const double x1 = poly_term(t, l1);
        cvec(k * nb + l1) += a * ri * x1 * ystar;
        for (int l2 = 0; l2 < nb; ++l2) A(k * nb + l1, k * nb + l2) += a * ri * x1 * poly_term(t, l2);
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("A") = A, Rcpp::Named("c") = cvec,
                            Rcpp::Named("Dnum") = Dnum, Rcpp::Named("Dden") = Dden,
                            Rcpp::Named("bhat") = bhat, Rcpp::Named("zcz") = zcz);
}

// Marginal (GLS) Normal Equations of the Fixed Effects for One Profile
//
// A = sum_i z_i w_i X_i' V_i^{-1} X_i and c = sum_i z_i w_i X_i' V_i^{-1} y_i,
// with V_i = Z_i D Z_i' + R_i evaluated at the current covariance
// parameters (Woodbury form). Maximizing the expected complete-data
// log-likelihood of the (y, u, class) formulation over the fixed effects
// for fixed covariance parameters gives the weighted GLS solution
// A^{-1} c; this is the first cycle of the alternating ECM algorithm used
// by robust_gmm(). Its curvature is that of the observed-data likelihood,
// which makes LASSO penalties on the fixed effects act at the right scale.
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List gmm_class_gls_cpp(const arma::vec& y, const arma::vec& time, const arma::ivec& outc,
                             const arma::ivec& starts, const arma::mat& D, const arma::vec& sig2,
                             int K, int degree, int nrand, const arma::vec& z, const arma::vec& w) {
  const int N = starts.n_elem - 1;
  const int Q = K * nrand;
  const int nb = degree + 1;
  const int P = K * nb;
  arma::mat L;
  if (Q > 0) L = chol_lower_safe(D);
  arma::mat A(P, P, arma::fill::zeros);
  arma::vec cvec(P, arma::fill::zeros);
  arma::vec zero_beta(P, arma::fill::zeros);
  PersonWork pw;
  for (int i = 0; i < N; ++i) {
    const int s0 = starts(i), s1 = starts(i + 1);
    const int ni = s1 - s0;
    const double a = z(i) * w(i);
    if (ni == 0 || a <= 0.0) continue;
    // with beta = 0, pw.r = y and pw.s = Zt' R^{-1} y
    person_work(pw, s0, s1, y, time, outc, zero_beta, L, sig2, degree, nrand, Q);
    arma::mat X(ni, P, arma::fill::zeros);
    for (int j = 0; j < ni; ++j) {
      const int row = s0 + j;
      const int k = outc(row);
      for (int l = 0; l < nb; ++l) X(j, k * nb + l) = poly_term(time(row), l);
    }
    arma::mat XR = (X.each_col() % pw.rinv).t();      // P x ni
    arma::mat XRX = XR * X;
    arma::vec XRy = XR * pw.r;
    if (Q > 0) {
      arma::mat XRZ = XR * pw.Zt;                     // P x Q
      arma::mat T1 = arma::solve(arma::trimatl(pw.Lm), XRZ.t());   // Lm^{-1} Z'R^{-1}X
      arma::vec t2 = arma::solve(arma::trimatl(pw.Lm), pw.s);      // Lm^{-1} Z'R^{-1}y
      XRX -= T1.t() * T1;
      XRy -= T1.t() * t2;
    }
    A += a * XRX;
    cvec += a * XRy;
  }
  return Rcpp::List::create(Rcpp::Named("A") = A, Rcpp::Named("c") = cvec);
}

// Residual-Variance Statistics for One Growth-Mixture Profile
//
// With the updated fixed effects beta and the conditional moments of the
// random effects from gmm_class_suffstats_cpp(), returns per outcome
//   num_k = sum_i sum_{j in k} [z_i w_i e_ij^2 + c_i zcz_ij],  den_k = sum c_i,
// where e_ij = y_ij - x_ij' beta - z_ij' bhat_i.
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List gmm_class_resid_cpp(const arma::vec& y, const arma::vec& time, const arma::ivec& outc,
                               const arma::ivec& starts, const arma::vec& beta, const arma::mat& bhat,
                               const arma::vec& zcz, int K, int degree, int nrand,
                               const arma::vec& z, const arma::vec& w, bool t_scatter) {
  const int N = starts.n_elem - 1;
  const int Q = K * nrand;
  const int nb = degree + 1;
  arma::vec num(K, arma::fill::zeros), den(K, arma::fill::zeros);
  for (int i = 0; i < N; ++i) {
    const double a = z(i) * w(i);
    const double cw = t_scatter ? z(i) : a;
    if (cw <= 0.0 && a <= 0.0) continue;
    for (int row = starts(i); row < starts(i + 1); ++row) {
      const int k = outc(row);
      const double t = time(row);
      double fit = 0.0;
      for (int l = 0; l < nb; ++l) fit += beta(k * nb + l) * poly_term(t, l);
      if (Q > 0) for (int l = 0; l < nrand; ++l) fit += bhat(i, k * nrand + l) * poly_term(t, l);
      const double e = y(row) - fit;
      num(k) += a * e * e + cw * zcz(row);
      den(k) += cw;
    }
  }
  return Rcpp::List::create(Rcpp::Named("num") = num, Rcpp::Named("den") = den);
}

// ===========================================================================
// MCMC
// ===========================================================================

namespace {

}  // namespace

// Gibbs Sampler for the Robust Growth Mixture Model
//
// One chain for the Gaussian (robust_type 0), Huber-weighted heuristic
// (robust_type 1) or multivariate-t (robust_type 2) growth mixture model on
// standardized outcomes. Each sweep:
//   1. allocates every person using the exact marginal (random effects and,
//      for t, latent scales integrated out) density of his/her data;
//   2. (t, estimate_nu) updates nu by random-walk Metropolis on log(nu)
//      with the latent scales integrated out (Gamma(2, 0.1) prior);
//   3. draws the latent scales u_i ~ Gamma((nu + n_i)/2, (nu + d_i)/2)
//      (random effects integrated out);
//   4. draws the fixed effects of each profile from their Gaussian
//      conditional with the random effects integrated out (then, 5, the
//      random effects b_i | beta, y_i, u_i: a blocked update of (beta, b)),
//      under a vague N(0, beta_prior_sd^2) prior, optionally
//      combined with Bayesian-Lasso (Laplace) priors (Park & Casella, 2008)
//      on the growth terms (rates rate_growth, G x P) and on the deviations
//      of each profile's coefficients from their across-profile centre
//      (rates rate_diff, G x P, or, with group_diff, one rate per outcome
//      rate_group, as in the Bayesian group Lasso of Kyung et al., 2010);
//      rates may differ by coefficient (adaptive Lasso);
//   5. draws the random effects given the fixed effects, and 6. the
//      random-effect covariance matrices (full, block-diagonal by outcome,
//      or diagonal; profile-specific or shared) by parameter expansion,
//      D = diag(alpha) Psi diag(alpha), under the scaled inverse-Wishart
//      prior (Psi ~ IW(q + 1, I), alpha ~ N(0, re_prior_scale^2); O'Malley
//      & Zaslavsky, 2008; Gelman et al., 2008), which removes the slow
//      mixing of small variance components;
//   7. draws the residual variances (inverse-gamma(1, 0.1) priors) and the
//      mixing proportions (Dirichlet(1, ..., 1)).
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List gmm_mcmc_chain_cpp(const arma::vec& y, const arma::vec& time, const arma::ivec& outc,
                              const arma::ivec& starts, int K, int degree, int nrand, int G, int mcmc_iter,
                              int robust_type, double alpha, double nu_init, bool estimate_nu,
                              bool re_equal, bool resid_equal, int re_structure,
                              const arma::mat& rate_growth, const arma::mat& rate_diff,
                              const arma::vec& rate_group, const arma::mat& group_scale, bool group_diff,
                              const arma::vec& beta_prior_sd, const arma::vec& re_prior_scale,
                              const arma::mat& beta0, const Rcpp::List& D0, const arma::mat& sig2_0,
                              const arma::vec& pi0) {
  const int N = starts.n_elem - 1;
  const int nb = degree + 1;
  const int P = K * nb;
  const int Q = K * nrand;
  const int dist = (robust_type == 2) ? 1 : 0;
  const int burnin = mcmc_iter / 2;

  arma::mat beta = beta0;                  // G x P
  std::vector<arma::mat> D(G);
  for (int g = 0; g < G; ++g) D[g] = Rcpp::as<arma::mat>(D0[g]);
  arma::mat sig2 = sig2_0;                 // G x K
  arma::vec pi_g = pi0 / arma::sum(pi0);
  double nu = nu_init;

  // random-effect covariance blocks
  std::vector<arma::uvec> blocks;
  if (Q > 0) {
    if (re_structure == 0) {
      blocks.push_back(arma::regspace<arma::uvec>(0, Q - 1));
    } else if (re_structure == 1) {
      for (int k = 0; k < K; ++k) blocks.push_back(arma::regspace<arma::uvec>(k * nrand, k * nrand + nrand - 1));
    } else {
      for (int q = 0; q < Q; ++q) blocks.push_back(arma::uvec({static_cast<arma::uword>(q)}));
    }
  }
  // Parameter-expanded random effects (Liu & Wu, 1999; Gelman et al., 2008):
  // b_i = alpha_g o xi_i, xi_i ~ N(0, Psi_g / u_i), so D_g = diag(alpha_g)
  // Psi_g diag(alpha_g) has a scaled inverse-Wishart prior (O'Malley &
  // Zaslavsky, 2008): Psi ~ IW(q + 1, I) per block, alpha_q ~ N(0, A_q^2).
  arma::mat px_scale(G, std::max(Q, 1), arma::fill::ones);
  std::vector<arma::mat> Psi(G);
  arma::mat xi(N, std::max(Q, 1), arma::fill::zeros);
  if (Q > 0) {
    for (int g = 0; g < G; ++g) {
      arma::vec sd0 = arma::sqrt(arma::clamp(D[g].diag(), 1e-8, arma::datum::inf));
      px_scale.row(g) = sd0.t();
      Psi[g] = D[g] / (sd0 * sd0.t());
    }
  }

  // Lasso state
  arma::uvec growth_mask(P, arma::fill::zeros);
  for (int k = 0; k < K; ++k) for (int l = 1; l < nb; ++l) growth_mask(k * nb + l) = 1;
  arma::mat tau_g(G, P, arma::fill::ones), tau_d(G, P, arma::fill::ones);
  arma::vec tau_group(K, arma::fill::ones);
  arma::vec centre = arma::mean(beta, 0).t();
  const bool use_diff = group_diff ? arma::any(rate_group > 0.0) : arma::any(arma::vectorise(rate_diff) > 0.0);
  for (int g = 0; g < G; ++g)
    for (int p = 0; p < P; ++p) {
      if (!(growth_mask(p) && rate_growth(g, p) > 0.0)) tau_g(g, p) = 1e10;
      if (!group_diff && rate_diff(g, p) <= 0.0) tau_d(g, p) = 1e10;
    }
  for (int k = 0; k < K; ++k) if (group_diff && rate_group(k) <= 0.0) tau_group(k) = 1e10;

  arma::vec nobs(N), dstar(N, arma::fill::zeros), uval(N, arma::fill::ones), wgt(N, arma::fill::ones);
  arma::uvec alloc(N, arma::fill::zeros);
  for (int i = 0; i < N; ++i) nobs(i) = starts(i + 1) - starts(i);
  arma::mat b(N, std::max(Q, 1), arma::fill::zeros);

  arma::mat beta_chain(mcmc_iter, G * P), sig2_chain(mcmc_iter, G * K), pi_chain(mcmc_iter, G);
  arma::mat D_chain(mcmc_iter, std::max(G * Q * Q, 1), arma::fill::zeros);
  arma::vec nu_chain(mcmc_iter);
  double step_nu = 0.3;
  int acc_nu = 0, try_nu = 0;
  PersonWork pw;

  for (int iter = 0; iter < mcmc_iter; ++iter) {
    Rcpp::checkUserInterrupt();
    std::vector<arma::mat> Lg(G);
    if (Q > 0) for (int g = 0; g < G; ++g) Lg[g] = chol_lower_safe(D[g]);

    // ---- 1. allocation -------------------------------------------------------
    for (int i = 0; i < N; ++i) {
      const int s0 = starts(i), s1 = starts(i + 1);
      const int ni = s1 - s0;
      if (ni == 0) {
        alloc(i) = sample_class(pi_g);
        dstar(i) = 0.0;
        continue;
      }
      arma::vec logp(G), dvec(G);
      for (int g = 0; g < G; ++g) {
        arma::vec bg = beta.row(g).t();
        arma::vec sg = sig2.row(g).t();
        person_work(pw, s0, s1, y, time, outc, bg, Lg[g], sg, degree, nrand, Q);
        double d, logdet;
        person_maha(pw, Q, d, logdet);
        dvec(g) = d;
        logp(g) = std::log(std::max(pi_g(g), 1e-300)) + log_const(dist, ni, logdet, nu) + log_kernel(dist, ni, d, nu);
      }
      arma::vec probs = arma::exp(logp - log_sum_exp_cpp(logp));
      const int gs = sample_class(probs);
      alloc(i) = gs;
      dstar(i) = dvec(gs);
    }

    // ---- 2. degrees of freedom (latent scales integrated out) ----------------
    if (robust_type == 2 && estimate_nu) {
      double prop = nu * std::exp(step_nu * R::rnorm(0.0, 1.0));
      try_nu++;
      if (prop >= 1.0 && prop <= 1000.0) {
        double log_r = log_marg_nu(prop, dstar, nobs) - log_marg_nu(nu, dstar, nobs) +
          (std::log(prop) - 0.1 * prop) - (std::log(nu) - 0.1 * nu) + std::log(prop) - std::log(nu);
        if (std::log(R::runif(0.0, 1.0)) < log_r) {
          nu = prop;
          acc_nu++;
        }
      }
    }

    // ---- 3. latent scales / robustness weights (random effects integrated out) --
    for (int i = 0; i < N; ++i) {
      const int ni = starts(i + 1) - starts(i);
      double u = 1.0, w = 1.0;
      if (robust_type == 2) {
        u = R::rgamma(0.5 * (nu + ni), 1.0 / (0.5 * (nu + dstar(i))));
        w = u;
      } else if (robust_type == 1 && ni > 0) {
        const double cut = R::qchisq(1.0 - alpha, ni, 1, 0);
        if (dstar(i) > cut) w = std::sqrt(cut / dstar(i));
      }
      uval(i) = u;
      wgt(i) = w;
    }

    // ---- 4. fixed effects, random effects integrated out (blocked update) ------
    // beta_g | u, allocation, D, R ~ N(A^{-1} c, A^{-1}) with
    // A = sum_i w_i X_i' V_i^{-1} X_i + prior precision; drawing beta
    // marginally over b_i and then b_i | beta (step 5) avoids the slow mixing
    // of the conditional (beta | b) Gibbs update for the intercepts.
    arma::vec n_counts(G, arma::fill::zeros);
    for (int i = 0; i < N; ++i) n_counts(alloc(i)) += 1.0;
    arma::vec zero_beta(P, arma::fill::zeros);
    for (int g = 0; g < G; ++g) {
      if (n_counts(g) == 0) continue;
      arma::mat Pm(P, P, arma::fill::zeros);
      arma::vec rhs(P, arma::fill::zeros);
      arma::vec sg = sig2.row(g).t();
      for (int i = 0; i < N; ++i) {
        if (static_cast<int>(alloc(i)) != g) continue;
        const int s0 = starts(i), s1 = starts(i + 1);
        const int ni = s1 - s0;
        if (ni == 0) continue;
        const double wi = wgt(i);
        person_work(pw, s0, s1, y, time, outc, zero_beta, Lg[g], sg, degree, nrand, Q);
        arma::mat X(ni, P, arma::fill::zeros);
        for (int j = 0; j < ni; ++j) {
          const int row = s0 + j;
          const int k = outc(row);
          for (int l = 0; l < nb; ++l) X(j, k * nb + l) = poly_term(time(row), l);
        }
        arma::mat XR = (X.each_col() % pw.rinv).t();
        arma::mat XVX = XR * X;
        arma::vec XVy = XR * pw.r;
        if (Q > 0) {
          arma::mat T1 = arma::solve(arma::trimatl(pw.Lm), (XR * pw.Zt).t());
          arma::vec t2 = arma::solve(arma::trimatl(pw.Lm), pw.s);
          XVX -= T1.t() * T1;
          XVy -= T1.t() * t2;
        }
        Pm += wi * XVX;
        rhs += wi * XVy;
      }
      for (int p = 0; p < P; ++p) {
        double prec = 1.0 / (beta_prior_sd(p) * beta_prior_sd(p));
        double pm = 0.0;
        if (rate_growth(g, p) > 0.0 && growth_mask(p)) prec += 1.0 / tau_g(g, p);
        if (use_diff) {
          const double td = group_diff ? tau_group(p / nb) / (group_scale(g, p) * group_scale(g, p)) : tau_d(g, p);
          prec += 1.0 / td;
          pm += centre(p) / td;
        }
        Pm(p, p) += prec;
        rhs(p) += pm;
      }
      Pm = 0.5 * (Pm + Pm.t());
      arma::mat Lp = chol_lower_safe(Pm);
      arma::vec mean = arma::solve(arma::trimatu(Lp.t()), arma::solve(arma::trimatl(Lp), rhs));
      arma::vec eps(P, arma::fill::randn);
      beta.row(g) = (mean + arma::solve(arma::trimatu(Lp.t()), eps)).t();
    }

    // ---- 5. random effects b_i | beta, u_i ------------------------------------------
    if (Q > 0) {
      for (int i = 0; i < N; ++i) {
        const int s0 = starts(i), s1 = starts(i + 1);
        const int ni = s1 - s0;
        const int gs = alloc(i);
        const double u = uval(i);
        arma::vec eps(Q, arma::fill::randn);
        if (ni == 0) {
          b.row(i) = (Lg[gs] * eps / std::sqrt(u)).t();
          continue;
        }
        arma::vec bg = beta.row(gs).t();
        arma::vec sg = sig2.row(gs).t();
        person_work(pw, s0, s1, y, time, outc, bg, Lg[gs], sg, degree, nrand, Q);
        arma::vec v = arma::solve(arma::trimatl(pw.Lm), pw.s);
        arma::vec Minv_s = arma::solve(arma::trimatu(pw.Lm.t()), v);
        arma::vec bh = Lg[gs] * Minv_s;
        arma::vec dev = arma::solve(arma::trimatu(pw.Lm.t()), eps);
        b.row(i) = (bh + Lg[gs] * dev / std::sqrt(u)).t();
      }
    }

    // ---- Lasso hyper-parameters --------------------------------------------------
    if (use_diff) {
      for (int p = 0; p < P; ++p) {
        double prec = 0.0, s = 0.0;
        for (int g = 0; g < G; ++g) {
          const double td = group_diff ? tau_group(p / nb) / (group_scale(g, p) * group_scale(g, p)) : tau_d(g, p);
          prec += 1.0 / td;
          s += beta(g, p) / td;
        }
        centre(p) = s / prec + R::rnorm(0.0, 1.0) / std::sqrt(prec);
      }
      if (group_diff) {
        for (int k = 0; k < K; ++k) {
          const double lam = rate_group(k);
          if (lam <= 0.0) { tau_group(k) = 1e10; continue; }
          double ss = 0.0;
          for (int g = 0; g < G; ++g)
            for (int l = 0; l < nb; ++l) {
              const double dv = group_scale(g, k * nb + l) * (beta(g, k * nb + l) - centre(k * nb + l));
              ss += dv * dv;
            }
          tau_group(k) = 1.0 / rinvgauss(lam / std::max(std::sqrt(ss), 1e-8), lam * lam);
        }
      } else {
        for (int g = 0; g < G; ++g)
          for (int p = 0; p < P; ++p)
            if (rate_diff(g, p) <= 0.0) tau_d(g, p) = 1e10; else tau_d(g, p) = 1.0 / rinvgauss(rate_diff(g, p) / std::max(std::abs(beta(g, p) - centre(p)), 1e-8),
                                          rate_diff(g, p) * rate_diff(g, p));
      }
    }
    for (int g = 0; g < G; ++g)
      for (int p = 0; p < P; ++p)
        if (growth_mask(p) && rate_growth(g, p) > 0.0)
          tau_g(g, p) = 1.0 / rinvgauss(rate_growth(g, p) / std::max(std::abs(beta(g, p)), 1e-8),
                                        rate_growth(g, p) * rate_growth(g, p));

    // ---- 6. random-effect covariances by parameter expansion ----------------------
    if (Q > 0) {
      // working random effects
      for (int i = 0; i < N; ++i) xi.row(i) = b.row(i) / px_scale.row(alloc(i));
      // (a) Psi | xi: inverse-Wishart per block, profile-specific or pooled
      std::vector<arma::mat> Sg(G, arma::mat(Q, Q, arma::fill::zeros));
      arma::vec ng(G, arma::fill::zeros);
      for (int i = 0; i < N; ++i) {
        const int g = alloc(i);
        arma::rowvec xr = xi.row(i);
        Sg[g] += wgt(i) * (xr.t() * xr);
        ng(g) += (robust_type == 1) ? wgt(i) : 1.0;
      }
      const int n_groups = re_equal ? 1 : G;
      for (int gg = 0; gg < n_groups; ++gg) {
        arma::mat S(Q, Q, arma::fill::zeros);
        double n = 0.0;
        if (re_equal) { for (int g = 0; g < G; ++g) { S += Sg[g]; n += ng(g); } }
        else { S = Sg[gg]; n = ng(gg); }
        arma::mat Pnew(Q, Q, arma::fill::zeros);
        for (size_t bl = 0; bl < blocks.size(); ++bl) {
          const arma::uvec& idx = blocks[bl];
          const double qb = idx.n_elem;
          arma::mat Sb = S.submat(idx, idx) + arma::eye<arma::mat>(idx.n_elem, idx.n_elem);
          Pnew.submat(idx, idx) = rinvwishart(qb + 1.0 + n, 0.5 * (Sb + Sb.t()));
        }
        Psi[gg] = Pnew;
      }
      // (b) alpha | xi, beta, u, sigma2: Gaussian regression of the residuals
      //     y - X beta on the columns z_q * xi_q
      for (int gg = 0; gg < n_groups; ++gg) {
        arma::mat Pa = arma::diagmat(1.0 / arma::square(re_prior_scale));
        arma::vec ra(Q, arma::fill::zeros);
        for (int i = 0; i < N; ++i) {
          const int g = alloc(i);
          if (!re_equal && g != gg) continue;
          const double wi = wgt(i);
          for (int row = starts(i); row < starts(i + 1); ++row) {
            const int k = outc(row);
            const double t = time(row);
            double fit = 0.0;
            for (int l = 0; l < nb; ++l) fit += beta(g, k * nb + l) * poly_term(t, l);
            const double r = y(row) - fit;
            const double prec = wi / sig2(g, k);
            arma::vec wv(Q, arma::fill::zeros);
            for (int l = 0; l < nrand; ++l) wv(k * nrand + l) = poly_term(t, l) * xi(i, k * nrand + l);
            Pa += prec * (wv * wv.t());
            ra += prec * wv * r;
          }
        }
        Pa = 0.5 * (Pa + Pa.t());
        arma::mat La = chol_lower_safe(Pa);
        arma::vec ma = arma::solve(arma::trimatu(La.t()), arma::solve(arma::trimatl(La), ra));
        arma::vec eps(Q, arma::fill::randn);
        arma::vec anew = ma + arma::solve(arma::trimatu(La.t()), eps);
        for (int q = 0; q < Q; ++q) if (std::abs(anew(q)) < 1e-8) anew(q) = (anew(q) < 0 ? -1e-8 : 1e-8);
        px_scale.row(gg) = anew.t();
      }
      if (re_equal) for (int g = 1; g < G; ++g) { px_scale.row(g) = px_scale.row(0); Psi[g] = Psi[0]; }
      for (int g = 0; g < G; ++g) {
        arma::vec a = px_scale.row(g).t();
        D[g] = Psi[g] % (a * a.t());
      }
      for (int i = 0; i < N; ++i) b.row(i) = xi.row(i) % px_scale.row(alloc(i));
    }

    // ---- 7. residual variances and mixing proportions ------------------------------
    {
      arma::mat sse(G, K, arma::fill::zeros), cnt(G, K, arma::fill::zeros);
      for (int i = 0; i < N; ++i) {
        const int g = alloc(i);
        const double wi = wgt(i);
        for (int row = starts(i); row < starts(i + 1); ++row) {
          const int k = outc(row);
          const double t = time(row);
          double fit = 0.0;
          for (int l = 0; l < nb; ++l) fit += beta(g, k * nb + l) * poly_term(t, l);
          if (Q > 0) for (int l = 0; l < nrand; ++l) fit += b(i, k * nrand + l) * poly_term(t, l);
          const double e = y(row) - fit;
          sse(g, k) += wi * e * e;
          cnt(g, k) += (robust_type == 1) ? wi : 1.0;
        }
      }
      if (resid_equal) {
        arma::rowvec sp = arma::sum(sse, 0), cp = arma::sum(cnt, 0);
        for (int k = 0; k < K; ++k) {
          const double v = rinvgamma_cpp(1.0 + 0.5 * cp(k), 0.1 + 0.5 * sp(k));
          for (int g = 0; g < G; ++g) sig2(g, k) = v;
        }
      } else {
        for (int g = 0; g < G; ++g)
          for (int k = 0; k < K; ++k) sig2(g, k) = rinvgamma_cpp(1.0 + 0.5 * cnt(g, k), 0.1 + 0.5 * sse(g, k));
      }
      arma::vec gd(G);
      for (int g = 0; g < G; ++g) gd(g) = R::rgamma(n_counts(g) + 1.0, 1.0);
      pi_g = gd / arma::sum(gd);
    }

    // ---- adapt the nu step during burn-in -----------------------------------------
    if (iter < burnin && (iter + 1) % 50 == 0 && try_nu > 0) {
      const double rate = static_cast<double>(acc_nu) / try_nu;
      if (rate > 0.45) step_nu *= 1.25; else if (rate < 0.3) step_nu /= 1.25;
      step_nu = std::min(std::max(step_nu, 0.02), 2.0);
      acc_nu = try_nu = 0;
    }

    // ---- store -------------------------------------------------------------------
    for (int g = 0; g < G; ++g) {
      for (int p = 0; p < P; ++p) beta_chain(iter, g * P + p) = beta(g, p);
      for (int k = 0; k < K; ++k) sig2_chain(iter, g * K + k) = sig2(g, k);
      if (Q > 0) for (int q = 0; q < Q * Q; ++q) D_chain(iter, g * Q * Q + q) = D[g](q);
      pi_chain(iter, g) = pi_g(g);
    }
    nu_chain(iter) = nu;
  }

  return Rcpp::List::create(
    Rcpp::Named("beta_chain") = beta_chain, Rcpp::Named("D_chain") = D_chain,
    Rcpp::Named("sig2_chain") = sig2_chain, Rcpp::Named("pi_chain") = pi_chain,
    Rcpp::Named("nu_chain") = nu_chain
  );
}
