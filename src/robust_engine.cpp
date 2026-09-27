#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "rlpa_utils.h"

using namespace Rcpp;

// ===========================================================================
// Internal helpers (not exported to R)
// ===========================================================================






// Convert an R list of (0-based) integer vectors into arma::uvec's.
static std::vector<arma::uvec> as_uvec_list(const Rcpp::List& L) {
  std::vector<arma::uvec> out(L.size());
  for (int k = 0; k < L.size(); ++k) {
    Rcpp::IntegerVector v = L[k];
    arma::uvec u(v.size());
    for (int j = 0; j < v.size(); ++j) u(j) = static_cast<arma::uword>(v[j]);
    out[k] = u;
  }
  return out;
}

// Indices in 0..p-1 that are NOT in `obs` (the missing variables of a pattern).
static arma::uvec complement_idx(const arma::uvec& obs, arma::uword p) {
  std::vector<bool> is_obs(p, false);
  for (arma::uword j = 0; j < obs.n_elem; ++j) is_obs[obs(j)] = true;
  std::vector<arma::uword> mis;
  for (arma::uword j = 0; j < p; ++j) if (!is_obs[j]) mis.push_back(j);
  return arma::uvec(mis);
}





// ===========================================================================
// Simple robust location (used by robust_mean())
// ===========================================================================

// Trimmed Centroid Based on Distance to the Coordinate-Wise Median
//
// An observation is included in the average only if its Euclidean distance
// to the coordinate-wise median of the data is at most `threshold`.
//
// @param X A numeric matrix (no missing values).
// @param threshold Maximum Euclidean distance to the coordinate-wise median.
// @return A numeric row vector; zeros (with a warning) if no observation
//   falls within `threshold` of the median.
// @keywords internal
// @noRd
// [[Rcpp::export]]
arma::rowvec robust_mean_cpp(arma::mat X, double threshold) {
  int n = X.n_rows;
  int p = X.n_cols;
  arma::rowvec center(p);
  for (int j = 0; j < p; j++) center(j) = arma::median(X.col(j));

  arma::rowvec sum = arma::zeros<arma::rowvec>(p);
  int valid_count = 0;
  for (int i = 0; i < n; i++) {
    arma::rowvec current_row = X.row(i);
    if (arma::norm(current_row - center) <= threshold) {
      sum += current_row;
      valid_count++;
    }
  }
  if (valid_count > 0) return sum / valid_count;
  Rcpp::warning("robust_mean_cpp: no observation within `threshold` of the median; returning a zero vector.");
  return arma::zeros<arma::rowvec>(p);
}

// ===========================================================================
// EM engine: E-step and M-step for a single profile
// ===========================================================================

// Per-Observation Log-Density and Mahalanobis Distance for One Profile
//
// Evaluates, for every row of X, the log-density of its *observed* entries
// under a Gaussian (dist = 0) or multivariate t (dist = 1) profile with
// location `mu` and scale matrix `Sigma` -- i.e. the exact observed-data
// (FIML) likelihood contribution -- together with the squared Mahalanobis
// distance on the observed entries. Rows are processed in blocks sharing
// the same missingness pattern, so each pattern's sub-covariance is
// factorized only once. Everything is computed on the log scale, so
// arbitrarily extreme observations never underflow.
//
// @param X Numeric matrix, NA allowed.
// @param mu Numeric row vector (length p).
// @param Sigma p x p covariance (Gaussian) or scale (t) matrix.
// @param pat_obs List of 0-based observed-variable indices, one per pattern.
// @param pat_rows List of 0-based row indices, one per pattern.
// @param dist 0 = Gaussian, 1 = multivariate t.
// @param nu Degrees of freedom (dist = 1 only).
// @return List with `logdens`, `maha`, `pobs` (numeric vectors of length n).
//   Rows with no observed variable get logdens = 0, maha = 0, pobs = 0.
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List class_estep_cpp(const arma::mat& X, const arma::rowvec& mu, const arma::mat& Sigma,
                           const Rcpp::List& pat_obs, const Rcpp::List& pat_rows,
                           int dist, double nu) {
  const arma::uword n = X.n_rows;
  std::vector<arma::uvec> obs = as_uvec_list(pat_obs);
  std::vector<arma::uvec> rows = as_uvec_list(pat_rows);

  arma::vec logdens(n, arma::fill::zeros), maha(n, arma::fill::zeros), pobs(n, arma::fill::zeros);
  arma::vec mu_c = mu.t();

  for (size_t k = 0; k < obs.size(); ++k) {
    const arma::uvec& o = obs[k];
    const arma::uvec& r = rows[k];
    const int po = o.n_elem;
    if (po == 0 || r.n_elem == 0) continue;

    arma::mat L = chol_lower_safe(Sigma.submat(o, o));
    double logdet = 2.0 * arma::sum(arma::log(L.diag()));
    double cst = log_const(dist, po, logdet, nu);

    arma::mat D = X.submat(r, o);
    D.each_row() -= mu_c.elem(o).t();
    arma::mat Y = arma::solve(arma::trimatl(L), D.t());   // po x n_k
    arma::rowvec d = arma::sum(arma::square(Y), 0);

    for (arma::uword ii = 0; ii < r.n_elem; ++ii) {
      arma::uword i = r(ii);
      maha(i) = d(ii);
      pobs(i) = po;
      logdens(i) = cst + log_kernel(dist, po, d(ii), nu);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("logdens") = logdens,
    Rcpp::Named("maha") = maha,
    Rcpp::Named("pobs") = pobs
  );
}

// Weighted M-Step for One Profile with Exact EM Treatment of Missing Values
//
// Given posterior membership probabilities `z` and robustness weights `w`
// (Huber weights, multivariate-t latent-scale expectations u, or all 1s), and
// the profile's *current* parameters (mu_cur, Sigma_cur), computes the
// updated mean and (unconstrained) covariance. Missing entries are handled
// by the standard EM for incomplete multivariate normal / t data
// (Ghahramani & Jordan, 1994; Liu & Rubin, 1995): each missing block is
// replaced by its conditional expectation given the observed block,
//   xhat_m = mu_m + S_mo S_oo^{-1} (x_o - mu_o),
// and the conditional covariance C_mm = S_mm - S_mo S_oo^{-1} S_om is added
// to the scatter matrix. This yields maximum-likelihood estimates under
// MAR, unlike pairwise available-case moments.
//
//   mean       = sum_i z_i w_i xhat_i / sum_i z_i w_i       (then LASSO soft-threshold)
//   covariance = sum_i [ z_i w_i (xhat_i - mean)(xhat_i - mean)' + c_i C_i ] / denom
//
// with c_i = z_i and denom = sum_i z_i when `t_scatter` is TRUE (the exact
// ECM update of a multivariate-t mixture), and c_i = z_i w_i and
// denom = sum_i z_i w_i otherwise (Gaussian when all w = 1; Huber-type
// weighted covariance when w are Huber weights).
//
// @param X Numeric matrix, NA allowed.
// @param z Posterior probabilities for this profile (length n).
// @param w Robustness weights (length n).
// @param mu_cur,Sigma_cur Current profile parameters (used for conditional
//   expectations of missing entries).
// @param pat_obs,pat_rows Missingness patterns (see class_estep_cpp()).
// @param lambda Non-negative soft-thresholding penalty for the mean.
// @param t_scatter Logical, see above.
// @return List with `mean` (row vector) and `covariance` (matrix).
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List class_mstep_cpp(const arma::mat& X, const arma::vec& z, const arma::vec& w,
                           const arma::rowvec& mu_cur, const arma::mat& Sigma_cur,
                           const Rcpp::List& pat_obs, const Rcpp::List& pat_rows,
                           double lambda, bool t_scatter) {
  const arma::uword p = X.n_cols;
  std::vector<arma::uvec> obs = as_uvec_list(pat_obs);
  std::vector<arma::uvec> rows = as_uvec_list(pat_rows);

  arma::vec a = z % w;
  arma::vec c = t_scatter ? z : a;
  double denom_mu = arma::sum(a);
  double denom_cov = t_scatter ? arma::sum(z) : denom_mu;

  if (denom_mu < 1e-10 || denom_cov < 1e-10) {
    // Profile has (numerically) no members: keep its current parameters.
    return Rcpp::List::create(Rcpp::Named("mean") = mu_cur, Rcpp::Named("covariance") = Sigma_cur);
  }

  arma::mat Xhat = X;
  arma::mat Csum(p, p, arma::fill::zeros);
  arma::vec mu_c = mu_cur.t();

  for (size_t k = 0; k < obs.size(); ++k) {
    const arma::uvec& o = obs[k];
    const arma::uvec& r = rows[k];
    if (r.n_elem == 0) continue;
    arma::uvec m = complement_idx(o, p);
    if (m.n_elem == 0) continue;

    double c_k = arma::sum(c.elem(r));

    if (o.n_elem == 0) {
      for (arma::uword ii = 0; ii < r.n_elem; ++ii) Xhat.row(r(ii)) = mu_cur;
      Csum += c_k * Sigma_cur;
      continue;
    }

    arma::mat L = chol_lower_safe(Sigma_cur.submat(o, o));
    arma::mat Smo = Sigma_cur.submat(m, o);
    // B = S_mo S_oo^{-1}, via two triangular solves on S_oo = L L'
    arma::mat tmp = arma::solve(arma::trimatl(L), Smo.t());
    arma::mat Bt = arma::solve(arma::trimatu(L.t()), tmp);    // (po x pm) = S_oo^{-1} S_om
    arma::mat Cmm = Sigma_cur.submat(m, m) - Smo * Bt;
    Cmm = 0.5 * (Cmm + Cmm.t());

    arma::mat D = X.submat(r, o);
    D.each_row() -= mu_c.elem(o).t();
    arma::mat Xm = D * Bt;                                    // n_k x pm
    Xm.each_row() += mu_c.elem(m).t();
    Xhat.submat(r, m) = Xm;

    Csum.submat(m, m) += c_k * Cmm;
  }

  arma::rowvec new_mu = (Xhat.t() * a).t() / denom_mu;

  if (lambda > 0.0) {
    for (arma::uword j = 0; j < p; ++j) {
      if (new_mu(j) > lambda) new_mu(j) -= lambda;
      else if (new_mu(j) < -lambda) new_mu(j) += lambda;
      else new_mu(j) = 0.0;
    }
  }

  arma::mat Xc = Xhat.each_row() - new_mu;
  arma::mat S = Xc.t() * (Xc.each_col() % a);
  arma::mat new_sigma = (S + Csum) / denom_cov;
  new_sigma = 0.5 * (new_sigma + new_sigma.t());

  return Rcpp::List::create(
    Rcpp::Named("mean") = new_mu,
    Rcpp::Named("covariance") = new_sigma
  );
}

// Pairwise Available-Case Weighted Moments (Initialization Only)
//
// Posterior-probability-weighted mean and covariance using, for each
// variable (pair), only the rows where it is observed. Not a maximum
// likelihood estimator under missingness: used only to obtain starting
// values before the first exact EM M-step (class_mstep_cpp()).
//
// @param X Numeric matrix, NA allowed.
// @param z Non-negative weights (length n).
// @return List with `mean` and `covariance`.
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List pairwise_moments_cpp(const arma::mat& X, const arma::vec& z) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::rowvec new_mu(p, arma::fill::zeros);
  arma::mat new_sigma(p, p, arma::fill::zeros);

  for (int j = 0; j < p; j++) {
    double sum_val = 0, sum_w = 0;
    for (int i = 0; i < n; i++) {
      if (std::isfinite(X(i, j))) {
        sum_val += z(i) * X(i, j);
        sum_w += z(i);
      }
    }
    new_mu(j) = (sum_w > 1e-10) ? (sum_val / sum_w) : 0.0;
  }
  for (int j = 0; j < p; j++) {
    for (int k = 0; k <= j; k++) {
      double sum_val = 0, sum_w = 0;
      for (int i = 0; i < n; i++) {
        if (std::isfinite(X(i, j)) && std::isfinite(X(i, k))) {
          sum_val += z(i) * (X(i, j) - new_mu(j)) * (X(i, k) - new_mu(k));
          sum_w += z(i);
        }
      }
      double v = (sum_w > 1e-10) ? (sum_val / sum_w) : ((j == k) ? 1.0 : 0.0);
      new_sigma(j, k) = v;
      new_sigma(k, j) = v;
    }
  }
  return Rcpp::List::create(Rcpp::Named("mean") = new_mu, Rcpp::Named("covariance") = new_sigma);
}

// Linear Sum Assignment (Hungarian Algorithm)
//
// Solves min_perm sum_i cost[i, perm[i]] exactly for a square cost matrix,
// in O(G^3). Used to align arbitrarily labeled profiles (MCMC draws,
// bootstrap refits) to a reference labeling.
//
// @param cost A square numeric matrix (non-finite entries are treated as a
//   very large cost).
// @return Integer vector (0-based): element i is the column assigned to row i.
// @keywords internal
// @noRd
// [[Rcpp::export]]
arma::uvec solve_lsap_cpp(arma::mat cost) {
  const int n = cost.n_rows;
  if (n == 0) return arma::uvec();
  double big = 1.0;
  for (arma::uword k = 0; k < cost.n_elem; ++k) if (std::isfinite(cost(k))) big = std::max(big, std::abs(cost(k)));
  for (arma::uword k = 0; k < cost.n_elem; ++k) if (!std::isfinite(cost(k))) cost(k) = 1e6 * big;

  const double INF = std::numeric_limits<double>::infinity();
  std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0);
  std::vector<int> pp(n + 1, 0), way(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    pp[0] = i;
    int j0 = 0;
    std::vector<double> minv(n + 1, INF);
    std::vector<char> used(n + 1, false);
    do {
      used[j0] = true;
      int i0 = pp[j0], j1 = 0;
      double delta = INF;
      for (int j = 1; j <= n; ++j) {
        if (!used[j]) {
          double cur = cost(i0 - 1, j - 1) - u[i0] - v[j];
          if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
          if (minv[j] < delta) { delta = minv[j]; j1 = j; }
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) { u[pp[j]] += delta; v[j] -= delta; }
        else minv[j] -= delta;
      }
      j0 = j1;
    } while (pp[j0] != 0);
    do {
      int j1 = way[j0];
      pp[j0] = pp[j1];
      j0 = j1;
    } while (j0);
  }
  arma::uvec ans(n);
  for (int j = 1; j <= n; ++j) ans(pp[j] - 1) = j - 1;
  return ans;
}

// ===========================================================================
// MCMC engine
// ===========================================================================

















// Draw the inverse scales a = 1/d of a covariance D R D given R, one
// coordinate at a time from its exact full conditional
//   p(a_j | .) propto a_j^k exp(-A a_j^2 - B a_j),
// with k = df + 2*alpha0 - 1, A = M_jj / 2 + beta0, B = sum_{l != j} M_jl a_l,
// M = R^{-1} o S (Hadamard product; summed over profiles when D is shared),
// under an InvGamma(alpha0 = 1, beta0 = 1) prior on each variance d_j^2 --
// the same prior used for the variances of the diagonal models.
static arma::vec draw_inverse_scales(const arma::mat& M, arma::vec a, double df) {
  const double k = df + 1.0;   // df + 2 * 1 - 1
  const arma::uword p = a.n_elem;
  for (arma::uword j = 0; j < p; ++j) {
    double B = arma::dot(M.row(j), a.t()) - M(j, j) * a(j);
    double A = 0.5 * M(j, j) + 1.0;
    auto logf = [k, A, B](double x) {
      double aa = std::exp(x);
      return (k + 1.0) * x - A * aa * aa - B * aa;   // includes the log-scale Jacobian
    };
    a(j) = std::exp(slice_sample(std::log(a(j)), logf));
  }
  return a;
}

// Log full conditional of a correlation matrix R (uniform LKJ(1) prior)
// given the standardized scatter T and the total degrees of freedom df:
//   -df/2 log|R| - tr(R^{-1} T) / 2.  Returns -Inf if R is not PD.
static double log_target_cor(const arma::mat& Rm, const arma::mat& T, double df) {
  arma::mat L;
  if (!arma::chol(L, Rm, "lower")) return -arma::datum::inf;
  double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::mat Linv = arma::inv(arma::trimatl(L));
  arma::mat Rinv = Linv.t() * Linv;
  return -0.5 * df * logdet - 0.5 * arma::accu(Rinv % T);
}

// One sweep of element-wise random-walk Metropolis updates of the
// off-diagonal entries of a correlation matrix (symmetric proposal, so the
// acceptance ratio is the target ratio; proposals leaving the set of
// positive-definite correlation matrices are rejected).
static arma::mat mh_correlation(arma::mat Rm, const arma::mat& T, double df, double step,
                                int& accepted, int& tried) {
  const arma::uword p = Rm.n_rows;
  double cur = log_target_cor(Rm, T, df);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword l = j + 1; l < p; ++l) {
      double prop = Rm(j, l) + step * R::rnorm(0.0, 1.0);
      tried++;
      if (std::abs(prop) >= 1.0) continue;
      arma::mat Rp = Rm;
      Rp(j, l) = prop;
      Rp(l, j) = prop;
      double lp = log_target_cor(Rp, T, df);
      if (std::isfinite(lp) && std::log(R::runif(0.0, 1.0)) < lp - cur) {
        Rm = Rp;
        cur = lp;
        accepted++;
      }
    }
  }
  return Rm;
}

struct PatCache {
  arma::mat L;      // chol of Sigma_oo
  double logdet;
  arma::mat Bt;     // Sigma_oo^{-1} Sigma_om  (po x pm)
  arma::mat Lc;     // chol of the conditional covariance of the missing block
};

// Gibbs Sampler for (Robust) Bayesian Latent Profile Analysis
//
// Runs one MCMC chain for a Gaussian (robust_type 0 or 1) or multivariate-t
// (robust_type 2) mixture with a Bayesian Lasso (Laplace) prior on the
// profile means, for the six variance-covariance models of robust_lpa().
//
// * Allocation uses the exact observed-data (FIML) marginal density of each
//   row under every profile (Gaussian or t), on the log scale.
// * Missing values are handled by data augmentation: after allocation,
//   each row's missing block is drawn from its conditional distribution
//   given the observed block (divided by the latent scale u_i for the t
//   model), so all complete-data conditionals below are exact.
// * robust_type 2 (multivariate t): the latent scale u_i is drawn from its
//   Gamma((nu + p_obs)/2, rate (nu + d_i)/2) full conditional, and enters
//   the mean/covariance conditionals as a precision weight -- a proper
//   Gibbs sampler for the t mixture (scale-mixture-of-normals). If
//   `estimate_nu`, nu is updated by a random-walk Metropolis step on log(nu)
//   under a Gamma(2, 0.1) prior, with the latent scales integrated out.
// * Models 4 and 5 (Sigma_g = D_g R_g D_g with correlation matrices R_g) use
//   exact slice-sampling updates of the scales and random-walk Metropolis
//   updates of the correlation matrices (uniform prior), with step sizes
//   adapted during burn-in only.
// * robust_type 1 (Huber): a heuristic that down-weights the sufficient
//   statistics of outlying observations with Huber weights computed from
//   their distance to the allocated profile; it does not target a
//   well-defined posterior, and is kept for backward compatibility.
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
Rcpp::List mcmc_chain_cpp(const arma::mat& X, int G, int model, int mcmc_iter, double prior_laplace,
                          int robust_type, double alpha, double nu_init, bool estimate_nu,
                          const Rcpp::List& init_mu, const Rcpp::List& init_sigma, const arma::vec& init_pi,
                          const Rcpp::List& pat_obs, const Rcpp::List& pat_rows) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  std::vector<arma::uvec> obs = as_uvec_list(pat_obs);
  std::vector<arma::uvec> rows = as_uvec_list(pat_rows);
  const int K = obs.size();
  std::vector<arma::uvec> mis(K);
  for (int k = 0; k < K; ++k) mis[k] = complement_idx(obs[k], p);

  const int dist = (robust_type == 2) ? 1 : 0;
  double nu = nu_init;

  std::vector<arma::rowvec> mu(G);
  std::vector<arma::mat> sigma(G);
  std::vector<arma::vec> tau(G);
  for (int g = 0; g < G; g++) {
    mu[g] = Rcpp::as<arma::rowvec>(init_mu[g]);
    sigma[g] = Rcpp::as<arma::mat>(init_sigma[g]);
    tau[g] = arma::ones<arma::vec>(p);
  }
  arma::vec pi_g = init_pi / arma::sum(init_pi);

  arma::vec cutoff(p + 1, arma::fill::zeros);
  for (int q = 1; q <= p; ++q) cutoff(q) = R::qchisq(1.0 - alpha, q, 1, 0);

  arma::mat Xc = X;
  arma::uvec alloc(n, arma::fill::zeros);
  arma::vec wgt(n, arma::fill::ones), uval(n, arma::fill::ones), dstar(n, arma::fill::zeros);
  arma::vec pobs_row(n, arma::fill::zeros);
  for (int k = 0; k < K; ++k) for (arma::uword ii = 0; ii < rows[k].n_elem; ++ii) pobs_row(rows[k](ii)) = obs[k].n_elem;

  Rcpp::List mu_chain(mcmc_iter), sigma_chain(mcmc_iter);
  arma::mat pi_chain(mcmc_iter, G);
  arma::vec nu_chain(mcmc_iter);
  int nu_accept = 0;

  std::vector<std::vector<PatCache>> cache(G, std::vector<PatCache>(K));

  // Models 4 and 5 are parameterized as Sigma_g = D_g R_g D_g with R_g a
  // correlation matrix (model 4: R_g = R shared; model 5: D_g = D shared).
  std::vector<arma::mat> Rcor(G);
  std::vector<arma::vec> dsc(G);
  for (int g = 0; g < G; g++) {
    dsc[g] = arma::sqrt(arma::clamp(sigma[g].diag(), 1e-10, arma::datum::inf));
    Rcor[g] = cov_to_cor(sigma[g]);
  }
  if (model == 4) for (int g = 1; g < G; g++) Rcor[g] = Rcor[0];
  if (model == 5) for (int g = 1; g < G; g++) dsc[g] = dsc[0];

  // Random-walk step sizes, adapted during burn-in only (first half of the
  // chain), so the retained draws come from a fixed, valid transition kernel.
  const int burnin = mcmc_iter / 2;
  double step_R = 0.1, step_nu = 0.3;
  int acc_R = 0, try_R = 0, acc_nu_win = 0, try_nu_win = 0;

  for (int iter = 0; iter < mcmc_iter; iter++) {
    Rcpp::checkUserInterrupt();

    // ---- per-profile, per-pattern factorizations -------------------------
    for (int g = 0; g < G; ++g) {
      for (int k = 0; k < K; ++k) {
        const arma::uvec& o = obs[k];
        const arma::uvec& m = mis[k];
        PatCache& pc = cache[g][k];
        if (o.n_elem > 0) {
          pc.L = chol_lower_safe(sigma[g].submat(o, o));
          pc.logdet = 2.0 * arma::sum(arma::log(pc.L.diag()));
        }
        if (m.n_elem > 0) {
          arma::mat Cmm;
          if (o.n_elem > 0) {
            arma::mat Smo = sigma[g].submat(m, o);
            arma::mat tmp = arma::solve(arma::trimatl(pc.L), Smo.t());
            pc.Bt = arma::solve(arma::trimatu(pc.L.t()), tmp);
            Cmm = sigma[g].submat(m, m) - Smo * pc.Bt;
          } else {
            Cmm = sigma[g];
          }
          pc.Lc = chol_lower_safe(0.5 * (Cmm + Cmm.t()));
        }
      }
    }

    // ---- allocation (u and missing values integrated out) ----------------
    for (int k = 0; k < K; ++k) {
      const arma::uvec& o = obs[k];
      const arma::uvec& r = rows[k];
      const int po = o.n_elem;
      for (arma::uword ii = 0; ii < r.n_elem; ++ii) {
        const arma::uword i = r(ii);
        if (po == 0) {
          alloc(i) = sample_class(pi_g);
          dstar(i) = 0.0;
          continue;
        }
        arma::vec xo(po);
        for (int j = 0; j < po; ++j) xo(j) = X(i, o(j));
        arma::vec logp(G), dvec(G);
        for (int g = 0; g < G; ++g) {
          arma::vec diff = xo - arma::vec(mu[g].elem(o));
          arma::vec y = arma::solve(arma::trimatl(cache[g][k].L), diff);
          double d = arma::dot(y, y);
          dvec(g) = d;
          logp(g) = std::log(std::max(pi_g(g), 1e-300)) +
            log_const(dist, po, cache[g][k].logdet, nu) + log_kernel(dist, po, d, nu);
        }
        arma::vec probs = arma::exp(logp - log_sum_exp_cpp(logp));
        int gs = sample_class(probs);
        alloc(i) = gs;
        dstar(i) = dvec(gs);
      }
    }

    // ---- degrees of freedom of the t model, u integrated out ----------------
    // Random-walk Metropolis on log(nu) targeting p(nu | allocation, mu,
    // Sigma, x_obs), i.e. the product of the (marginal) multivariate-t
    // densities of the observed data under a Gamma(2, 0.1) prior; the latent
    // scales are drawn afterwards given the new nu. Marginalizing u avoids
    // the strong nu-u dependence that makes the naive nu | u update mix
    // very slowly.
    if (robust_type == 2 && estimate_nu) {
      double prop = nu * std::exp(step_nu * R::rnorm(0.0, 1.0));
      try_nu_win++;
      if (prop >= 1.0 && prop <= 1000.0) {
        double log_r = log_marg_nu(prop, dstar, pobs_row) - log_marg_nu(nu, dstar, pobs_row) +
          (std::log(prop) - 0.1 * prop) - (std::log(nu) - 0.1 * nu) +   // Gamma(2, 0.1) prior
          std::log(prop) - std::log(nu);                              // log-scale Jacobian
        if (std::log(R::runif(0.0, 1.0)) < log_r) {
          nu = prop;
          nu_accept++;
          acc_nu_win++;
        }
      }
    }

    // ---- latent scales / robustness weights, data augmentation -------------
    for (int k = 0; k < K; ++k) {
      const arma::uvec& o = obs[k];
      const arma::uvec& m = mis[k];
      const arma::uvec& r = rows[k];
      const int po = o.n_elem;
      for (arma::uword ii = 0; ii < r.n_elem; ++ii) {
        const arma::uword i = r(ii);
        const int gs = alloc(i);
        const double d_star = dstar(i);

        double u = 1.0, w = 1.0;
        if (robust_type == 2) {
          u = R::rgamma(0.5 * (nu + po), 1.0 / (0.5 * (nu + d_star)));
          w = u;
        } else if (robust_type == 1 && po > 0 && d_star > cutoff(po)) {
          w = std::sqrt(cutoff(po) / d_star);
        }

        if (m.n_elem > 0) {
          arma::vec mean_m = arma::vec(mu[gs].elem(m));
          if (po > 0) {
            arma::vec xo(po);
            for (int j = 0; j < po; ++j) xo(j) = X(i, o(j));
            mean_m += cache[gs][k].Bt.t() * (xo - arma::vec(mu[gs].elem(o)));
          }
          arma::vec zn(m.n_elem, arma::fill::randn);
          arma::vec xm = mean_m + cache[gs][k].Lc * zn / std::sqrt(u);
          for (arma::uword j = 0; j < m.n_elem; ++j) Xc(i, m(j)) = xm(j);
        }
        wgt(i) = w;
        uval(i) = u;
      }
    }

    // ---- profile means (Bayesian Lasso) and scatter matrices --------------
    arma::vec n_counts(G, arma::fill::zeros), df_g(G, arma::fill::zeros);
    std::vector<arma::mat> S(G, arma::mat(p, p, arma::fill::zeros));

    for (int g = 0; g < G; g++) {
      arma::uvec idx = arma::find(alloc == static_cast<arma::uword>(g));
      const double n_g = idx.n_elem;
      n_counts(g) = n_g;
      if (n_g == 0) continue;   // empty profile: keep its current mean

      arma::mat Xg = Xc.rows(idx);
      arma::vec wg = wgt.elem(idx);
      double sw = arma::sum(wg);
      arma::vec sx = Xg.t() * wg;

      arma::mat Sigma_inv = inv_spd_safe(sigma[g]);
      arma::vec D_tau_inv = 1.0 / arma::clamp(tau[g], 1e-8, arma::datum::inf);
      arma::mat V_post = inv_spd_safe(sw * Sigma_inv + arma::diagmat(D_tau_inv));
      arma::vec m_post = V_post * (Sigma_inv * sx);
      arma::vec zn(p, arma::fill::randn);
      arma::rowvec new_mu = (m_post + chol_lower_safe(V_post) * zn).t();
      mu[g] = new_mu;

      for (int j = 0; j < p; j++) {
        double mu_j_abs = std::max(std::abs(new_mu(j)), 1e-6);
        tau[g](j) = 1.0 / rinvgauss(prior_laplace / mu_j_abs, prior_laplace * prior_laplace);
      }

      arma::mat Dm = Xg.each_row() - new_mu;
      S[g] = Dm.t() * (Dm.each_col() % wg);
      // Degrees of freedom of the covariance update: the raw count for the
      // Gaussian and t models (for t, u_i enters only through the scatter
      // matrix, exactly as in the scale-mixture representation); the
      // effective Huber-weighted size for the heuristic Huber sampler.
      df_g(g) = (robust_type == 1) ? sw : n_g;
    }

    // ---- mixing proportions: Dirichlet(n_g + 1) full conditional -----------
    {
      arma::vec gamma_draws(G);
      for (int g = 0; g < G; g++) gamma_draws(g) = R::rgamma(n_counts(g) + 1.0, 1.0);
      pi_g = gamma_draws / arma::sum(gamma_draws);
    }

    // ---- covariance matrices -------------------------------------------------
    const double total_df = arma::sum(df_g);
    arma::mat S_pool(p, p, arma::fill::zeros);
    for (int g = 0; g < G; g++) S_pool += S[g];

    if (model == 1) {
      // Variable-specific variances shared across profiles (as in the EM
      // engine): sigma2_j ~ IG(1 + total_df / 2, 1 + sum_g S_g[j, j] / 2).
      arma::vec pooled_var(p);
      for (int j = 0; j < p; j++) pooled_var(j) = rinvgamma_cpp(0.5 * total_df + 1.0, 0.5 * S_pool(j, j) + 1.0);
      for (int g = 0; g < G; g++) sigma[g] = arma::diagmat(pooled_var);
    } else if (model == 2) {
      for (int g = 0; g < G; g++) {
        arma::mat new_sig(p, p, arma::fill::zeros);
        for (int j = 0; j < p; j++) new_sig(j, j) = rinvgamma_cpp(0.5 * df_g(g) + 1.0, 0.5 * S[g](j, j) + 1.0);
        sigma[g] = new_sig;
      }
    } else if (model == 3) {
      arma::mat pooled_sig = rinvwishart(total_df + p + 1.0, S_pool + arma::eye<arma::mat>(p, p));
      for (int g = 0; g < G; g++) sigma[g] = pooled_sig;
    } else if (model == 4) {
      // Sigma_g = D_g R D_g: exact slice-sampling updates of each profile's
      // scales given R, then Metropolis updates of the shared correlation
      // matrix R given the scales (uniform prior on correlation matrices).
      arma::mat Rinv = inv_spd_safe(Rcor[0]);
      arma::mat Tsum(p, p, arma::fill::zeros);
      for (int g = 0; g < G; g++) {
        arma::vec a = draw_inverse_scales(Rinv % S[g], 1.0 / dsc[g], df_g(g));
        dsc[g] = 1.0 / a;
        Tsum += S[g] / (dsc[g] * dsc[g].t());
      }
      if (p > 1) Rcor[0] = mh_correlation(Rcor[0], Tsum, total_df, step_R, acc_R, try_R);
      for (int g = 0; g < G; g++) {
        Rcor[g] = Rcor[0];
        sigma[g] = Rcor[0] % (dsc[g] * dsc[g].t());
      }
    } else if (model == 5) {
      // Sigma_g = D R_g D: exact slice-sampling update of the shared scales
      // given the correlation matrices, then Metropolis updates of each
      // profile's correlation matrix given the scales.
      arma::mat M(p, p, arma::fill::zeros);
      for (int g = 0; g < G; g++) M += inv_spd_safe(Rcor[g]) % S[g];
      arma::vec a = draw_inverse_scales(M, 1.0 / dsc[0], total_df);
      arma::vec d = 1.0 / a;
      for (int g = 0; g < G; g++) {
        dsc[g] = d;
        if (p > 1) Rcor[g] = mh_correlation(Rcor[g], S[g] / (d * d.t()), df_g(g), step_R, acc_R, try_R);
        sigma[g] = Rcor[g] % (d * d.t());
      }
    } else if (model == 6) {
      for (int g = 0; g < G; g++) sigma[g] = rinvwishart(df_g(g) + p + 1.0, S[g] + arma::eye<arma::mat>(p, p));
    }

    // ---- adapt random-walk step sizes during burn-in ---------------------------
    if (iter < burnin && (iter + 1) % 50 == 0) {
      if (try_R > 0) {
        double rate = static_cast<double>(acc_R) / try_R;
        if (rate > 0.35) step_R *= 1.25; else if (rate < 0.2) step_R /= 1.25;
        step_R = std::min(std::max(step_R, 1e-3), 0.5);
      }
      if (try_nu_win > 0) {
        double rate = static_cast<double>(acc_nu_win) / try_nu_win;
        if (rate > 0.45) step_nu *= 1.25; else if (rate < 0.3) step_nu /= 1.25;
        step_nu = std::min(std::max(step_nu, 0.02), 2.0);
      }
      acc_R = try_R = acc_nu_win = try_nu_win = 0;
    }

    // ---- store ----------------------------------------------------------------
    Rcpp::List mu_iter(G), sigma_iter(G);
    for (int g = 0; g < G; g++) {
      mu_iter[g] = mu[g];
      sigma_iter[g] = sigma[g];
    }
    mu_chain[iter] = mu_iter;
    sigma_chain[iter] = sigma_iter;
    pi_chain.row(iter) = pi_g.t();
    nu_chain(iter) = nu;
  }

  return Rcpp::List::create(
    Rcpp::Named("mu_chain") = mu_chain,
    Rcpp::Named("sigma_chain") = sigma_chain,
    Rcpp::Named("pi_chain") = pi_chain,
    Rcpp::Named("nu_chain") = nu_chain,
    Rcpp::Named("nu_acceptance") = (mcmc_iter > 0) ? static_cast<double>(nu_accept) / mcmc_iter : 0.0
  );
}
