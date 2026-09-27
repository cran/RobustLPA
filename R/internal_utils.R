#' Null-Default Operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Missingness Patterns of a Data Matrix
#'
#' Groups the rows of \code{X} by their pattern of observed variables, so the
#' C++ engine can factorize each pattern's sub-covariance once per profile
#' instead of once per row.
#'
#' @param X A numeric matrix (\code{NA} allowed).
#' @return A list with \code{obs} (list of 0-based observed-column indices,
#'   one per pattern), \code{rows} (list of 0-based row indices, one per
#'   pattern), and \code{pobs} (integer vector, number of observed variables
#'   per row).
#' @keywords internal
#' @noRd
.missing_patterns <- function(X) {
  obs_mat <- !is.na(X)
  if (ncol(X) == 0L) stop("`X` must have at least one column.")
  key <- apply(obs_mat, 1, function(r) paste(as.integer(r), collapse = ""))
  ukeys <- unique(key)
  pat_id <- match(key, ukeys)
  obs <- lapply(ukeys, function(k) as.integer(which(strsplit(k, "", fixed = TRUE)[[1]] == "1") - 1L))
  rows <- lapply(seq_along(ukeys), function(k) as.integer(which(pat_id == k) - 1L))
  list(obs = obs, rows = rows, pobs = as.integer(rowSums(obs_mat)))
}

#' Row-Wise Log-Sum-Exp
#'
#' Numerically stable \code{log(rowSums(exp(M)))}, used to turn per-profile
#' log-densities into log-likelihood contributions and posterior
#' probabilities without underflow (so even extreme outliers get valid
#' posterior probabilities that sum to 1).
#'
#' @param M A numeric matrix.
#' @return A numeric vector of length \code{nrow(M)}.
#' @keywords internal
#' @noRd
.row_logsumexp <- function(M) {
  if (ncol(M) == 1L) return(M[, 1])
  m <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  m[!is.finite(m)] <- 0
  m + log(rowSums(exp(M - m)))
}

#' Posterior Probabilities from Per-Profile Log-Densities
#'
#' @param logf An \code{n x G} matrix of \code{log(pi_g) + log f_g(x_i)}.
#' @return A list with \code{z} (posterior probabilities, rows summing to 1)
#'   and \code{loglik} (the observed-data log-likelihood).
#' @keywords internal
#' @noRd
.posterior_from_logf <- function(logf) {
  lse <- .row_logsumexp(logf)
  z <- exp(logf - lse)
  bad <- !is.finite(lse) | !is.finite(rowSums(z))
  if (any(bad)) z[bad, ] <- 1 / ncol(logf)
  list(z = z, loglik = sum(lse[is.finite(lse)]))
}

#' Map the User-Facing Robustness Arguments to an Internal Method Name
#' @keywords internal
#' @noRd
.robust_method_name <- function(robust, robust_method) {
  if (!isTRUE(robust)) "none" else robust_method
}

#' Robustness Weights from Squared Mahalanobis Distances
#'
#' \itemize{
#'   \item \code{"none"}: all weights 1 (classical Gaussian estimation).
#'   \item \code{"huber"}: 1 inside the \code{1 - alpha} chi-squared cutoff
#'     for the row's number of observed variables, \code{sqrt(cutoff / d)}
#'     beyond it.
#'   \item \code{"t"}: the E-step expectation of the multivariate-t latent
#'     scale, \code{(nu + p_obs) / (nu + d)} (Liu & Rubin, 1995).
#' }
#'
#' @param maha Squared Mahalanobis distances on the observed entries.
#' @param pobs Number of observed variables per row.
#' @param method One of \code{"none"}, \code{"huber"}, \code{"t"}.
#' @param alpha Huber significance level.
#' @param nu Degrees of freedom of the t model.
#' @return A numeric vector of weights.
#' @keywords internal
#' @noRd
.robust_weights <- function(maha, pobs, method, alpha = 0.05, nu = 4) {
  if (method == "none") return(rep(1, length(maha)))
  if (method == "t") return((nu + pobs) / (nu + maha))
  cutoff <- stats::qchisq(1 - alpha, df = pmax(pobs, 1))
  w <- ifelse(maha > cutoff, sqrt(cutoff / pmax(maha, 1e-300)), 1)
  w[pobs == 0] <- 1
  w
}

#' ECM Update of a Common Degrees-of-Freedom Parameter (t Mixture)
#'
#' Solves the conditional-maximization equation for \code{nu} of a
#' multivariate-t mixture with a common \code{nu} across profiles
#' (McLachlan & Peel, 2000, Section 7.5; with row-specific observed
#' dimensions under missing data, Liu & Rubin, 1995), using the E-step
#' quantities computed at the previous value \code{nu_old}. The solution is
#' bounded to \code{[1, 200]}; 200 effectively corresponds to a Gaussian.
#'
#' @param z An \code{n x G} matrix of posterior probabilities.
#' @param maha_mat An \code{n x G} matrix of squared Mahalanobis distances.
#' @param pobs Integer vector, observed variables per row.
#' @param nu_old Current value of \code{nu}.
#' @return The updated \code{nu}.
#' @keywords internal
#' @noRd
.update_nu_ecm <- function(z, maha_mat, pobs, nu_old) {
  keep <- pobs > 0
  if (!any(keep)) return(nu_old)
  zk <- z[keep, , drop = FALSE]
  dk <- maha_mat[keep, , drop = FALSE]
  pk <- pobs[keep]
  u <- (nu_old + pk) / (nu_old + dk)
  const <- 1 + mean(rowSums(zk * (log(u) - u))) +
    mean(digamma((nu_old + pk) / 2) - log((nu_old + pk) / 2))
  f <- function(v) -digamma(v / 2) + log(v / 2) + const
  lo <- 1
  hi <- 200
  if (f(hi) > 0) return(hi)
  if (f(lo) < 0) return(lo)
  stats::uniroot(f, c(lo, hi), tol = 1e-8)$root
}

#' Per-Profile Log-Densities for a Fitted Model at Given Parameters
#'
#' @param X Numeric matrix (\code{NA} allowed).
#' @param pat Result of \code{.missing_patterns(X)}.
#' @param means,covariances Lists of length \code{G}.
#' @param proportions Numeric vector of length \code{G}.
#' @param method \code{"none"}, \code{"huber"} or \code{"t"}.
#' @param nu Degrees of freedom (t only).
#' @return A list with \code{logf} (\code{n x G} matrix of
#'   \code{log(pi_g) + log f_g}) and \code{maha} (\code{n x G}).
#' @keywords internal
#' @noRd
.mixture_logf <- function(X, pat, means, covariances, proportions, method, nu = NULL) {
  G <- length(means)
  n <- nrow(X)
  dist <- if (identical(method, "t")) 1L else 0L
  nu_val <- if (dist == 1L) nu else 1
  logf <- matrix(0, n, G)
  maha <- matrix(0, n, G)
  for (g in seq_len(G)) {
    es <- class_estep_cpp(X, as.numeric(means[[g]]), covariances[[g]], pat$obs, pat$rows, dist, nu_val)
    logf[, g] <- es$logdens + log(max(proportions[g], 1e-300))
    maha[, g] <- es$maha
  }
  list(logf = logf, maha = maha)
}

#' Log-Likelihood of New Data Under a Fitted Model
#'
#' Observed-data (FIML) log-likelihood of \code{newX} under the parameters of
#' a fitted \code{robust_lpa} object (Gaussian mixture for classical and
#' Huber fits, multivariate-t mixture for \code{robust_method = "t"} fits).
#' Used for cross-validated LASSO tuning in
#' \code{\link{estimate_profiles_robust}}.
#'
#' @param fit A \code{robust_lpa} object.
#' @param newX A numeric matrix with the same columns as the fitted data.
#' @return A single number.
#' @keywords internal
#' @noRd
.loglik_newdata <- function(fit, newX) {
  newX <- as.matrix(newX)
  storage.mode(newX) <- "double"
  pat <- .missing_patterns(newX)
  lf <- .mixture_logf(newX, pat, fit$means, fit$covariances, fit$proportions,
                      fit$robust_method %||% "none", fit$nu)
  .posterior_from_logf(lf$logf)$loglik
}

#' Simulate From a Fitted Gaussian or Multivariate-t Mixture
#'
#' @param n Number of observations.
#' @param props,means,covs Mixing proportions, mean vectors and covariance
#'   (or t scale) matrices.
#' @param method \code{"t"} simulates multivariate-t components with
#'   \code{nu} degrees of freedom; anything else simulates Gaussian ones.
#' @param nu Degrees of freedom (t only).
#' @return An \code{n x p} numeric matrix.
#' @keywords internal
#' @noRd
.simulate_mixture <- function(n, props, means, covs, method = "none", nu = NULL) {
  G <- length(props)
  p <- length(means[[1]])
  cls <- if (G == 1) rep(1L, n) else sample.int(G, size = n, replace = TRUE, prob = props)
  out <- matrix(NA_real_, nrow = n, ncol = p)
  for (g in seq_len(G)) {
    idx <- which(cls == g)
    if (length(idx) == 0) next
    L <- chol(.force_pd(covs[[g]]))
    Z <- matrix(stats::rnorm(length(idx) * p), nrow = length(idx), ncol = p) %*% L
    if (identical(method, "t")) {
      Z <- Z / sqrt(stats::rchisq(length(idx), df = nu) / nu)
    }
    out[idx, ] <- sweep(Z, 2, as.numeric(means[[g]]), "+")
  }
  out
}

#' ECME Update of a Common Degrees-of-Freedom Parameter (t Mixture)
#'
#' Maximizes the observed-data (FIML) log-likelihood of the t mixture over a
#' common \code{nu}, holding the other parameters at their just-updated
#' values (the ECME variant of Liu & Rubin, 1994/1995). This keeps the
#' algorithm monotone and converges much faster than the ECM update of
#' \code{nu} when the likelihood is flat in \code{nu}. The search is on
#' \code{log(nu)} over \code{[1, 200]}; the current value is kept if it is
#' at least as good as the optimizer's solution.
#'
#' @param X,pat Data and its missingness patterns.
#' @param mu,sigma,pi_g Current profile parameters.
#' @param nu_cur Current \code{nu}.
#' @return The updated \code{nu}.
#' @keywords internal
#' @noRd
.update_nu_ecme <- function(X, pat, mu, sigma, pi_g, nu_cur) {
  obj <- function(lv) {
    .posterior_from_logf(.mixture_logf(X, pat, mu, sigma, pi_g, "t", exp(lv))$logf)$loglik
  }
  opt <- stats::optimize(obj, interval = log(c(1, 200)), maximum = TRUE, tol = 1e-3)
  if (obj(log(nu_cur)) >= opt$objective) nu_cur else exp(opt$maximum)
}
