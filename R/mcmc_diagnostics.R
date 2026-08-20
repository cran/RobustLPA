#' Compute Gelman-Rubin R-hat and Effective Sample Size for an MCMC Fit
#'
#' Computes the classic multi-chain Gelman-Rubin potential scale reduction
#' statistic (\eqn{\hat{R}}; Gelman & Rubin, 1992) and the classic
#' (autocorrelation/spectral-density-based) effective sample size for every
#' scalar parameter of a \code{\link{robust_lpa}} MCMC fit, via
#' \code{coda::gelman.diag()} and \code{coda::effectiveSize()}. Burn-in (the
#' first half of each chain, as used elsewhere in \code{\link{robust_lpa}})
#' is discarded before computing either diagnostic; both are computed on the
#' same \code{[iterations, chains, parameters]} array built by
#' \code{\link{.reshape_mcmc_draws}}, so parameter names and ordering are
#' always consistent with \code{\link{plot_mcmc_chains}}.
#'
#' \eqn{\hat{R}} compares the between-chain and within-chain variance of each
#' parameter; values close to 1 indicate the chains have converged to the
#' same target distribution, while values above about 1.1 are the classic
#' rule-of-thumb threshold for suspecting non-convergence (Gelman & Rubin,
#' 1992; Gelman et al., 2013, \emph{Bayesian Data Analysis}). It is a
#' between-chain diagnostic by construction and is therefore undefined
#' (returned as \code{NA}) when only a single chain is available.
#'
#' Effective sample size (ESS) estimates how many independent draws the
#' (autocorrelated) MCMC draws are worth; a low ESS relative to the total
#' number of post-burn-in draws indicates high autocorrelation and less
#' precise posterior summaries.
#'
#' @param chains A list of length \code{n_chains}, each element as returned
#'   by \code{robust_mcmc_cpp()} (see \code{\link{.reshape_mcmc_draws}}).
#' @param mcmc_iter Integer, the total number of iterations per chain.
#' @param burnin Integer, the number of leading iterations per chain to
#'   discard as burn-in before computing diagnostics.
#' @return A data.frame with columns \code{Parameter}, \code{Rhat}, and
#'   \code{ESS}, one row per scalar parameter (see
#'   \code{\link{.mcmc_param_names}}); or \code{NULL} (with a
#'   \code{warning()}) if the \pkg{coda} package is not installed, or if
#'   fewer than 2 post-burn-in iterations are available per chain.
#' @references
#'   Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#'   using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'   \doi{10.1214/ss/1177011136}
#' @keywords internal
#' @noRd
.compute_mcmc_diagnostics <- function(chains, mcmc_iter, burnin) {
  valid_iters <- (burnin + 1):mcmc_iter
  n_valid <- length(valid_iters)
  n_chains <- length(chains)

  if (n_valid < 2) {
    warning(
      "Fewer than 2 post-burn-in MCMC iterations are available; convergence diagnostics ",
      "(R-hat, ESS) cannot be computed. Increase `mcmc_iter`."
    )
    return(NULL)
  }

  if (!requireNamespace("coda", quietly = TRUE)) {
    warning(
      "The 'coda' package is required to compute MCMC convergence diagnostics (Gelman-Rubin ",
      "R-hat, effective sample size). Install it with install.packages('coda'); returning NULL ",
      "for `mcmc_diagnostics`."
    )
    return(NULL)
  }

  draws_array <- .reshape_mcmc_draws(chains, iters = valid_iters)
  param_names <- dimnames(draws_array)[[3]]

  mcmc_list <- coda::mcmc.list(
    lapply(seq_len(n_chains), function(k) coda::mcmc(draws_array[, k, , drop = TRUE]))
  )

  ess <- tryCatch({
    es <- coda::effectiveSize(mcmc_list)
    as.numeric(es[param_names])
  }, error = function(e) {
    warning("Failed to compute effective sample size via 'coda': ", conditionMessage(e))
    rep(NA_real_, length(param_names))
  })

  if (n_chains < 2) {
    rhat <- rep(NA_real_, length(param_names))
  } else {
    rhat <- tryCatch({
      gd <- coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
      point_est <- gd$psrf[, "Point est."]
      names(point_est) <- rownames(gd$psrf)
      as.numeric(point_est[param_names])
    }, error = function(e) {
      warning(
        "coda::gelman.diag() failed (", conditionMessage(e), "); falling back to a manual ",
        "Gelman-Rubin calculation."
      )
      as.numeric(.rhat_manual(draws_array)[param_names])
    })
  }

  data.frame(
    Parameter = param_names,
    Rhat = rhat,
    ESS = ess,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Manual Fallback Gelman-Rubin R-hat Calculation
#'
#' A direct implementation of the classic (non-rank-normalized, non-split)
#' Gelman-Rubin potential scale reduction statistic (Gelman & Rubin, 1992),
#' used as a fallback by \code{\link{.compute_mcmc_diagnostics}} if
#' \code{coda::gelman.diag()} errors for a given draws array (e.g. because
#' \pkg{coda} declines to invert a near-singular multi-parameter covariance
#' -- which cannot happen here since diagnostics are requested with
#' \code{multivariate = FALSE}, but this remains a defensive fallback for any
#' other \pkg{coda}-internal failure).
#'
#' @param draws_array A \code{[iterations, chains, parameters]} numeric array
#'   as returned by \code{\link{.reshape_mcmc_draws}} (already restricted to
#'   post-burn-in iterations).
#' @return A named numeric vector (one value per parameter). \code{NA} for
#'   any parameter whose within-chain variance is numerically zero (R-hat is
#'   undefined for a parameter that never moves).
#' @references
#'   Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#'   using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'   \doi{10.1214/ss/1177011136}
#' @keywords internal
#' @noRd
.rhat_manual <- function(draws_array) {
  n_iter <- dim(draws_array)[1]
  n_chains <- dim(draws_array)[2]
  n_params <- dim(draws_array)[3]
  param_names <- dimnames(draws_array)[[3]]

  rhat <- rep(NA_real_, n_params)
  names(rhat) <- param_names

  for (pidx in seq_len(n_params)) {
    chain_means <- numeric(n_chains)
    chain_vars <- numeric(n_chains)
    for (c_idx in seq_len(n_chains)) {
      x <- draws_array[, c_idx, pidx]
      chain_means[c_idx] <- mean(x)
      chain_vars[c_idx] <- stats::var(x)
    }
    grand_mean <- mean(chain_means)
    B <- (n_iter / (n_chains - 1)) * sum((chain_means - grand_mean)^2)
    W <- mean(chain_vars)

    if (!is.finite(W) || W < 1e-12) next

    var_hat <- ((n_iter - 1) / n_iter) * W + B / n_iter
    rhat[pidx] <- sqrt(var_hat / W)
  }

  rhat
}
