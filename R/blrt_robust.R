#' Bootstrapped Likelihood Ratio Test for Robust LPA
#'
#' Compares a robust LPA model with \code{G} profiles against a null model
#' with \code{G - 1} profiles using parametric bootstrapping (Nylund et al.,
#' 2007): the null model is fit to the observed data, data are simulated
#' from it, and both the null and alternative models are refit to each
#' simulated dataset to build a reference distribution for the likelihood
#' ratio test statistic under \code{H0}. Supports FIML simulation conditions
#' (the missingness pattern of the observed data is replicated in every
#' simulated dataset).
#'
#' @param data A matrix or data.frame.
#' @param G The number of profiles for the alternative hypothesis (compared against \code{G - 1}).
#' @param model An integer (1 to 6) specifying the variance-covariance parameterization (see \code{\link{robust_lpa}}).
#' @param engine String, either \code{"EM"} (default) or \code{"MCMC"}; passed
#'   through to every internal call to \code{\link{robust_lpa}} so that the
#'   observed and bootstrap-refit models use the same estimation engine.
#' @param n_samples Number of bootstrap samples. Default is 50 for speed; 200+ is recommended for publications.
#' @param n_starts Number of starts for the EM algorithm execution (ignored when \code{engine = "MCMC"}).
#' @param cores Integer, number of CPU cores to use to run the \code{n_samples}
#'   bootstrap replicates in parallel (default \code{1}, sequential). Each
#'   replicate fits two independent \code{\link{robust_lpa}} models (null and
#'   alternative) on its own simulated dataset, so replicates are
#'   "embarrassingly parallel". Uses the same \code{parallel::mclapply()} /
#'   \code{parallel::makeCluster()} PSOCK backend as \code{\link{robust_lpa}}'s
#'   own \code{cores} argument (see there for details); falls back to
#'   sequential execution with a \code{warning()} if the \pkg{parallel}
#'   package is unavailable. If you also pass \code{cores} through \code{...}
#'   to \code{\link{robust_lpa}} (to parallelize each replicate's EM restarts
#'   / MCMC chains too), keep the product of the two \code{cores} values at or
#'   below your machine's core count to avoid oversubscription; for most uses
#'   it is simplest to parallelize only at this (bootstrap-replicate) level.
#' @param ... Additional arguments passed on to \code{\link{robust_lpa}}
#'   (e.g. \code{max_iter}, \code{tol}, \code{mcmc_iter}, \code{n_chains},
#'   \code{prior_laplace}, \code{robust}, \code{alpha}) -- for instance, pass
#'   \code{robust = FALSE} here to run the BLRT with classical (non-robust)
#'   estimation throughout, for either engine. Note that with
#'   \code{engine = "MCMC"} every one of the \code{2 * (n_samples + 1)}
#'   internal \code{\link{robust_lpa}} calls this function makes will run
#'   \code{n_chains} chains each; consider passing a smaller \code{n_chains}
#'   and/or \code{mcmc_iter} than the \code{\link{robust_lpa}} defaults, and/or
#'   using \code{cores > 1} above, to keep the BLRT's bootstrap loop tractable.
#' @return A list containing:
#'   \describe{
#'     \item{LRT_Observed}{The observed likelihood ratio test statistic (non-negative).}
#'     \item{Bootstrap_LRTs}{Numeric vector of the successfully-fit bootstrap replicates.}
#'     \item{p_value}{The empirical p-value, computed with the standard
#'       "+1" small-sample correction (\code{(sum(Bootstrap_LRTs >= LRT_Observed) + 1) / (length(Bootstrap_LRTs) + 1)}),
#'       which avoids reporting an (impossible) exact p-value of 0 from a finite bootstrap.}
#'     \item{Bootstrap_Failures}{Integer, how many of the \code{n_samples}
#'       bootstrap replicates failed to converge and were excluded from
#'       \code{Bootstrap_LRTs} / \code{p_value}.}
#'   }
#' @references
#'   Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
#'   number of classes in latent class analysis and growth mixture modeling:
#'   A Monte Carlo simulation study. \emph{Structural Equation Modeling},
#'   14(4), 535-569. \doi{10.1080/10705510701575396}
#' @examples
#' # Fast demonstration of the robust BLRT: is a 2nd profile justified over 1?
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' blrt_res <- suppressWarnings(blrt_robust(x, G = 2, model = 1, n_samples = 2, n_starts = 3))
#' # Print the summary of the results
#' blrt_res
#' @export
blrt_robust <- function(data, G, model = 6, engine = "EM", n_samples = 50, n_starts = 2, cores = 1, ...) {

  # ---- input validation -----------------------------------------------
  if (!is.numeric(G) || length(G) != 1 || G <= 1 || G != round(G)) {
    stop("`G` must be a single integer greater than 1 to perform the Bootstrapped Likelihood Ratio Test.")
  }
  if (!(length(model) == 1 && model %in% 1:6)) {
    stop("`model` must be a single integer between 1 and 6.")
  }
  if (!(length(engine) == 1 && engine %in% c("EM", "MCMC"))) {
    stop("`engine` must be either 'EM' or 'MCMC'.")
  }
  if (!is.numeric(n_samples) || n_samples < 1) stop("`n_samples` must be at least 1.")
  if (!is.numeric(n_starts) || n_starts < 1) stop("`n_starts` must be at least 1.")
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive integer.")
  }

  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)
  dots <- list(...)

  message(paste0("Executing Robust BLRT: Comparing ", G - 1, " vs ", G, " Profiles (Model ", model, ")..."))

  # Estimation of empirical models on real data
  mod_null <- do.call(robust_lpa, c(list(data = X, G = G - 1, model = model, engine = engine, n_starts = n_starts), dots))
  mod_alt  <- do.call(robust_lpa, c(list(data = X, G = G,     model = model, engine = engine, n_starts = n_starts), dots))

  if (is.null(mod_null) || is.null(mod_alt)) {
    stop("One of the baseline models failed to converge. Unable to compute BLRT.")
  }

  # Calculation of the observed empirical log-likelihood ratio
  lrt_obs <- -2 * (mod_null$fit$LogLik - mod_alt$fit$LogLik)
  if (lrt_obs < 0) lrt_obs <- 0

  # Parametric stochastic generator based on null model parameters
  simulate_null_mixture <- function(n, p, props, means, covs) {
    sim_X <- matrix(NA, nrow = n, ncol = p)
    assigned_classes <- sample(1:(G - 1), size = n, replace = TRUE, prob = props)

    for (i in 1:n) {
      c_i <- assigned_classes[i]
      # Cholesky decomposition to generate stable multivariate distributions
      z_norm <- rnorm(p)
      L <- chol(covs[[c_i]] + diag(1e-6, p))
      sim_X[i, ] <- means[[c_i]] + as.vector(t(L) %*% z_norm)
    }
    return(sim_X)
  }

  # One bootstrap replicate: simulate data under H0 from the observed null
  # model, then refit both the null and alternative models on it. Wrapped in
  # its own tryCatch (rather than relying on the caller) so that an
  # unexpected error on a parallel worker returns NA for that replicate
  # instead of aborting the whole `.run_parallel()` dispatch.
  run_one_bootstrap <- function(b) {
    if (cores == 1) message(paste0("Bootstrap Sample ", b, " / ", n_samples))

    tryCatch({
      sim_data <- simulate_null_mixture(n, p, mod_null$proportions, mod_null$means, mod_null$covariances)

      # If the initial data contained NA, replicate the same missingness structure for FIML
      if (any(is.na(X))) {
        sim_data[is.na(X)] <- NA
      }

      # Fit of models on simulated data
      fit_b_null <- tryCatch(do.call(robust_lpa, c(list(data = sim_data, G = G - 1, model = model, engine = engine, n_starts = n_starts), dots)), error = function(e) NULL)
      fit_b_alt  <- tryCatch(do.call(robust_lpa, c(list(data = sim_data, G = G,     model = model, engine = engine, n_starts = n_starts), dots)), error = function(e) NULL)

      if (!is.null(fit_b_null) && !is.null(fit_b_alt)) {
        val <- -2 * (fit_b_null$fit$LogLik - fit_b_alt$fit$LogLik)
        if (val < 0) val <- 0
        return(val)
      }
      NA_real_
    }, error = function(e) NA_real_)
  }

  # Bootstrap resampling cycle (sequential if `cores == 1`, parallel otherwise;
  # see `cores`). Order is preserved regardless of backend.
  boot_results <- .run_parallel(
    cores, n_samples, run_one_bootstrap,
    export_vars = c("X", "n", "p", "G", "model", "engine", "n_starts", "n_samples", "mod_null", "dots", "simulate_null_mixture"),
    export_env = environment()
  )
  if (cores > 1) {
    message(sprintf("Completed %d bootstrap replicate(s) (parallel, cores = %d).", n_samples, cores))
  }
  lrt_boot <- vapply(boot_results, function(v) if (is.null(v)) NA_real_ else v, numeric(1))

  # Cleaning vectors from potential internal convergence failures
  n_failed <- sum(is.na(lrt_boot))
  lrt_boot <- lrt_boot[!is.na(lrt_boot)]
  
  if (length(lrt_boot) == 0) {
    stop("All ", n_samples, " bootstrap replicates failed to converge. Unable to compute a BLRT p-value. ",
         "Try increasing `n_starts`, simplifying `model`, or inspecting the null model fit for near-degenerate profiles.")
  }
  if (n_failed > 0) {
    warning(sprintf("%d of %d bootstrap replicates failed to converge and were excluded from the p-value calculation.", n_failed, n_samples))
  }
  
  # Empirical p-value with the standard "+1" small-sample correction (see
  # e.g. Davison & Hinkley, 1997), which avoids ever reporting an exact
  # p-value of 0 -- a value that is not attainable from a finite bootstrap.
  p_value <- (sum(lrt_boot >= lrt_obs) + 1) / (length(lrt_boot) + 1)
  message("BLRT Execution Successfully Completed.")
  
  return(list(
    LRT_Observed = lrt_obs,
    Bootstrap_LRTs = lrt_boot,
    p_value = p_value,
    Bootstrap_Failures = n_failed
  ))
}