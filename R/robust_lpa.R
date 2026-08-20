#' Fit a Single Robust Latent Profile Analysis Model
#'
#' Estimates a Latent Profile Analysis (Gaussian mixture) model that is
#' robust to multivariate outliers (via Huber down-weighting) and to missing
#' data (via a Full Information Maximum Likelihood, FIML, available-case
#' treatment), using either an EM or an MCMC (Bayesian Lasso) engine. Both
#' engines share the same \code{robust}/\code{alpha} Huber down-weighting
#' mechanism, on by default (see the "Robust estimation" section below). The
#' MCMC engine runs \code{n_chains} independent chains (4 by default) and
#' reports classic multi-chain convergence diagnostics (Gelman-Rubin
#' \eqn{\hat{R}} and effective sample size) in \code{$mcmc_diagnostics}. Set
#' \code{cores > 1} to run the EM engine's random restarts, or the MCMC
#' engine's chains, in parallel.
#'
#' @param data A matrix or data.frame of observations (numeric columns only;
#'   \code{NA} is allowed and triggers the FIML code path).
#' @param G The number of latent profiles to extract (a single positive integer).
#' @param model An integer (1 to 6) specifying the variance-covariance
#'   parameterization, following the same numbering convention as tidyLPA /
#'   mclust:
#'   \describe{
#'     \item{1}{Equal variances across profiles, covariances fixed to 0 (diagonal,
#'       shared across profiles -- each variable keeps its own variance level,
#'       constrained to be the same in every profile; this is \emph{not} a
#'       single isotropic/spherical variance shared across variables too,
#'       despite that being a common shorthand for this model elsewhere).}
#'     \item{2}{Varying variances across profiles, covariances fixed to 0 (diagonal, profile-specific).}
#'     \item{3}{Equal variances and equal covariances across profiles (one shared full covariance matrix).}
#'     \item{4}{Varying variances, equal covariance \emph{structure} (shared correlation matrix, profile-specific variances).}
#'     \item{5}{Equal variances, varying covariances (shared variances, profile-specific correlation matrices).}
#'     \item{6}{Fully unconstrained: each profile has its own variances and covariances (default).}
#'   }
#' @param engine String. Either \code{"EM"} (default) or \code{"MCMC"}.
#' @param max_iter Maximum number of EM iterations. Ignored when \code{engine = "MCMC"}.
#' @param tol Tolerance for EM convergence (on the observed-data log-likelihood). Ignored when \code{engine = "MCMC"}.
#' @param n_starts Number of random EM initializations; the fit with the
#'   highest log-likelihood across starts is returned. Ignored when \code{engine = "MCMC"}.
#' @param lambda Non-negative soft-thresholding (LASSO-type) penalty applied
#'   to the profile means, via direct per-coordinate soft-thresholding of the
#'   (Huber- and posterior-probability-weighted) mean at every M-step -- see
#'   \code{\link{robust_m_step}}/\code{robust_update_cpp()}. This shrinks each
#'   mean component toward zero on the scale of the (as supplied) \code{data};
#'   it is only statistically meaningful as a sparsity-inducing penalty on
#'   centered/scaled data (so that zero corresponds to "no departure from the
#'   grand mean"), and, because the threshold is a flat amount rather than one
#'   rescaled by each profile's variance/sample size, it is a computationally
#'   convenient approximation to (not an exact coordinate-wise solution of)
#'   the corresponding L1-penalized weighted log-likelihood except when
#'   profile covariances are close to the identity -- as they will be on
#'   standardized data with roughly independent variables. Standardize
#'   \code{data} first if you intend to use \code{lambda > 0}; a
#'   \code{warning()} is raised if \code{data} does not look centered/scaled.
#'   Default \code{0} (no shrinkage).
#' @param mcmc_iter Number of iterations per MCMC chain (see \code{n_chains}).
#'   The first half of each chain is discarded as burn-in before computing
#'   posterior summaries and convergence diagnostics. Ignored when \code{engine = "EM"}.
#' @param prior_laplace Positive numeric, the Laplace (Bayesian Lasso)
#'   shrinkage/rate hyperparameter for the profile means under the MCMC
#'   engine (denoted \eqn{\lambda} in Park & Casella, 2008; \emph{not} a
#'   dispersion "scale" in the usual sense -- \strong{larger} values induce
#'   \strong{more} shrinkage of the means toward zero, analogous to a bigger
#'   \code{lambda} in the EM engine). Ignored when \code{engine = "EM"}.
#' @param robust Logical. If \code{TRUE} (default), robust (outlier
#'   down-weighted) estimation is used regardless of \code{engine}: Huber
#'   weights, computed from per-observation squared Mahalanobis distances,
#'   down-weight outliers when accumulating the sufficient statistics used
#'   for each profile's mean/covariance -- in the EM engine's M-step (see
#'   \code{\link{robust_m_step}}) and, in the same spirit, in the MCMC
#'   engine's Gibbs updates (see "Robust estimation" below). If \code{FALSE},
#'   this down-weighting is skipped in both engines, which reduce to
#'   classical (non-robust) maximum-likelihood / Bayesian estimation for a
#'   Gaussian mixture -- useful if you want to fit a standard LPA/GMM with
#'   the same interface, missing-data handling, and model parameterizations
#'   as the rest of this package.
#' @param alpha Significance level for the Huber down-weighting threshold
#'   (observations with squared Mahalanobis distance beyond the
#'   \code{1 - alpha} chi-squared quantile are down-weighted). Smaller values
#'   down-weight fewer, more extreme points; larger values down-weight more
#'   aggressively. Default \code{0.05}. Ignored when \code{robust = FALSE}.
#' @param n_chains Number of independent MCMC chains to run (default
#'   \code{4}). Running multiple chains from independent starting values is
#'   what makes the Gelman-Rubin \eqn{\hat{R}} diagnostic possible (it is a
#'   between- vs. within-chain comparison and is undefined for a single
#'   chain). Posterior summaries (\code{means}, \code{covariances},
#'   \code{proportions}) pool the post-burn-in draws of all chains together.
#'   Ignored when \code{engine = "EM"}.
#' @param cores Integer, number of CPU cores to use for parallel estimation
#'   \emph{within this single} \code{robust_lpa()} call (default \code{1},
#'   sequential). For \code{engine = "EM"}, parallelizes across the
#'   \code{n_starts} random restarts. For \code{engine = "MCMC"},
#'   parallelizes across the \code{n_chains} chains. Uses
#'   \code{parallel::mclapply()} (forking) on Unix-alikes, and a
#'   \code{parallel::makeCluster()} PSOCK cluster -- with independent
#'   per-worker RNG streams via \code{parallel::clusterSetRNGStream()} -- on
#'   Windows; falls back to sequential execution with a \code{warning()} if
#'   the \pkg{parallel} package is unavailable. This is a different axis of
#'   parallelism from, and should not usually be combined with (nested
#'   parallelism can oversubscribe your CPU), the \code{cores} argument of
#'   \code{\link{estimate_profiles_robust}}, which instead parallelizes
#'   across \code{G}/\code{model} combinations.
#' @return A list with S3 class \code{"robust_lpa"} (see
#'   \code{\link{print.robust_lpa}} and \code{\link{summary.robust_lpa}} for
#'   concise and detailed views of the fit) containing:
#'   \describe{
#'     \item{engine}{The estimation engine used.}
#'     \item{means}{A list of length \code{G} with the estimated profile means.}
#'     \item{covariances}{A list of length \code{G} with the estimated profile covariance matrices.}
#'     \item{proportions}{Numeric vector of length \code{G} with the estimated mixing proportions.}
#'     \item{probabilities}{An \code{n x G} matrix of posterior profile-membership probabilities.}
#'     \item{fit}{A one-row data.frame with \code{Model}, \code{Profiles}, \code{LogLik}, \code{Parameters}, \code{AIC}, \code{BIC}, \code{SABIC}, \code{Entropy}, \code{Min_Size}, and \code{Max_Size}.}
#'     \item{assignments}{Integer vector of length \code{n} with the most likely profile for each observation.}
#'     \item{mcmc_draws}{(MCMC engine only) a list with \code{chains} (the raw per-chain draws, as consumed by \code{\link{plot_mcmc_chains}}), \code{n_chains}, \code{mcmc_iter}, and \code{burnin}.}
#'     \item{mcmc_diagnostics}{(MCMC engine only) a data.frame with one row per scalar parameter (\code{Parameter}, \code{Rhat}, \code{ESS}); see "MCMC convergence diagnostics" below. \code{NULL} if the \pkg{coda} package is unavailable or there are too few post-burn-in iterations.}
#'     \item{data}{The numeric matrix actually fit (\code{data} coerced via \code{as.matrix()}).}
#'     \item{call_args}{A named list of every argument controlling this fit (\code{G}, \code{model}, \code{engine}, ...), for internal reuse -- e.g. \code{\link{bch_robust}}'s \code{correction = "bootstrap"} refits this exact specification on resampled data via \code{do.call(robust_lpa, modifyList(call_args, list(data = new_data)))}.}
#'   }
#' @section Robust estimation:
#' Both engines share the same outlier-down-weighting idea, adapted to how
#' each one accumulates information:
#' \itemize{
#'   \item \strong{EM}: at every M-step, each profile's mean/covariance is
#'     recomputed from Huber-down-weighted, posterior-probability-weighted
#'     observations (\code{\link{robust_m_step}}).
#'   \item \strong{MCMC}: at every Gibbs sweep, after observations are
#'     allocated to profiles, a Huber weight is computed for each observation
#'     from its squared Mahalanobis distance to its \emph{currently assigned}
#'     profile's previous-sweep mean/covariance (the same \code{alpha}
#'     chi-squared cutoff used by the EM engine). These weights down-weight
#'     the sufficient statistics that drive that sweep's mean/covariance
#'     updates, so an outlying observation contributes a smaller effective
#'     sample size to its profile's posterior. As in the EM engine, the
#'     mixing-proportion update is unaffected (it uses the raw allocation
#'     counts) and the allocation step itself is not down-weighted. Setting
#'     \code{robust = FALSE} recovers a standard (non-robust) Gibbs sampler
#'     for the same Bayesian-Lasso Gaussian mixture. A row of \code{data}
#'     with no observed variables at all is still (re)allocated to a profile
#'     every sweep, drawn from the current mixing proportions (its posterior
#'     given zero data), rather than being left permanently stuck in a
#'     single profile.
#' }
#' @section MCMC convergence diagnostics:
#' When \code{engine = "MCMC"}, \code{$mcmc_diagnostics} reports, for every
#' scalar parameter (\code{"mu[g,j]"}, \code{"sigma[g,j]"}, \code{"pi[g]"}):
#' \itemize{
#'   \item \strong{Rhat}: the classic Gelman-Rubin potential scale reduction
#'     statistic (Gelman & Rubin, 1992), comparing between- and within-chain
#'     variance on the post-burn-in draws. Values noticeably above 1.1 are
#'     the classic rule-of-thumb warning sign of non-convergence, and trigger
#'     a \code{warning()}. \code{NA} when \code{n_chains = 1} (undefined for
#'     a single chain).
#'   \item \strong{ESS}: the classic (autocorrelation/spectral-density-based)
#'     effective sample size, i.e. how many independent draws the
#'     autocorrelated post-burn-in draws (pooled across chains) are worth.
#' }
#' Computing these requires the \pkg{coda} package; \code{$mcmc_diagnostics}
#' is \code{NULL} (with a \code{warning()}) if it is not installed.
#' @seealso \code{\link{estimate_profiles_robust}} to fit and compare many
#'   \code{G} / \code{model} combinations at once, \code{\link{blrt_robust}}
#'   for a bootstrapped likelihood ratio test to choose the number of
#'   profiles, and \code{\link{plot_mcmc_chains}} to inspect MCMC chains.
#' @references
#'   Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#'   using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'   \doi{10.1214/ss/1177011136}
#'
#'   Park, T., & Casella, G. (2008). The Bayesian Lasso. \emph{Journal of the
#'   American Statistical Association}, 103(482), 681-686.
#'   \doi{10.1198/016214508000000337}
#' @examples
#' # Fast demonstration on the bundled `neuro_data` dataset (standardized,
#' # as recommended -- see `lambda` and `prior_laplace` above).
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "Attention", "Executive_Functions",
#'                                      "RT_Stroop", "RT_TMT")]))
#' # Huber-weighted robust estimation can occasionally emit a log-likelihood
#' # decrease warning as an expected side effect of down-weighting outliers
#' # mid-fit (see the "Robust estimation" section above); wrapped in
#' # suppressWarnings() below for a clean example, not because it signals a
#' # problem with the fit.
#' fit <- suppressWarnings(robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 30))
#' fit           # print.robust_lpa(): concise overview of the fit
#' summary(fit)  # summary.robust_lpa(): profile means/sizes and full fit indices
#'
#' # `neuro_data` injects extra, variable-magnitude outliers into RT_Stroop
#' # and RT_TMT for a subset of the Pathological group; compare robust
#' # (default) vs. classical (robust = FALSE) estimation of the profile means.
#' fit_robust <- suppressWarnings(robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 30))
#' fit_classical <- robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 30, robust = FALSE)
#' sapply(fit_robust$means, `[`, "RT_Stroop")
#' sapply(fit_classical$means, `[`, "RT_Stroop")
#'
#' \donttest{
#' # MCMC engine: 4 chains (default), with a small mcmc_iter for speed, and
#' # the resulting Gelman-Rubin / ESS convergence diagnostics.
#' fit_mcmc <- suppressWarnings(robust_lpa(x, G = 2, model = 6, engine = "MCMC",
#'                                          mcmc_iter = 200, n_chains = 4))
#' summary(fit_mcmc)  # includes $mcmc_diagnostics (Rhat / ESS) in the printout
#' }
#' @export
robust_lpa <- function(data, G, model = 6, engine = "EM", max_iter = 100, tol = 1e-6,
                       n_starts = 5, lambda = 0, mcmc_iter = 2000, prior_laplace = 0.1,
                       robust = TRUE, alpha = 0.05, n_chains = 4, cores = 1) {
  
  # ---- input validation -----------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a matrix or data.frame.")
  }
  X <- as.matrix(data)
  storage.mode(X) <- "double"
  if (!is.numeric(X)) stop("`data` must contain only numeric (or coercible-to-numeric) columns.")
  
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(n) || n < 2) stop("`data` must have at least 2 rows.")
  if (is.null(p) || p < 1) stop("`data` must have at least 1 column.")
  if (!is.numeric(G) || length(G) != 1 || G < 1 || G != round(G)) {
    stop("`G` must be a single positive integer.")
  }
  if (!(length(model) == 1 && model %in% 1:6)) {
    stop("`model` must be a single integer between 1 and 6.")
  }
  if (!(length(engine) == 1 && engine %in% c("EM", "MCMC"))) {
    stop("`engine` must be either 'EM' or 'MCMC'.")
  }
  if (!is.numeric(max_iter) || max_iter < 1) stop("`max_iter` must be a positive integer.")
  if (!is.numeric(tol) || tol <= 0) stop("`tol` must be a positive number.")
  if (!is.numeric(n_starts) || n_starts < 1) stop("`n_starts` must be at least 1.")
  if (!is.numeric(lambda) || lambda < 0) stop("`lambda` must be non-negative.")
  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) stop("`robust` must be a single TRUE/FALSE value.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("`alpha` must be a single number strictly between 0 and 1.")
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive integer.")
  }
  if (G > n) stop("`G` cannot exceed the number of observations.")
  if (engine == "MCMC") {
    if (!is.numeric(mcmc_iter) || mcmc_iter < 2) stop("`mcmc_iter` must be at least 2.")
    if (!is.numeric(prior_laplace) || prior_laplace <= 0) stop("`prior_laplace` must be positive.")
    if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains < 1 || n_chains != round(n_chains)) {
      stop("`n_chains` must be a single positive integer.")
    }
  }
  if (lambda > 0) {
    # `lambda` soft-thresholds the raw-scale mean (see @param lambda); flag
    # data that clearly is not centered/scaled so the shrinkage is
    # interpretable, without being so strict that legitimately-scaled data
    # (SD not exactly 1, mean not exactly 0) triggers spurious warnings.
    col_means <- colMeans(X, na.rm = TRUE)
    col_sds <- apply(X, 2, stats::sd, na.rm = TRUE)
    looks_unscaled <- any(col_sds < 0.5 | col_sds > 2, na.rm = TRUE) ||
      any(abs(col_means) > 0.5 * col_sds, na.rm = TRUE)
    if (isTRUE(looks_unscaled)) {
      warning(
        "`lambda` > 0 (LASSO soft-thresholding) shrinks each mean component by a flat amount on ",
        "the raw scale of `data`, and is only statistically meaningful as a sparsity penalty on ",
        "centered/scaled data. `data` does not look centered/scaled (column means/SDs are far ",
        "from 0/1); consider fitting on `scale(data)` instead, or interpret the shrinkage with caution."
      )
    }
  }

  # ---- parallel dispatch (EM restarts / MCMC chains) --------------------
  # `.run_parallel()` (defined once in R/parallel_utils.R and shared with
  # `blrt_robust()`, which parallelizes its bootstrap replicates the same
  # way) falls back to sequential `lapply()` if `cores == 1`, the `parallel`
  # package is unavailable, or there is only one unit of work to run anyway.

  # Snapshot of every argument that controls how this model was fit, stored
  # on the returned object (see `call_args`/`data` in `@return`) so that
  # other functions can refit the identical specification on different data
  # via `do.call(robust_lpa, modifyList(fit$call_args, list(data = new_data)))`
  # -- used internally by `bch_robust(correction = "bootstrap")` to refit the
  # model on bootstrap resamples without asking the caller to repeat every
  # argument.
  call_args <- list(
    G = G, model = model, engine = engine, max_iter = max_iter, tol = tol,
    n_starts = n_starts, lambda = lambda, mcmc_iter = mcmc_iter,
    prior_laplace = prior_laplace, robust = robust, alpha = alpha,
    n_chains = n_chains, cores = cores
  )

  has_na <- any(is.na(X))
  var_names <- colnames(X)

  if (engine == "MCMC") {
    # Run `n_chains` independent chains from independent (identical, since
    # the sampler always initializes at mu = 0 / sigma = I) starting states;
    # independence comes from each chain drawing its own random sequence
    # from R's RNG stream, not from distinct starting values. This is what
    # makes the between- vs. within-chain Gelman-Rubin comparison in
    # .compute_mcmc_diagnostics() meaningful (see the "MCMC convergence
    # diagnostics" section of this function's documentation). With
    # `cores > 1`, the chains run in parallel (see `?robust_lpa`'s `cores`).
    run_one_chain <- function(chain_id) {
      robust_mcmc_cpp(
        X = X, G = G, model = model, mcmc_iter = mcmc_iter,
        prior_laplace = prior_laplace, robust = robust, alpha = alpha
      )
    }
    chains_raw <- .run_parallel(
      cores, n_chains, run_one_chain,
      export_vars = c("X", "G", "model", "mcmc_iter", "prior_laplace", "robust", "alpha"),
      export_env = environment()
    )

    # Align every chain's (and every iteration's) profile labels to a common
    # ordering before anything downstream touches them -- raw multi-chain
    # draws are NOT directly comparable/poolable otherwise (see
    # .relabel_mcmc_chains()'s documentation: "label switching"). Must happen
    # before burn-in is discarded and before diagnostics/pooling below, since
    # both would otherwise silently blend different real-world profiles
    # together (symptom: pooled profile means that look nearly identical to
    # each other, and a classification entropy near 0).
    chains_raw <- .relabel_mcmc_chains(chains_raw, G)

    burnin <- floor(mcmc_iter / 2)
    valid_iters <- (burnin + 1):mcmc_iter
    n_valid <- length(valid_iters)

    # Pool the post-burn-in draws of every chain together for the posterior
    # point estimates (means, covariances, mixing proportions). This is only
    # appropriate once the chains have (approximately) converged to the same
    # target distribution -- see mcmc_diagnostics / the Rhat warning below.
    pi_draws <- do.call(rbind, lapply(chains_raw, function(ch) ch$pi_chain[valid_iters, , drop = FALSE]))
    pi_g <- colMeans(pi_draws)

    mu <- list()
    sigma <- list()
    n_pooled <- n_valid * n_chains

    for (g in 1:G) {
      mu_g <- matrix(0, nrow = n_pooled, ncol = p)
      sigma_g <- matrix(0, nrow = p, ncol = p)
      row_i <- 0

      for (chain_id in seq_len(n_chains)) {
        ch <- chains_raw[[chain_id]]
        for (iter in valid_iters) {
          row_i <- row_i + 1
          mu_g[row_i, ] <- as.numeric(ch$mu_chain[[iter]][[g]])
          sigma_g <- sigma_g + ch$sigma_chain[[iter]][[g]]
        }
      }

      mu[[g]] <- colMeans(mu_g)
      sigma[[g]] <- sigma_g / n_pooled
    }

    # Reattach the original variable names, lost when means/covariances pass
    # through the (unnamed, arma-based) C++/MCMC engine, so that downstream
    # use (printing, plot_robust_lpa(), indexing model$means[[g]]["var"], etc.)
    # shows meaningful labels instead of anonymous positions.
    if (!is.null(var_names)) {
      for (g in 1:G) {
        names(mu[[g]]) <- var_names
        dimnames(sigma[[g]]) <- list(var_names, var_names)
      }
    }

    density_matrix <- matrix(0, nrow = n, ncol = G)

    if (!has_na) {
      for (g in 1:G) {
        density_matrix[, g] <- dmvnorm_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
      }
    } else {
      for (g in 1:G) {
        density_matrix[, g] <- dmvnorm_fiml_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
      }
    }

    row_sums <- rowSums(density_matrix)
    row_sums[row_sums < 1e-300] <- 1e-300

    log_lik <- sum(log(row_sums))
    z <- density_matrix / row_sums

    fi <- .compute_fit_indices(mu, log_lik, z, n, p, G, model)
    .warn_if_degenerate_profile(fi$fit_indices, engine)

    diagnostics <- .compute_mcmc_diagnostics(chains_raw, mcmc_iter = mcmc_iter, burnin = burnin)
    if (!is.null(diagnostics) && isTRUE(any(diagnostics$Rhat > 1.1, na.rm = TRUE))) {
      warning(
        "Some MCMC parameters have a Gelman-Rubin R-hat > 1.1, suggesting the chains have not ",
        "fully converged. Consider increasing `mcmc_iter`. See `model$mcmc_diagnostics`."
      )
    }

    out <- list(
      engine = engine, means = mu, covariances = sigma,
      proportions = pi_g, probabilities = z, fit = fi$fit_indices,
      assignments = fi$assignments,
      mcmc_draws = list(chains = chains_raw, n_chains = n_chains, mcmc_iter = mcmc_iter, burnin = burnin),
      mcmc_diagnostics = diagnostics,
      data = X, call_args = call_args
    )
    class(out) <- "robust_lpa"
    return(out)

  } else if (engine == "EM") {
    # One random-restart EM run. Returns everything needed to (a) pick the
    # best start by log-likelihood and (b) build the final fitted model,
    # without depending on any state outside its own arguments/closure --
    # required so it can be dispatched to parallel workers via
    # `.run_parallel()` (see `cores`) exactly as easily as run sequentially.
    # `warnings_out` accumulates log-likelihood-decrease messages (see
    # below) instead of calling `warning()` directly, so that (i) they are
    # never silently dropped when this runs on a forked/PSOCK worker (where
    # a bare `warning()` call may not propagate back to the caller), and
    # (ii) every start's warnings are still surfaced, not just the winning
    # start's, exactly matching the original sequential behavior.
    run_one_start <- function(start) {
      pi_g <- rep(1 / G, G)
      mu <- list()
      sigma <- list()
      warnings_out <- character(0)

      z_init <- matrix(runif(n * G), nrow = n, ncol = G)
      z_init <- z_init / rowSums(z_init)

      for (g in 1:G) {
        init_results <- robust_m_step(data = X, z = z_init[, g], alpha = alpha, lambda = lambda, robust = robust)
        mu[[g]] <- init_results$mean
        sigma[[g]] <- .force_pd(init_results$covariance)
      }

      log_lik <- -Inf

      # Track the best (highest log-lik) iterate seen along this start's
      # trajectory, and return THAT instead of wherever the loop happens to
      # end up. A correctly-specified EM M-step should never decrease the
      # observed-data log-likelihood by more than numerical noise, but two
      # sources can legitimately cause small, and occasionally sustained,
      # decreases: models 4 and 5 use an approximate constrained M-step (see
      # .apply_covariance_model()), and .force_pd()'s positive-definiteness
      # safety net can perturb any model's covariance away from the
      # unconstrained M-step optimum when it is near-singular -- which can
      # itself happen repeatedly in a row while a profile's covariance is
      # collapsing toward a spurious near-singular solution (see
      # .force_pd()'s documentation). Keeping the best iterate is a
      # strictly-safe improvement regardless of the cause. `n_decreases` /
      # `max_single_decrease` accumulate across the whole trajectory so a
      # SINGLE consolidated warning can be issued after the loop, instead of
      # one warning per offending iteration (which floods the console when a
      # run decreases for many iterations in a row).
      best_log_lik <- -Inf
      best_mu <- mu; best_sigma <- sigma; best_pi_g <- pi_g; best_z <- NULL
      n_decreases <- 0L
      max_single_decrease <- 0
      consecutive_decreases <- 0L
      divergence_iter <- NA_integer_
      final_iter <- 0L

      for (iter in 1:max_iter) {
        final_iter <- iter
        density_matrix <- matrix(0, nrow = n, ncol = G)

        if (!has_na) {
          for (g in 1:G) {
            density_matrix[, g] <- dmvnorm_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
          }
        } else {
          for (g in 1:G) {
            density_matrix[, g] <- dmvnorm_fiml_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
          }
        }

        row_sums <- rowSums(density_matrix)
        row_sums[row_sums < 1e-300] <- 1e-300

        new_log_lik <- sum(log(row_sums))
        z <- density_matrix / row_sums

        # The threshold below (1e-4, versus a naive 1e-6) absorbs ordinary
        # numerical noise so bookkeeping is reserved for decreases actually
        # worth knowing about; log_lik starts at -Inf, so this never counts
        # the first iteration.
        if (new_log_lik < log_lik - 1e-4) {
          n_decreases <- n_decreases + 1L
          max_single_decrease <- max(max_single_decrease, log_lik - new_log_lik)
          consecutive_decreases <- consecutive_decreases + 1L
        } else {
          consecutive_decreases <- 0L
        }

        if (new_log_lik > best_log_lik) {
          best_log_lik <- new_log_lik
          best_mu <- mu; best_sigma <- sigma; best_pi_g <- pi_g; best_z <- z
        }

        converged <- abs(new_log_lik - log_lik) < tol
        log_lik <- new_log_lik
        if (converged) break

        # Sustained decreases (as opposed to a single small blip) mean the
        # loop has stopped making progress -- typically a profile's
        # covariance oscillating around a near-singular collapse (see
        # .force_pd()). Stop wasting iterations on it; the best iterate
        # tracked above is still returned.
        if (consecutive_decreases >= 5L) {
          divergence_iter <- iter
          break
        }

        raw_sigmas <- list()
        for (g in 1:G) {
          m_step_results <- robust_m_step(data = X, z = z[, g], alpha = alpha, lambda = lambda, robust = robust)
          mu[[g]] <- m_step_results$mean
          raw_sigmas[[g]] <- m_step_results$covariance
          pi_g[g] <- mean(z[, g])
        }

        sigma <- .apply_covariance_model(raw_sigmas, pi_g, model, G, p)
      }

      if (n_decreases > 0L) {
        cause <- if (model %in% c(4, 5)) {
          "the approximate constrained M-step used for models 4 and 5"
        } else if (model %in% c(3, 6)) {
          "the positive-definiteness safety net (.force_pd()) perturbing a near-singular covariance"
        } else if (robust) {
          "the robust (Huber-weighted) M-step, which is an approximate one-step update and does not carry the same exact monotonic-increase guarantee as classical (non-robust) EM"
        } else {
          "ordinary numerical noise near convergence"
        }
        divergence_note <- if (!is.na(divergence_iter)) {
          sprintf(
            " Stopped early at iteration %d after %d consecutive decreases (likely a profile's covariance collapsing toward a spurious near-singular solution); the best iterate found (log-lik = %.6f) was kept instead.",
            divergence_iter, 5L, best_log_lik
          )
        } else {
          sprintf(" The best iterate found (log-lik = %.6f) was kept instead of the final one.", best_log_lik)
        }
        warnings_out <- c(warnings_out, sprintf(
          "Start %d: log-likelihood decreased in %d of %d iteration(s) (largest single decrease %.6f) for model %d. This can happen with %s; consider comparing against a different `model`, increasing `n_starts`, or checking for an implausibly small profile size if this persists.%s",
          start, n_decreases, final_iter, max_single_decrease, model, cause, divergence_note
        ))
      }

      list(log_lik = best_log_lik, mu = best_mu, sigma = best_sigma, pi_g = best_pi_g, z = best_z, warnings = warnings_out)
    }

    start_results <- .run_parallel(
      cores, n_starts, run_one_start,
      export_vars = c("X", "n", "p", "G", "model", "max_iter", "tol", "lambda", "robust", "alpha", "has_na"),
      export_env = environment()
    )

    for (res in start_results) {
      for (w in res$warnings) warning(w, call. = FALSE)
    }

    log_liks <- vapply(start_results, function(r) r$log_lik, numeric(1))
    if (!any(is.finite(log_liks))) {
      # Every start ended in a degenerate fit (all-zero densities every
      # iteration). The original implementation silently returned NULL here
      # (its `best_log_lik > -Inf` selection check can never fire in this
      # case); an explicit, actionable error is preferable to a silent NULL
      # that would otherwise surface later as a confusing "$ operator
      # invalid" error wherever the caller first uses the result.
      stop(
        "All ", n_starts, " EM start(s) produced a non-finite log-likelihood (model ", model,
        ", G = ", G, "). This usually indicates a numerically degenerate fit (e.g. too few ",
        "observations for the requested number of profiles/parameters); try a different `model`, ",
        "fewer profiles, or more `n_starts`."
      )
    }
    best <- start_results[[which.max(log_liks)]]
    mu <- best$mu; sigma <- best$sigma; pi_g <- best$pi_g; z <- best$z; best_log_lik <- best$log_lik

    # Reattach original variable names (see the analogous step in the
    # MCMC branch above for why this is needed).
    if (!is.null(var_names)) {
      for (g in 1:G) {
        names(mu[[g]]) <- var_names
        dimnames(sigma[[g]]) <- list(var_names, var_names)
      }
    }
    fi <- .compute_fit_indices(mu, best_log_lik, z, n, p, G, model)
    .warn_if_degenerate_profile(fi$fit_indices, engine)
    best_model <- list(engine = engine, means = mu, covariances = sigma, proportions = pi_g,
                       probabilities = z, fit = fi$fit_indices, assignments = fi$assignments,
                       data = X, call_args = call_args)
    class(best_model) <- "robust_lpa"
    return(best_model)
  } else {
    stop("Unsupported engine. Please choose either 'EM' or 'MCMC'.")
  }
}

#' Force a Matrix to be Symmetric Positive Definite
#'
#' Symmetrizes a matrix and, if its smallest eigenvalue is below a floor,
#' adds a ridge to the diagonal so that it becomes safely positive definite.
#' Used as a numerical safety net after every covariance update for models
#' 3-6 in \code{\link{robust_lpa}} (models 1-2 use a separate, simpler
#' per-variable floor -- see \code{\link{.apply_covariance_model}}).
#'
#' The floor is set \emph{relative to the matrix's own average eigenvalue}
#' (\code{min_ratio} times \code{mean(eigenvalues)}), not a fixed absolute
#' constant. An earlier version used a fixed constant (effectively ~1e-5),
#' which is a negligible floor relative to the variance scale of typical
#' (e.g. standardized, variance ~= 1) data -- too weak to stop a profile's
#' covariance from collapsing toward a near-singular, near-zero-variance
#' matrix that fits a tiny, near-coincident handful of points almost
#' perfectly. That is a classic spurious/degenerate local maximum of the
#' unconstrained Gaussian mixture likelihood (see e.g. McLachlan & Peel,
#' 2000, \emph{Finite Mixture Models}, Sec. 3.10): the resulting "profile"
#' can have both an implausibly high log-likelihood contribution and a
#' vanishingly small size, and is a common way multivariate EM for
#' unconstrained-covariance mixtures goes wrong in any implementation, not
#' specific to this package. A relative floor scales automatically with
#' whatever units \code{data} is on and meaningfully resists this collapse,
#' at the cost of (very slightly) regularizing genuinely tight-but-real
#' clusters toward less extreme covariances. This does not eliminate the
#' possibility of a spurious solution -- no finite floor can, since it is a
#' structural property of unconstrained-covariance finite mixture ML -- but
#' substantially raises how tight a cluster has to be before it can trigger
#' one. If you see an implausibly small profile size for model 3-6, compare
#' against a smaller \code{G}, more \code{n_starts}, or a more constrained
#' \code{model} (1-2, which do not go through this floor).
#'
#' @param mat A square numeric matrix.
#' @param min_ratio Numeric, the minimum eigenvalue as a fraction of
#'   \code{mat}'s own average eigenvalue. Default \code{1e-3}.
#' @return A symmetric positive-definite matrix of the same dimension.
#' @keywords internal
#' @noRd
.force_pd <- function(mat, min_ratio = 1e-3) {
  mat <- (mat + t(mat)) / 2
  ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  floor_ev <- max(mean(ev) * min_ratio, 1e-8)
  min_ev <- min(ev)
  if (min_ev < floor_ev) {
    mat <- mat + diag(floor_ev - min_ev, nrow(mat))
  }
  mat
}

#' Apply a Variance-Covariance Model to a Set of Raw Profile Covariances
#'
#' Given the (unconstrained) per-profile covariance matrices produced by the
#' M-step (\code{\link{robust_m_step}}), builds the constrained
#' variance-covariance matrices implied by \code{model} (see
#' \code{\link{robust_lpa}} for the definition of models 1 to 6).
#'
#' For models 4 and 5 (mixed variance/covariance constraints), the shared
#' component is built as a shared \emph{correlation} matrix and then rescaled
#' by each profile's own or the pooled variances, mirroring the construction
#' used by the MCMC engine (\code{robust_mcmc_cpp()}). This is a
#' well-defined, always-positive-semi-definite-before-regularization
#' construction; a naive "swap the diagonal of the pooled covariance matrix"
#' approach does not preserve positive-definiteness and does not correspond
#' to the maximizer of the constrained expected complete-data likelihood.
#'
#' @param raw_sigmas A list of length \code{G} with unconstrained covariance matrices.
#' @param pi_g Numeric vector of length \code{G} with the current mixing proportions.
#' @param model An integer between 1 and 6.
#' @param G Integer, the number of profiles.
#' @param p Integer, the number of variables.
#' @return A list of length \code{G} with the constrained covariance matrices.
#' @keywords internal
#' @noRd
.apply_covariance_model <- function(raw_sigmas, pi_g, model, G, p) {
  sigma <- vector("list", G)
  
  if (model == 1) {
    # Equal variances across profiles, covariances fixed to 0 (diagonal, shared
    # across profiles; not a single isotropic/spherical variance across variables).
    pooled_diag <- numeric(p)
    for (g in 1:G) pooled_diag <- pooled_diag + pi_g[g] * diag(raw_sigmas[[g]])
    pooled_diag <- pmax(pooled_diag, 1e-8)
    for (g in 1:G) sigma[[g]] <- diag(pooled_diag, p)
    
  } else if (model == 2) {
    # Varying variances across profiles, covariances fixed to 0.
    for (g in 1:G) sigma[[g]] <- diag(pmax(diag(raw_sigmas[[g]]), 1e-8), p)
    
  } else if (model == 3) {
    # Equal variances and equal covariances across profiles (one shared,
    # full covariance matrix).
    pooled_sigma <- matrix(0, nrow = p, ncol = p)
    for (g in 1:G) pooled_sigma <- pooled_sigma + pi_g[g] * raw_sigmas[[g]]
    pooled_sigma <- .force_pd(pooled_sigma)
    for (g in 1:G) sigma[[g]] <- pooled_sigma
    
  } else if (model == 4) {
    # Varying variances, equal covariance *structure*: the correlation
    # matrix is shared across profiles, but each profile keeps its own
    # variances. Build the shared correlation matrix from the pi-weighted
    # pooled covariance, then rescale it by each profile's own variances.
    pooled_sigma <- matrix(0, nrow = p, ncol = p)
    for (g in 1:G) pooled_sigma <- pooled_sigma + pi_g[g] * raw_sigmas[[g]]
    pooled_sigma <- .force_pd(pooled_sigma)
    pooled_sd <- sqrt(pmax(diag(pooled_sigma), 1e-8))
    R_pool <- pooled_sigma / outer(pooled_sd, pooled_sd)
    
    for (g in 1:G) {
      sd_g <- sqrt(pmax(diag(raw_sigmas[[g]]), 1e-8))
      sigma[[g]] <- .force_pd(R_pool * outer(sd_g, sd_g))
    }
    
  } else if (model == 5) {
    # Equal variances across profiles, varying covariances: each profile
    # keeps its own correlation matrix, but all profiles share the same
    # (pooled) variances.
    pooled_diag <- numeric(p)
    for (g in 1:G) pooled_diag <- pooled_diag + pi_g[g] * diag(raw_sigmas[[g]])
    pooled_sd <- sqrt(pmax(pooled_diag, 1e-8))
    
    for (g in 1:G) {
      sd_g <- sqrt(pmax(diag(raw_sigmas[[g]]), 1e-8))
      R_g <- raw_sigmas[[g]] / outer(sd_g, sd_g)
      sigma[[g]] <- .force_pd(R_g * outer(pooled_sd, pooled_sd))
    }
    
  } else if (model == 6) {
    # Fully unconstrained: each profile has its own variances and covariances.
    for (g in 1:G) sigma[[g]] <- .force_pd(raw_sigmas[[g]])
    
  } else {
    stop("`model` must be an integer between 1 and 6.")
  }
  
  sigma
}

#' Compute Fit Indices and Hard Assignments for a Fitted Mixture
#'
#' Shared helper used by both the EM and MCMC branches of
#' \code{\link{robust_lpa}} to compute AIC/BIC/SABIC, classification entropy,
#' profile sizes, and modal (hard) assignments from a fitted model's means
#' and posterior probabilities. Centralizing this logic keeps the two
#' engines' fit indices consistent by construction.
#'
#' @param mu A list of length \code{G} with the estimated profile means.
#' @param log_lik The observed-data log-likelihood of the fitted model.
#' @param z An \code{n x G} matrix of posterior profile-membership probabilities.
#' @param n Integer, the number of observations.
#' @param p Integer, the number of variables.
#' @param G Integer, the number of profiles.
#' @param model An integer between 1 and 6.
#' @return A list with \code{fit_indices} (a one-row data.frame) and
#'   \code{assignments} (an integer vector of length \code{n}).
#' @keywords internal
#' @noRd
.compute_fit_indices <- function(mu, log_lik, z, n, p, G, model) {
  if (G == 1) {
    K_base <- switch(as.character(model),
                     "1" = p + p, "2" = p + p, "3" = p + (p * (p + 1) / 2),
                     "4" = p + (p * (p + 1) / 2), "5" = p + (p * (p + 1) / 2), "6" = p + (p * (p + 1) / 2))
  } else {
    K_base <- switch(as.character(model),
                     "1" = (G * p) + p + (G - 1), "2" = (G * p) + (G * p) + (G - 1),
                     "3" = (G * p) + (p * (p + 1) / 2) + (G - 1), "4" = (G * p) + (G * p) + (p * (p - 1) / 2) + (G - 1),
                     "5" = (G * p) + p + (G * p * (p - 1) / 2) + (G - 1), "6" = (G * p) + (G * (p * (p + 1) / 2)) + (G - 1))
  }
  
  zeroed_means <- sum(sapply(mu, function(m) sum(abs(m) < 1e-5)))
  K_eff <- max(1, K_base - zeroed_means)
  
  AIC_val <- -2 * log_lik + 2 * K_eff
  BIC_val <- -2 * log_lik + K_eff * log(n)
  SABIC_val <- -2 * log_lik + K_eff * log((n + 2) / 24)
  
  z_safe <- z
  z_safe[z_safe < 1e-15] <- 1e-15
  entropy <- if (G == 1) 1 else 1 - (sum(-z * log(z_safe)) / (n * log(G)))
  
  assignments <- max.col(z)
  sizes <- table(factor(assignments, levels = 1:G))
  
  list(
    fit_indices = data.frame(
      Model = model, Profiles = G, LogLik = log_lik, Parameters = K_eff,
      AIC = AIC_val, BIC = BIC_val, SABIC = SABIC_val, Entropy = entropy,
      Min_Size = min(sizes) / n, Max_Size = max(sizes) / n
    ),
    assignments = assignments
  )
}

#' Warn When the Smallest Fitted Profile Looks Implausibly Small
#'
#' A transparency check run on every \code{\link{robust_lpa}} fit: an
#' implausibly small smallest-profile size (well below what a genuine small
#' subgroup would usually look like) is the standard symptom of a
#' spurious/degenerate solution -- most often a profile's covariance having
#' collapsed toward a tiny, near-coincident handful of points (see
#' \code{\link{.force_pd}}'s documentation), which is structurally more
#' likely for models 3-6 (which allow off-diagonal/full covariance) than for
#' models 1-2 (diagonal only). This does not attempt to distinguish a
#' genuine small group from a spurious one -- that requires judgment (compare
#' across `model`, more `n_starts`, fewer profiles) -- it only flags the
#' pattern so it isn't missed, especially by code that automatically selects
#' a "best" fit by AIC/BIC (e.g. \code{\link{plot_robust_lpa}} on an
#' \code{\link{estimate_profiles_robust}} result), where a spurious solution
#' can otherwise look deceptively attractive (an artificially high
#' log-likelihood from the near-singular covariance).
#'
#' @param fit_indices The one-row fit-indices data.frame from \code{\link{.compute_fit_indices}}.
#' @param engine String, \code{"EM"} or \code{"MCMC"} (for the warning text only).
#' @return \code{NULL}, invisibly; called for its \code{warning()} side effect.
#' @keywords internal
#' @noRd
.warn_if_degenerate_profile <- function(fit_indices, engine) {
  if (isTRUE(fit_indices$Min_Size < 0.02)) {
    warning(sprintf(
      "The smallest fitted profile holds only %.1f%% of the observations (model %d, %d profile(s), %s engine). This can be a genuine small group, but is also the classic symptom of a spurious/degenerate solution (see ?robust_lpa, \".force_pd()\" in the source) -- especially for model 3-6 with few `n_starts`. Before trusting this fit, compare against a more constrained `model` (1-2), more `n_starts`, or fewer profiles.",
      100 * fit_indices$Min_Size, fit_indices$Model, fit_indices$Profiles, engine
    ), call. = FALSE)
  }
  invisible(NULL)
}