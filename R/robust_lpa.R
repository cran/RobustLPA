#' Fit a Single Robust Latent Profile Analysis Model
#'
#' Estimates a Latent Profile Analysis (finite mixture) model that is robust
#' to multivariate outliers and handles missing data by full-information
#' maximum likelihood, using either an EM or an MCMC (Bayesian Lasso)
#' engine. Two robust estimators are available (see "Robust estimation"
#' below): Huber-type down-weighting (\code{robust_method = "huber"}, the
#' default) and a mixture of multivariate t distributions
#' (\code{robust_method = "t"}), which is a proper likelihood-based model.
#' The MCMC engine runs \code{n_chains} independent chains (4 by default)
#' and reports Gelman-Rubin \eqn{\hat{R}}, effective sample size and WAIC.
#' Set \code{cores > 1} to run the EM engine's random restarts, or the MCMC
#' engine's chains, in parallel.
#'
#' @param data A matrix or data.frame of observations (numeric columns only;
#'   \code{NA} is allowed, see "Missing data").
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
#'   Both engines implement exactly the same six parameterizations. Models 4
#'   and 5 have no closed-form M-step; the EM engine fits them by conditional
#'   maximization steps that never decrease the likelihood.
#' @param engine String. Either \code{"EM"} (default) or \code{"MCMC"}.
#' @param max_iter Maximum number of EM iterations. With \code{engine =
#'   "MCMC"}, used for the preliminary EM fit that initializes the chains.
#' @param tol Tolerance for EM convergence (absolute change of the
#'   observed-data log-likelihood between iterations).
#' @param n_starts Number of EM initializations; the fit with the highest
#'   log-likelihood across starts is returned. With \code{engine = "MCMC"},
#'   the preliminary EM fit uses \code{min(n_starts, 5)} starts.
#' @param lambda Non-negative soft-thresholding (LASSO-type) penalty applied
#'   to the profile means, via direct per-coordinate soft-thresholding of the
#'   (robustness- and posterior-probability-weighted) mean at every M-step.
#'   This shrinks each mean component toward zero on the scale of the (as
#'   supplied) \code{data}; it is only statistically meaningful as a
#'   sparsity-inducing penalty on centered/scaled data, and, because the
#'   threshold is a flat amount rather than one rescaled by each profile's
#'   variance/sample size, it is a computationally convenient approximation
#'   to (not an exact coordinate-wise solution of) the corresponding
#'   L1-penalized log-likelihood except when profile covariances are close
#'   to the identity. Standardize \code{data} first if you intend to use
#'   \code{lambda > 0}; a \code{warning()} is raised if \code{data} does not
#'   look centered/scaled. Default \code{0} (no shrinkage). Ignored by the
#'   MCMC engine (see \code{prior_laplace}).
#' @param mcmc_iter Number of iterations per MCMC chain (see \code{n_chains}).
#'   The first half of each chain is discarded as burn-in before computing
#'   posterior summaries and convergence diagnostics. Ignored when \code{engine = "EM"}.
#' @param prior_laplace Positive numeric, the Laplace (Bayesian Lasso)
#'   shrinkage/rate hyperparameter for the profile means under the MCMC
#'   engine (denoted \eqn{\lambda} in Park & Casella, 2008; \strong{larger}
#'   values induce \strong{more} shrinkage of the means toward zero).
#'   Ignored when \code{engine = "EM"}.
#' @param robust Logical. If \code{TRUE} (default), robust estimation is
#'   used with the method selected by \code{robust_method}. If \code{FALSE},
#'   both engines reduce to classical (non-robust) Gaussian-mixture
#'   estimation with the same interface, missing-data handling and
#'   parameterizations.
#' @param alpha Significance level for the Huber down-weighting threshold
#'   (observations with squared Mahalanobis distance beyond the
#'   \code{1 - alpha} chi-squared quantile are down-weighted). Default
#'   \code{0.05}. Used only when \code{robust = TRUE} and
#'   \code{robust_method = "huber"}.
#' @param n_chains Number of independent MCMC chains to run (default
#'   \code{4}). Ignored when \code{engine = "EM"}.
#' @param cores Integer, number of CPU cores to use for parallel estimation
#'   \emph{within this single} \code{robust_lpa()} call (default \code{1}).
#'   For \code{engine = "EM"}, parallelizes across the \code{n_starts} random
#'   restarts; for \code{engine = "MCMC"}, across the \code{n_chains} chains.
#'   Uses \code{parallel::mclapply()} on Unix-alikes and a PSOCK cluster on
#'   Windows. Results are identical for any value of \code{cores} given the
#'   same \code{set.seed()} (see "Reproducibility").
#' @param robust_method Either \code{"huber"} (default) or \code{"t"}. See
#'   "Robust estimation". Ignored when \code{robust = FALSE}.
#' @param nu Degrees of freedom of the multivariate-t components when
#'   \code{robust_method = "t"}. \code{NULL} (default) estimates a single
#'   \code{nu} shared by all profiles (EM: ECME update maximizing the
#'   observed-data log-likelihood over \code{[1, 200]}; MCMC: Metropolis step
#'   under a Gamma(2, 0.1) prior).
#'   A positive number fixes it (e.g. \code{nu = 4}).
#' @param init Initialization of the EM starts: \code{"kmeans"} (default;
#'   each start runs k-means with different random centers on standardized,
#'   mean-imputed data) or \code{"random"} (random soft partitions, as in
#'   version 1.0.0). Random soft partitions start from almost identical
#'   profiles, so they need several \code{n_starts} and can end in an empty
#'   profile; \code{"kmeans"} is recommended.
#' @return A list with S3 class \code{"robust_lpa"} (see
#'   \code{\link{print.robust_lpa}} and \code{\link{summary.robust_lpa}})
#'   containing:
#'   \describe{
#'     \item{engine}{The estimation engine used.}
#'     \item{robust_method}{\code{"none"}, \code{"huber"} or \code{"t"}.}
#'     \item{means}{A list of length \code{G} with the estimated profile means.}
#'     \item{covariances}{A list of length \code{G} with the estimated profile
#'       covariance matrices (for \code{robust_method = "t"}: the scale
#'       matrices; the covariance of a t component is \code{nu / (nu - 2)}
#'       times the scale matrix for \code{nu > 2}).}
#'     \item{proportions}{Numeric vector of length \code{G} with the estimated mixing proportions.}
#'     \item{nu}{The (estimated or fixed) t degrees of freedom; \code{NULL}
#'       unless \code{robust_method = "t"}. Posterior median for the MCMC engine.}
#'     \item{probabilities}{An \code{n x G} matrix of posterior profile-membership probabilities.}
#'     \item{weights}{Numeric vector of length \code{n}: each observation's
#'       robustness weight, averaged over profiles with its posterior
#'       probabilities; small values flag outliers. Huber weights lie in
#'       (0, 1] (1 = full weight). With \code{robust_method = "t"} they are the
#'       E-step weights \eqn{(\nu + p_i)/(\nu + \delta_i)} (\eqn{p_i}
#'       observed variables, \eqn{\delta_i} squared Mahalanobis distance),
#'       which exceed 1 for observations close to the profile mean. All 1
#'       when \code{robust = FALSE}.}
#'     \item{fit}{A one-row data.frame with \code{Model}, \code{Profiles},
#'       \code{LogLik}, \code{Parameters}, \code{AIC}, \code{BIC}, \code{SABIC},
#'       \code{Entropy}, \code{Min_Size}, \code{Max_Size} and, for the MCMC
#'       engine, \code{WAIC}. \code{Parameters} counts every free parameter;
#'       only with the EM engine and \code{lambda > 0} are the mean components
#'       shrunk exactly to zero left out.}
#'     \item{assignments}{Integer vector of length \code{n} with the most likely profile for each observation.}
#'     \item{converged, iterations}{(EM engine) whether the selected start
#'       met \code{tol} before \code{max_iter}, and how many iterations it used.}
#'     \item{mcmc_draws}{(MCMC engine only) a list with \code{chains} (the raw
#'       per-chain draws, relabeled to a common profile ordering),
#'       \code{n_chains}, \code{mcmc_iter}, \code{burnin}, and \code{nu_estimated}.}
#'     \item{mcmc_diagnostics}{(MCMC engine only) a data.frame with one row per scalar parameter (\code{Parameter}, \code{Rhat}, \code{ESS}).}
#'     \item{waic}{(MCMC engine only) a list with \code{WAIC}, \code{lppd} and \code{p_waic}.}
#'     \item{data}{The numeric matrix actually fit.}
#'     \item{call_args}{A named list of every argument controlling this fit,
#'       used to refit the identical specification on new data (e.g. by
#'       \code{\link{bch_robust}}'s bootstrap correction).}
#'   }
#' @section Robust estimation:
#' \describe{
#'   \item{\code{robust_method = "huber"}}{At every M-step each observation's
#'     contribution to a profile's mean and covariance is down-weighted by a
#'     Huber weight computed from its squared Mahalanobis distance to that
#'     profile's \emph{current} (already robust) estimates, so the
#'     down-weighting accumulates across iterations and converges to an
#'     iteratively reweighted M-estimator. This is an estimating-equation
#'     approach, not a likelihood: the reported \code{LogLik} (and hence
#'     AIC/BIC/SABIC and \code{\link{blrt_robust}}) is the Gaussian mixture
#'     log-likelihood evaluated at the robust estimates, which is a useful
#'     heuristic but not the maximized objective. Huber-type estimators also
#'     have a limited breakdown point: a large, compact cluster of outliers
#'     can still mask itself.}
#'   \item{\code{robust_method = "t"}}{Each profile is a multivariate t
#'     distribution (McLachlan & Peel, 1998; Peel & McLachlan, 2000), whose
#'     heavier tails automatically down-weight outlying observations through
#'     the latent-scale weights \code{(nu + p) / (nu + d)}. This is a proper
#'     likelihood-based model fitted by ECME, so \code{LogLik}, AIC/BIC, the
#'     BLRT and the BCH method are all used as intended; it is also markedly
#'     more resistant than Huber weighting to gross outliers. \code{nu} is
#'     estimated by default and counts as one extra parameter.}
#' }
#' In the MCMC engine, \code{"t"} is implemented as an exact sampler
#' (latent Gamma scales; \code{nu} updated with the latent scales integrated
#' out), whereas \code{"huber"} remains a heuristic that down-weights
#' sufficient statistics and does not target a well-defined posterior;
#' prefer \code{"t"} for Bayesian inference.
#' @section Missing data:
#' Rows with missing values contribute their exact observed-data likelihood
#' (the marginal density of their observed entries) to the E-step, and the
#' M-step uses the exact EM treatment of incomplete multivariate data:
#' missing entries are replaced by their conditional expectations given the
#' observed ones and the conditional covariance is added to the scatter
#' matrix (Ghahramani & Jordan, 1994; Liu & Rubin, 1995). The resulting
#' estimates are maximum likelihood under missing-at-random (MAR)
#' missingness. The MCMC engine uses the equivalent data-augmentation step
#' (missing entries are drawn from their conditional distribution every
#' sweep). Rows with no observed variables are allowed and are allocated
#' according to the mixing proportions only.
#' @section MCMC details:
#' The chains are initialized from a preliminary EM fit (same model and
#' robust method), with the profile means of each chain independently
#' perturbed by half a within-profile standard deviation so that the chains
#' start from dispersed points (as the Gelman-Rubin diagnostic assumes).
#' Priors: Laplace (Bayesian Lasso) on the means, Dirichlet(1, ..., 1) on
#' the mixing proportions, inverse-gamma(1, 1) on variances (models 1, 2 and
#' the scales of models 4-5), inverse-Wishart(p + 1, I) on full covariance
#' matrices (models 3 and 6), a uniform distribution on correlation matrices
#' (models 4-5), and Gamma(2, 0.1) on \code{nu}. All conditionals are
#' sampled exactly except the correlation matrices of models 4-5, which use
#' random-walk Metropolis steps (step sizes adapted during burn-in only).
#' Label switching is resolved by relabeling every draw to the EM solution
#' (pivotal reordering) with an exact optimal assignment on standardized
#' mean distances. \code{$mcmc_diagnostics} reports the Gelman-Rubin
#' \eqn{\hat{R}} and the effective sample size of every scalar parameter;
#' \code{$fit$WAIC} is the widely applicable information criterion
#' (Watanabe, 2010), computed from up to 1000 posterior draws and preferable
#' to the plug-in AIC/BIC for comparing Bayesian fits.
#' @section Reproducibility:
#' Every random quantity (EM starts, MCMC chains, bootstrap replicates) is
#' generated from a seed drawn from R's random number stream before any
#' parallel dispatch, so \code{set.seed(1); robust_lpa(..., cores = 1)} and
#' \code{set.seed(1); robust_lpa(..., cores = 4)} give identical results on
#' the same machine. Across operating systems, compilers or linear-algebra
#' libraries, results may differ in the last digits (and, when two solutions
#' are almost equally good, in the selected start or the order of the
#' profiles).
#' @seealso \code{\link{estimate_profiles_robust}} to fit and compare many
#'   \code{G} / \code{model} combinations at once, \code{\link{blrt_robust}}
#'   for a bootstrapped likelihood ratio test, \code{\link{bch_robust}} to
#'   relate profiles to distal outcomes, and \code{\link{plot_mcmc_chains}}
#'   to inspect MCMC chains.
#' @references
#'   Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#'   using multiple sequences. \emph{Statistical Science}, 7(4), 457-472.
#'   \doi{10.1214/ss/1177011136}
#'
#'   Ghahramani, Z., & Jordan, M. I. (1994). Supervised learning from
#'   incomplete data via an EM approach. \emph{Advances in Neural
#'   Information Processing Systems}, 6, 120-127.
#'
#'   Liu, C., & Rubin, D. B. (1995). ML estimation of the t distribution
#'   using EM and its extensions, ECM and ECME. \emph{Statistica Sinica},
#'   5(1), 19-39.
#'
#'   Peel, D., & McLachlan, G. J. (2000). Robust mixture modelling using the
#'   t distribution. \emph{Statistics and Computing}, 10(4), 339-348.
#'   \doi{10.1023/A:1008981510081}
#'
#'   Park, T., & Casella, G. (2008). The Bayesian Lasso. \emph{Journal of the
#'   American Statistical Association}, 103(482), 681-686.
#'   \doi{10.1198/016214508000000337}
#'
#'   Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#'   and widely applicable information criterion in singular learning
#'   theory. \emph{Journal of Machine Learning Research}, 11, 3571-3594.
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "Attention", "Executive_Functions",
#'                                      "RT_Stroop", "RT_TMT")]))
#' set.seed(1)
#' fit <- robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 50)
#' fit           # print.robust_lpa(): concise overview of the fit
#' summary(fit)  # summary.robust_lpa(): profile means/sizes and fit indices
#'
#' # Multivariate-t mixture: a likelihood-based robust alternative
#' fit_t <- robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 50,
#'                     robust_method = "t")
#' fit_t$nu
#' head(sort(fit_t$weights))  # smallest weights = most outlying observations
#'
#' # Compare with classical (non-robust) estimation
#' fit_classical <- robust_lpa(x, G = 2, model = 6, n_starts = 2, max_iter = 50,
#'                             robust = FALSE)
#' sapply(fit_t$means, `[`, "RT_Stroop")
#' sapply(fit_classical$means, `[`, "RT_Stroop")
#'
#' \donttest{
#' # MCMC engine: 4 chains (default), with a small mcmc_iter for speed
#' fit_mcmc <- robust_lpa(x, G = 2, model = 6, engine = "MCMC",
#'                        robust_method = "t", mcmc_iter = 300, n_chains = 4)
#' summary(fit_mcmc)  # includes Rhat / ESS ranges and WAIC
#' }
#' @export
robust_lpa <- function(data, G, model = 6, engine = "EM", max_iter = 100, tol = 1e-6,
                       n_starts = 5, lambda = 0, mcmc_iter = 2000, prior_laplace = 0.1,
                       robust = TRUE, alpha = 0.05, n_chains = 4, cores = 1,
                       robust_method = c("huber", "t"), nu = NULL,
                       init = c("kmeans", "random")) {

  # ---- input validation -----------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a matrix or data.frame.")
  }
  X <- as.matrix(data)
  if (!is.numeric(X) && !is.logical(X)) {
    stop("`data` must contain only numeric (or coercible-to-numeric) columns.")
  }
  storage.mode(X) <- "double"
  robust_method <- match.arg(robust_method)
  init <- match.arg(init)

  n <- nrow(X)
  p <- ncol(X)
  if (is.null(n) || n < 2) stop("`data` must have at least 2 rows.")
  if (is.null(p) || p < 1) stop("`data` must have at least 1 column.")
  if (!is.numeric(G) || length(G) != 1 || G < 1 || G != round(G)) {
    stop("`G` must be a single positive integer.")
  }
  G <- as.integer(G)
  if (!(length(model) == 1 && model %in% 1:6)) {
    stop("`model` must be a single integer between 1 and 6.")
  }
  model <- as.integer(model)
  if (!(length(engine) == 1 && engine %in% c("EM", "MCMC"))) {
    stop("`engine` must be either 'EM' or 'MCMC'.")
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) stop("`max_iter` must be a positive integer.")
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) stop("`tol` must be a positive number.")
  if (!is.numeric(n_starts) || length(n_starts) != 1 || n_starts < 1) stop("`n_starts` must be at least 1.")
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) stop("`lambda` must be non-negative.")
  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) stop("`robust` must be a single TRUE/FALSE value.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("`alpha` must be a single number strictly between 0 and 1.")
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive integer.")
  }
  if (!is.null(nu) && (!is.numeric(nu) || length(nu) != 1 || !is.finite(nu) || nu <= 0)) {
    stop("`nu` must be NULL (estimate it) or a single positive number.")
  }
  if (G > n) stop("`G` cannot exceed the number of observations.")
  if (all(is.na(X))) stop("`data` contains no observed values.")
  if (engine == "MCMC") {
    if (!is.numeric(mcmc_iter) || length(mcmc_iter) != 1 || mcmc_iter < 2) stop("`mcmc_iter` must be at least 2.")
    if (!is.numeric(prior_laplace) || length(prior_laplace) != 1 || prior_laplace <= 0) stop("`prior_laplace` must be positive.")
    if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains < 1 || n_chains != round(n_chains)) {
      stop("`n_chains` must be a single positive integer.")
    }
  }
  if (lambda > 0) {
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

  call_args <- list(
    G = G, model = model, engine = engine, max_iter = max_iter, tol = tol,
    n_starts = n_starts, lambda = lambda, mcmc_iter = mcmc_iter,
    prior_laplace = prior_laplace, robust = robust, alpha = alpha,
    n_chains = n_chains, cores = cores, robust_method = robust_method,
    nu = nu, init = init
  )

  method <- .robust_method_name(robust, robust_method)
  pat <- .missing_patterns(X)
  var_names <- colnames(X)
  nu_estimated <- identical(method, "t") && is.null(nu)
  extra_params <- as.integer(nu_estimated)

  if (engine == "MCMC") {
    return(.robust_lpa_mcmc(
      X = X, G = G, model = model, max_iter = max_iter, tol = tol, n_starts = n_starts,
      mcmc_iter = mcmc_iter, prior_laplace = prior_laplace, robust = robust, alpha = alpha,
      n_chains = n_chains, cores = cores, robust_method = robust_method, nu = nu, init = init,
      method = method, pat = pat, var_names = var_names, nu_estimated = nu_estimated,
      extra_params = extra_params, call_args = call_args
    ))
  }

  # ---- EM engine --------------------------------------------------------
  # One EM run from one initialization. Self-contained so that it can be
  # dispatched unchanged to parallel workers (see `.run_parallel()`); any
  # log-likelihood-decrease messages are returned rather than raised, so they
  # are never lost on a worker.
  #
  # For the likelihood-based fits ("none" = Gaussian, "t" = multivariate t),
  # every EM/ECM step should increase the observed-data log-likelihood; the
  # best iterate is tracked and returned, and sustained decreases (which can
  # only come from the approximate M-steps of models 4-5, LASSO
  # soft-thresholding, or the positive-definiteness safety net) are
  # reported. The Huber fit is an iteratively reweighted M-estimator whose
  # fixed point is *not* a maximizer of the Gaussian log-likelihood, so the
  # log-likelihood is not monotone by design: the final (converged) iterate
  # is returned and no decrease warnings are issued.
  monotone <- !identical(method, "huber")

  run_one_start <- function(start) {
    z0 <- .initial_z(X, G, init)
    pi_g <- colMeans(z0)
    mu <- vector("list", G)
    raw <- vector("list", G)
    for (g in seq_len(G)) {
      ms <- robust_m_step(X, z0[, g], alpha = alpha, lambda = lambda, robust = FALSE, patterns = pat)
      mu[[g]] <- as.numeric(ms$mean)
      raw[[g]] <- ms$covariance
    }
    sigma <- .apply_covariance_model(raw, pi_g, model, G, p)
    nu_cur <- if (identical(method, "t")) (nu %||% 10) else NULL

    log_lik <- -Inf
    best <- list(log_lik = -Inf, mu = mu, sigma = sigma, pi_g = pi_g, z = NULL, nu = nu_cur)
    n_decreases <- 0L
    max_single_decrease <- 0
    consecutive_decreases <- 0L
    divergence_iter <- NA_integer_
    final_iter <- 0L
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      final_iter <- iter
      lf <- .mixture_logf(X, pat, mu, sigma, pi_g, method, nu_cur)
      post <- .posterior_from_logf(lf$logf)
      z <- post$z
      new_log_lik <- post$loglik

      if (monotone && new_log_lik < log_lik - 1e-4) {
        n_decreases <- n_decreases + 1L
        max_single_decrease <- max(max_single_decrease, log_lik - new_log_lik)
        consecutive_decreases <- consecutive_decreases + 1L
      } else {
        consecutive_decreases <- 0L
      }

      if (!monotone || new_log_lik > best$log_lik) {
        best <- list(log_lik = new_log_lik, mu = mu, sigma = sigma, pi_g = pi_g, z = z, nu = nu_cur)
      }

      if (abs(new_log_lik - log_lik) < tol) {
        converged <- TRUE
        break
      }
      log_lik <- new_log_lik

      if (monotone && consecutive_decreases >= 5L) {
        divergence_iter <- iter
        break
      }

      # ---- M-step (weights from the current parameters, i.e. the ones that
      # produced z; exact EM treatment of missing entries) ----------------
      for (g in seq_len(G)) {
        w_g <- .robust_weights(lf$maha[, g], pat$pobs, method, alpha, nu_cur %||% 4)
        ms <- robust_m_step(X, z[, g], alpha = alpha, lambda = lambda,
                            robust = !identical(method, "none"),
                            mu = mu[[g]], sigma = sigma[[g]],
                            robust_method = if (identical(method, "t")) "t" else "huber",
                            nu = nu_cur %||% 4, weights = w_g, patterns = pat)
        mu[[g]] <- as.numeric(ms$mean)
        raw[[g]] <- ms$covariance
        pi_g[g] <- mean(z[, g])
      }
      sigma <- .apply_covariance_model(raw, pi_g, model, G, p, current = sigma)
      if (nu_estimated) nu_cur <- .update_nu_ecme(X, pat, mu, sigma, pi_g, nu_cur)
    }

    warnings_out <- character(0)
    if (monotone && n_decreases > 0L) {
      cause <- if (lambda > 0) {
        "the LASSO soft-thresholding step, which is an approximate penalized update"
      } else if (model %in% c(3, 6)) {
        "the positive-definiteness safety net (.force_pd()) perturbing a near-singular covariance"
      } else {
        "ordinary numerical noise near convergence"
      }
      divergence_note <- if (!is.na(divergence_iter)) {
        sprintf(
          " Stopped early at iteration %d after 5 consecutive decreases; the best iterate found (log-lik = %.6f) was kept instead.",
          divergence_iter, best$log_lik
        )
      } else {
        sprintf(" The best iterate found (log-lik = %.6f) was kept.", best$log_lik)
      }
      warnings_out <- sprintf(
        "Start %d: log-likelihood decreased in %d of %d iteration(s) (largest single decrease %.6f) for model %d. This can happen with %s; consider comparing against a different `model`, increasing `n_starts`, or checking for an implausibly small profile size if this persists.%s",
        start, n_decreases, final_iter, max_single_decrease, model, cause, divergence_note
      )
    }

    c(best, list(warnings = warnings_out, converged = converged, iterations = final_iter))
  }

  start_results <- .run_parallel(
    cores, n_starts, run_one_start,
    export_vars = c("X", "n", "p", "G", "model", "max_iter", "tol", "lambda", "robust",
                    "alpha", "method", "pat", "nu", "nu_estimated", "init", "monotone"),
    export_env = environment()
  )

  log_liks <- vapply(start_results, function(r) r$log_lik, numeric(1))
  if (!any(is.finite(log_liks))) {
    stop(
      "All ", n_starts, " EM start(s) produced a non-finite log-likelihood (model ", model,
      ", G = ", G, "). This usually indicates a numerically degenerate fit (e.g. too few ",
      "observations for the requested number of profiles/parameters); try a different `model`, ",
      "fewer profiles, or more `n_starts`."
    )
  }
  best_idx <- which.max(log_liks)
  best <- start_results[[best_idx]]
  for (w in best$warnings) warning(w, call. = FALSE)

  .finalize_fit(
    X = X, pat = pat, mu = best$mu, sigma = best$sigma, pi_g = best$pi_g, nu_hat = best$nu,
    method = method, alpha = alpha, model = model, G = G, engine = "EM",
    var_names = var_names, extra_params = extra_params, call_args = call_args,
    extra = list(converged = best$converged, iterations = best$iterations)
  )
}

#' Initial Posterior Probabilities for One EM Start
#'
#' @param X Numeric data matrix (\code{NA} allowed).
#' @param G Number of profiles.
#' @param init \code{"kmeans"} or \code{"random"}.
#' @return An \code{n x G} matrix whose rows sum to 1.
#' @keywords internal
#' @noRd
.initial_z <- function(X, G, init) {
  n <- nrow(X)
  if (G == 1) return(matrix(1, n, 1))
  if (identical(init, "kmeans")) {
    Xi <- X
    cm <- colMeans(Xi, na.rm = TRUE)
    cm[!is.finite(cm)] <- 0
    if (anyNA(Xi)) {
      idx <- which(is.na(Xi), arr.ind = TRUE)
      Xi[idx] <- cm[idx[, 2]]
    }
    sds <- apply(Xi, 2, stats::sd)
    sds[!is.finite(sds) | sds < 1e-12] <- 1
    Xi <- sweep(sweep(Xi, 2, cm, "-"), 2, sds, "/")
    km <- tryCatch(
      suppressWarnings(stats::kmeans(Xi, centers = G, nstart = 1, iter.max = 50)),
      error = function(e) NULL
    )
    if (!is.null(km) && length(unique(km$cluster)) == G) {
      z <- matrix(0.1 / (G - 1), n, G)
      z[cbind(seq_len(n), km$cluster)] <- 0.9
      return(z)
    }
  }
  z <- matrix(stats::runif(n * G), n, G)
  z / rowSums(z)
}

#' Assemble a Fitted robust_lpa Object From Final Parameters
#'
#' Shared by both engines: recomputes the posterior probabilities,
#' log-likelihood, robustness weights and fit indices at the final
#' parameters, reattaches variable names, and builds the S3 object.
#'
#' @keywords internal
#' @noRd
.finalize_fit <- function(X, pat, mu, sigma, pi_g, nu_hat, method, alpha, model, G, engine,
                          var_names, extra_params, call_args, extra = list()) {
  n <- nrow(X)
  p <- ncol(X)
  mu <- lapply(mu, as.numeric)
  lf <- .mixture_logf(X, pat, mu, sigma, pi_g, method, nu_hat)
  post <- .posterior_from_logf(lf$logf)
  z <- post$z

  W <- matrix(1, n, G)
  if (!identical(method, "none")) {
    for (g in seq_len(G)) W[, g] <- .robust_weights(lf$maha[, g], pat$pobs, method, alpha, nu_hat %||% 4)
  }
  obs_weights <- rowSums(z * W)

  if (!is.null(var_names)) {
    for (g in seq_len(G)) {
      names(mu[[g]]) <- var_names
      dimnames(sigma[[g]]) <- list(var_names, var_names)
    }
  }

  penalized <- identical(engine, "EM") && isTRUE(call_args$lambda > 0)
  fi <- .compute_fit_indices(mu, post$loglik, z, n, p, G, model, extra_params = extra_params,
                             penalized = penalized)
  .warn_if_degenerate_profile(fi$fit_indices, engine)

  out <- c(
    list(
      engine = engine, robust_method = method, means = mu, covariances = sigma,
      proportions = as.numeric(pi_g), nu = if (identical(method, "t")) nu_hat else NULL,
      probabilities = z, weights = obs_weights, fit = fi$fit_indices,
      assignments = fi$assignments
    ),
    extra,
    list(data = X, call_args = call_args)
  )
  class(out) <- "robust_lpa"
  out
}

#' MCMC Branch of robust_lpa()
#'
#' @keywords internal
#' @noRd
.robust_lpa_mcmc <- function(X, G, model, max_iter, tol, n_starts, mcmc_iter, prior_laplace,
                             robust, alpha, n_chains, cores, robust_method, nu, init,
                             method, pat, var_names, nu_estimated, extra_params, call_args) {
  n <- nrow(X)
  p <- ncol(X)

  # Preliminary EM fit (same model and robust method): starting values for
  # every chain and pivot for relabeling.
  init_fit <- suppressWarnings(robust_lpa(
    X, G = G, model = model, engine = "EM", max_iter = max_iter, tol = tol,
    n_starts = max(1, min(n_starts, 5)), lambda = 0, robust = robust, alpha = alpha,
    cores = 1, robust_method = robust_method, nu = nu, init = init
  ))
  init_mu <- lapply(init_fit$means, as.numeric)
  init_sigma <- lapply(init_fit$covariances, unname)
  init_pi <- pmax(init_fit$proportions, 1e-3)
  pooled_var <- Reduce(`+`, Map(function(S, w) w * diag(S), init_sigma, init_fit$proportions))
  pivot_scale <- sqrt(pmax(pooled_var, 1e-8))
  nu_start <- if (identical(method, "t")) (nu %||% init_fit$nu %||% 10) else 30
  robust_type <- switch(method, none = 0L, huber = 1L, t = 2L)

  run_one_chain <- function(chain_id) {
    mu0 <- lapply(seq_len(G), function(g) {
      init_mu[[g]] + stats::rnorm(p, 0, 0.5 * sqrt(pmax(diag(init_sigma[[g]]), 1e-8)))
    })
    mcmc_chain_cpp(
      X = X, G = G, model = model, mcmc_iter = mcmc_iter, prior_laplace = prior_laplace,
      robust_type = robust_type, alpha = alpha, nu_init = nu_start, estimate_nu = nu_estimated,
      init_mu = mu0, init_sigma = init_sigma, init_pi = init_pi,
      pat_obs = pat$obs, pat_rows = pat$rows
    )
  }
  chains_raw <- .run_parallel(
    cores, n_chains, run_one_chain,
    export_vars = c("X", "G", "p", "model", "mcmc_iter", "prior_laplace", "robust_type", "alpha",
                    "nu_start", "nu_estimated", "init_mu", "init_sigma", "init_pi", "pat"),
    export_env = environment()
  )
  for (k in seq_along(chains_raw)) chains_raw[[k]]$nu_estimated <- nu_estimated

  chains_raw <- .relabel_mcmc_chains(chains_raw, G, pivot_means = init_mu, scale = pivot_scale)

  burnin <- floor(mcmc_iter / 2)
  valid_iters <- (burnin + 1):mcmc_iter
  n_pooled <- length(valid_iters) * n_chains

  pi_draws <- do.call(rbind, lapply(chains_raw, function(ch) ch$pi_chain[valid_iters, , drop = FALSE]))
  pi_g <- colMeans(pi_draws)

  mu <- vector("list", G)
  sigma <- vector("list", G)
  for (g in seq_len(G)) {
    mu_sum <- numeric(p)
    sigma_sum <- matrix(0, p, p)
    for (ch in chains_raw) {
      for (it in valid_iters) {
        mu_sum <- mu_sum + as.numeric(ch$mu_chain[[it]][[g]])
        sigma_sum <- sigma_sum + ch$sigma_chain[[it]][[g]]
      }
    }
    mu[[g]] <- mu_sum / n_pooled
    sigma[[g]] <- sigma_sum / n_pooled
  }
  nu_hat <- if (identical(method, "t")) {
    stats::median(unlist(lapply(chains_raw, function(ch) ch$nu_chain[valid_iters])))
  } else {
    NULL
  }

  waic <- .compute_waic(chains_raw, valid_iters, X, pat, method)

  out <- .finalize_fit(
    X = X, pat = pat, mu = mu, sigma = sigma, pi_g = pi_g, nu_hat = nu_hat,
    method = method, alpha = alpha, model = model, G = G, engine = "MCMC",
    var_names = var_names, extra_params = extra_params, call_args = call_args
  )
  out$fit$WAIC <- waic$WAIC

  diagnostics <- .compute_mcmc_diagnostics(chains_raw, mcmc_iter = mcmc_iter, burnin = burnin)
  if (!is.null(diagnostics) && isTRUE(any(diagnostics$Rhat > 1.1, na.rm = TRUE))) {
    warning(
      "Some MCMC parameters have a Gelman-Rubin R-hat > 1.1, suggesting the chains have not ",
      "fully converged. Consider increasing `mcmc_iter`. See `model$mcmc_diagnostics`.",
      call. = FALSE
    )
  }

  out$mcmc_draws <- list(chains = chains_raw, n_chains = n_chains, mcmc_iter = mcmc_iter,
                         burnin = burnin, nu_estimated = nu_estimated)
  out$mcmc_diagnostics <- diagnostics
  out$waic <- waic
  # keep `data` and `call_args` last, as in EM fits
  out <- out[c(setdiff(names(out), c("data", "call_args")), "data", "call_args")]
  class(out) <- "robust_lpa"
  out
}

#' Widely Applicable Information Criterion for an MCMC Fit
#'
#' @param chains Relabeled chains.
#' @param valid_iters Post-burn-in iteration indices.
#' @param X Data matrix.
#' @param pat Missingness patterns.
#' @param method Robust method (density family).
#' @param max_draws Maximum number of posterior draws used (thinned evenly).
#' @return A list with \code{WAIC}, \code{lppd}, \code{p_waic}.
#' @keywords internal
#' @noRd
.compute_waic <- function(chains, valid_iters, X, pat, method, max_draws = 1000) {
  draws <- expand.grid(iter = valid_iters, chain = seq_along(chains))
  if (nrow(draws) > max_draws) {
    draws <- draws[unique(round(seq(1, nrow(draws), length.out = max_draws))), , drop = FALSE]
  }
  S <- nrow(draws)
  ll <- matrix(NA_real_, S, nrow(X))
  for (s in seq_len(S)) {
    ch <- chains[[draws$chain[s]]]
    it <- draws$iter[s]
    lf <- .mixture_logf(X, pat, ch$mu_chain[[it]], ch$sigma_chain[[it]], ch$pi_chain[it, ],
                        method, ch$nu_chain[it])
    ll[s, ] <- .row_logsumexp(lf$logf)
  }
  m <- apply(ll, 2, max)
  lppd <- sum(m + log(colMeans(exp(sweep(ll, 2, m, "-")))))
  p_waic <- if (S > 1) sum(apply(ll, 2, stats::var)) else NA_real_
  list(WAIC = -2 * (lppd - p_waic), lppd = lppd, p_waic = p_waic, n_draws = S)
}
