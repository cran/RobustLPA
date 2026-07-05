#' Bootstrapped Likelihood Ratio Test for Robust LPA
#'
#' Compares a robust LPA model with G profiles against a null model with G-1 profiles
#' using parametric bootstrapping. Supports FIML simulation conditions.
#'
#' @param data A matrix or data.frame.
#' @param G The number of profiles for the alternative hypothesis (compared against G-1).
#' @param model An integer (1 to 6) specifying the variance-covariance parameterization.
#' @param n_samples Number of bootstrap samples. Default is 50 for speed; 200+ is recommended for publications.
#' @param n_starts Number of starts for the EM algorithm execution.
#' @return A list containing the observed LRT statistic, the vector of bootstrap replicates, and the empirical p-value.
#' @examples
#' # Fast demonstration of the robust BLRT
#' data(iris)
#' blrt_res <- blrt_robust(iris[1:30, 1:2], G = 2, model = 1, n_samples = 2, n_starts = 1)
#' # Print the summary of the results
#' blrt_res
#' @export
blrt_robust <- function(data, G, model = 6, n_samples = 50, n_starts = 2) {
  if (G <= 1) stop("G must be greater than 1 to perform the Bootstrapped Likelihood Ratio Test.")

  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  message(paste0("Executing Robust BLRT: Comparing ", G-1, " vs ", G, " Profiles (Model ", model, ")..."))

  # Estimation of empirical models on real data
  mod_null <- robust_lpa(data = X, G = G - 1, model = model, n_starts = n_starts)
  mod_alt  <- robust_lpa(data = X, G = G, model = model, n_starts = n_starts)

  if (is.null(mod_null) || is.null(mod_alt)) {
    stop("One of the baseline models failed to converge. Unable to compute BLRT.")
  }

  # Calculation of the observed empirical log-likelihood ratio
  lrt_obs <- -2 * (mod_null$fit$LogLik - mod_alt$fit$LogLik)
  if (lrt_obs < 0) lrt_obs <- 0

  lrt_boot <- numeric(n_samples)

  # Parametric stochastic generator based on null model parameters
  simulate_null_mixture <- function(n, p, props, means, covs) {
    sim_X <- matrix(NA, nrow = n, ncol = p)
    assigned_classes <- sample(1:(G-1), size = n, replace = TRUE, prob = props)

    for (i in 1:n) {
      c_i <- assigned_classes[i]
      # Cholesky decomposition to generate stable multivariate distributions
      z_norm <- rnorm(p)
      L <- chol(covs[[c_i]] + diag(1e-6, p))
      sim_X[i, ] <- means[[c_i]] + as.vector(t(L) %*% z_norm)
    }
    return(sim_X)
  }

  # Bootstrap resampling cycle
  for (b in 1:n_samples) {
    message(paste0("Bootstrap Sample ", b, " / ", n_samples))

    # Generation of data under H0
    sim_data <- simulate_null_mixture(n, p, mod_null$proportions, mod_null$means, mod_null$covariances)

    # If the initial data contained NA, replicate the same missingness structure for FIML
    if (any(is.na(X))) {
      sim_data[is.na(X)] <- NA
    }

    # Fit of models on simulated data
    fit_b_null <- tryCatch(robust_lpa(data = sim_data, G = G - 1, model = model, n_starts = n_starts), error = function(e) NULL)
    fit_b_alt  <- tryCatch(robust_lpa(data = sim_data, G = G, model = model, n_starts = n_starts), error = function(e) NULL)

    if (!is.null(fit_b_null) && !is.null(fit_b_alt)) {
      val <- -2 * (fit_b_null$fit$LogLik - fit_b_alt$fit$LogLik)
      lrt_boot[b] <- if (val < 0) 0 else val
    } else {
      lrt_boot[b] <- NA
    }
  }

  # Cleaning vectors from potential internal convergence failures
  lrt_boot <- lrt_boot[!is.na(lrt_boot)]

  # Calculation of empirical p-value
  p_value <- sum(lrt_boot >= lrt_obs) / length(lrt_boot)
  message("BLRT Execution Successfully Completed.")

  return(list(
    LRT_Observed = lrt_obs,
    Bootstrap_LRTs = lrt_boot,
    p_value = p_value
  ))
}
