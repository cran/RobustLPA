#' Perform a Robust M-Step for Latent Profile Analysis
#'
#' @param data A matrix or data.frame of observations.
#' @param z A numeric vector containing the posterior probabilities of belonging to this cluster.
#' @param alpha Significance level for the Chi-squared outlier cutoff (default is 0.05).
#' @return A list containing the robust 'mean' vector and the robust 'covariance' matrix.
#' @examples
#' data(iris)
#' # Simulate initial probabilities for a single profile
#' z_init <- runif(30)
#' z_init <- z_init / sum(z_init)
#'
#' m_step_res <- robust_m_step(iris[1:30, 1:2], z = z_init)
#'
#' # Print the robust mean calculated in the M-step
#' m_step_res$mean
#' @export
robust_m_step <- function(data, z, alpha = 0.05) {
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  # Helper function to guarantee a matrix is positive-definite
  force_pd <- function(mat) {
    mat <- (mat + t(mat)) / 2 # Ensure perfect symmetry
    ev <- eigen(mat, symmetric = TRUE)$values
    min_ev <- min(ev)
    if (min_ev < 1e-5) mat <- mat + diag(abs(min_ev) + 1e-4, nrow(mat))
    return(mat)
  }

  has_na <- any(is.na(X))

  if (!has_na) {
    # -----------------------------------------------
    # FAST PATH WITHOUT MISSING DATA (Standard C++)
    # -----------------------------------------------
    initial_weights <- rep(1, n)
    baseline <- robust_update_cpp(X = X, z = z, w = initial_weights)

    chi_sq_cutoff <- qchisq(1 - alpha, df = p)
    distances <- mahalanobis_cpp(X = X, mu = baseline$mean, Sigma = baseline$covariance)
    robust_weights <- huber_weights_cpp(squared_dists = distances, chi_sq_cutoff = chi_sq_cutoff)

    robust_results <- robust_update_cpp(X = X, z = z, w = robust_weights)
    return(list(mean = as.numeric(robust_results$mean), covariance = force_pd(robust_results$covariance)))

  } else {
    # -----------------------------------------------
    # PATH WITH MISSING DATA (FIML Optimized in C++)
    # -----------------------------------------------
    initial_weights <- rep(1, n)

    # FIML Baseline via C++
    baseline <- robust_update_fiml_cpp(X = X, z = z, w = initial_weights)
    baseline_cov_pd <- force_pd(baseline$covariance)

    # FIML Mahalanobis and Huber via C++
    distances <- mahalanobis_fiml_cpp(X = X, mu = baseline$mean, Sigma = baseline_cov_pd)
    robust_weights <- huber_weights_fiml_cpp(X = X, squared_dists = distances, alpha = alpha)

    # FIML Robust Update via C++
    robust_results <- robust_update_fiml_cpp(X = X, z = z, w = robust_weights)

    return(list(mean = as.numeric(robust_results$mean), covariance = force_pd(robust_results$covariance)))
  }
}
