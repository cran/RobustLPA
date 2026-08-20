#' Auxiliary M-Step Function for Robust Estimation
#'
#' Computes a single profile's posterior-probability-weighted mean and
#' covariance matrix, with optional Huber down-weighting of multivariate
#' outliers and optional soft-thresholding (LASSO-type shrinkage) of the
#' mean. Used internally by \code{\link{robust_lpa}} at every EM iteration,
#' for every profile; not intended to be called directly by end users.
#'
#' When \code{robust = TRUE} (the default), robustness weights are computed
#' from Mahalanobis distances to a preliminary (non-robust) weighted
#' mean/covariance estimated directly from \code{data}. For complete data
#' this preliminary estimate is the ordinary weighted mean/covariance; for
#' data with missing values it is obtained via pairwise-available-case FIML
#' estimation (\code{robust_update_fiml_cpp()} with all robustness
#' weights set to 1), so that the Huber cutoff is computed against a
#' Mahalanobis metric that actually reflects the scale and correlation
#' structure of the data rather than an arbitrary fixed matrix.
#'
#' When \code{robust = FALSE}, the Huber down-weighting step is skipped
#' entirely (all robustness weights are fixed at 1) and the profile's mean
#' and covariance are computed in a single ordinary posterior-probability-
#' weighted pass -- i.e. classical (non-robust) EM estimation for a Gaussian
#' mixture, with the same FIML available-case handling of missing data. This
#' also skips the preliminary-estimate pass, so \code{robust = FALSE} is
#' somewhat faster in addition to being non-robust.
#'
#' @param data A numeric matrix, possibly containing \code{NA} values.
#' @param z Posterior probabilities for a given cluster (length \code{nrow(data)}).
#' @param alpha Significance level for the Huber threshold (chi-squared
#'   cutoff): observations whose squared Mahalanobis distance exceeds the
#'   \code{1 - alpha} chi-squared quantile are down-weighted. Smaller
#'   \code{alpha} down-weights fewer, more extreme points; larger \code{alpha}
#'   down-weights more aggressively. Default \code{0.05}. Ignored when
#'   \code{robust = FALSE}.
#' @param lambda Non-negative LASSO penalty applied to the mean vector via
#'   soft-thresholding. Only meaningful on centered/scaled data (shrinkage is
#'   toward zero on the raw variable scale); see \code{\link{robust_lpa}}.
#'   Default \code{0} (no shrinkage).
#' @param robust Logical. If \code{TRUE} (default), use Huber down-weighting
#'   of multivariate outliers as described above. If \code{FALSE}, estimate
#'   the mean and covariance with ordinary (non-robust) posterior-probability
#'   weights only.
#' @return A list with \code{mean} (numeric row vector) and \code{covariance}
#'   (a positive-(semi)definite matrix).
#' @keywords internal
robust_m_step <- function(data, z, alpha = 0.05, lambda = 0, robust = TRUE) {
  has_na <- any(is.na(data))
  n <- nrow(data)
  p <- ncol(data)
  
  if (n == 0 || p == 0) stop("`data` must have at least one row and one column.")
  if (length(z) != n) stop("`z` must have length equal to `nrow(data)`.")
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be strictly between 0 and 1.")
  if (lambda < 0) stop("`lambda` must be non-negative.")
  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) stop("`robust` must be a single TRUE/FALSE value.")
  
  sum_z <- max(sum(z), 1e-6)
  valid_weights <- z / sum_z
  
  if (!robust) {
    # Classical (non-robust) weighted estimation: uniform robustness weights
    # (w = 1) in a single pass, still using FIML available-case handling for
    # missing data. This recovers standard (non-robust) EM estimation for a
    # Gaussian mixture.
    w_flat <- rep(1, n)
    if (!has_na) {
      return(robust_update_cpp(data, z, w_flat, lambda))
    } else {
      return(robust_update_fiml_cpp(data, z, w_flat, lambda))
    }
  }
  
  if (!has_na) {
    mu_init <- colSums(data * valid_weights)
    diff <- sweep(data, 2, mu_init, "-")
    sigma_init <- t(diff) %*% (diff * valid_weights)
    sigma_init <- (sigma_init + t(sigma_init)) / 2
    
    ev <- eigen(sigma_init, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-5) {
      sigma_init <- sigma_init + diag(abs(min(ev)) + 1e-4, p)
    }
    
    dists <- mahalanobis_cpp(data, mu_init, sigma_init)
    cutoff <- qchisq(1 - alpha, df = p)
    w <- huber_weights_cpp(dists, cutoff)
    return(robust_update_cpp(data, z, w, lambda))
  } else {
    # Preliminary (non-robust, w = 1) pairwise-available-case FIML mean and
    # covariance. This replaces a previous version of this function that
    # used a fixed identity matrix here, which made the Mahalanobis
    # distances (and therefore the Huber weights) blind to the actual scale
    # and correlation structure of the data whenever variables were not
    # already standardized and uncorrelated. Using an actual data-driven
    # preliminary covariance keeps the missing-data path consistent with the
    # complete-data path above.
    prelim <- robust_update_fiml_cpp(data, z, w = rep(1, n), lambda = 0)
    mu_init <- prelim$mean
    sigma_init <- prelim$covariance
    sigma_init <- (sigma_init + t(sigma_init)) / 2
    
    ev <- eigen(sigma_init, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-5) {
      sigma_init <- sigma_init + diag(abs(min(ev)) + 1e-4, p)
    }
    
    dists <- mahalanobis_fiml_cpp(data, mu_init, sigma_init)
    w <- huber_weights_fiml_cpp(data, dists, alpha)
    return(robust_update_fiml_cpp(data, z, w, lambda))
  }
}