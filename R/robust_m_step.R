#' Auxiliary M-Step Function for Robust Estimation
#'
#' Computes a single profile's updated mean and (unconstrained) covariance
#' matrix from its posterior membership probabilities, with optional
#' robustness weighting and optional soft-thresholding (LASSO-type
#' shrinkage) of the mean. Used internally by \code{\link{robust_lpa}} at
#' every EM iteration, for every profile; not intended to be called directly
#' by end users.
#'
#' @section Robustness weights:
#' The robustness weights are computed from each observation's squared
#' Mahalanobis distance (on its observed entries) to the profile's
#' \emph{current} mean/covariance, i.e. the parameters that produced the
#' posterior probabilities \code{z}. Because these current parameters are
#' themselves the robust estimates from the previous iteration, the
#' down-weighting compounds across EM iterations and converges to a
#' fixed point (an iteratively reweighted M-estimator), instead of being
#' recomputed each time from a non-robust starting point (which lets
#' outliers mask themselves by inflating the covariance they are measured
#' against).
#' \describe{
#'   \item{\code{robust_method = "huber"}}{Huber weights: 1 inside the
#'     \code{1 - alpha} chi-squared quantile (with degrees of freedom equal
#'     to the number of observed variables of the row), \code{sqrt(cutoff / d)}
#'     beyond it. The covariance is the Huber-weighted scatter matrix.}
#'   \item{\code{robust_method = "t"}}{The E-step expectation of the latent
#'     scale of a multivariate-t distribution, \code{(nu + p_obs) / (nu + d)};
#'     together with the covariance update below this is the exact ECM
#'     M-step of a multivariate-t mixture (McLachlan & Peel, 2000).}
#' }
#'
#' @section Missing data:
#' Missing entries are handled by the exact EM treatment of incomplete
#' multivariate data (Ghahramani & Jordan, 1994; Liu & Rubin, 1995): each
#' missing block is replaced by its conditional expectation given the
#' observed entries under the profile's current parameters, and the
#' corresponding conditional covariance is added to the scatter matrix.
#' This gives maximum-likelihood estimates under missing-at-random (MAR)
#' missingness. When no current parameters are supplied (initialization),
#' pairwise available-case moments are used as the starting point for that
#' single step.
#'
#' @param data A numeric matrix, possibly containing \code{NA} values.
#' @param z Posterior probabilities for a given profile (length \code{nrow(data)}).
#' @param alpha Significance level for the Huber threshold. Default
#'   \code{0.05}. Used only when \code{robust = TRUE} and
#'   \code{robust_method = "huber"}.
#' @param lambda Non-negative LASSO penalty applied to the mean vector via
#'   soft-thresholding. Only meaningful on centered/scaled data; see
#'   \code{\link{robust_lpa}}. Default \code{0} (no shrinkage).
#' @param robust Logical. If \code{FALSE}, all robustness weights are 1
#'   (classical Gaussian EM).
#' @param mu,sigma The profile's current mean vector and covariance matrix.
#'   If \code{NULL} (default), a preliminary estimate is computed from
#'   \code{z} (pairwise available-case moments, all weights 1).
#' @param robust_method Either \code{"huber"} (default) or \code{"t"}; see
#'   the "Robustness weights" section.
#' @param nu Degrees of freedom of the multivariate t (used only when
#'   \code{robust_method = "t"}). Default \code{4}.
#' @param weights Optional precomputed robustness weights (length
#'   \code{nrow(data)}); if supplied they are used as-is. Used internally to
#'   avoid recomputing distances already available from the E-step.
#' @param patterns Optional precomputed missingness-pattern structure
#'   (internal use).
#' @return A list with \code{mean} (numeric row vector), \code{covariance}
#'   (a symmetric matrix) and \code{weights} (the robustness weights used).
#' @references
#'   Ghahramani, Z., & Jordan, M. I. (1994). Supervised learning from
#'   incomplete data via an EM approach. \emph{Advances in Neural Information
#'   Processing Systems}, 6, 120-127.
#'
#'   Liu, C., & Rubin, D. B. (1995). ML estimation of the t distribution using
#'   EM and its extensions, ECM and ECME. \emph{Statistica Sinica}, 5(1), 19-39.
#'
#'   McLachlan, G. J., & Peel, D. (2000). \emph{Finite Mixture Models}. Wiley.
#'   \doi{10.1002/0471721182}
#' @keywords internal
robust_m_step <- function(data, z, alpha = 0.05, lambda = 0, robust = TRUE,
                          mu = NULL, sigma = NULL, robust_method = c("huber", "t"),
                          nu = 4, weights = NULL, patterns = NULL) {
  robust_method <- match.arg(robust_method)
  data <- as.matrix(data)
  storage.mode(data) <- "double"
  n <- nrow(data)
  p <- ncol(data)

  if (n == 0 || p == 0) stop("`data` must have at least one row and one column.")
  if (length(z) != n) stop("`z` must have length equal to `nrow(data)`.")
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be strictly between 0 and 1.")
  if (lambda < 0) stop("`lambda` must be non-negative.")
  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) stop("`robust` must be a single TRUE/FALSE value.")

  pat <- patterns %||% .missing_patterns(data)
  method <- .robust_method_name(robust, robust_method)

  if (is.null(mu) || is.null(sigma)) {
    prelim <- pairwise_moments_cpp(data, as.numeric(z))
    mu <- as.numeric(prelim$mean)
    sigma <- .force_pd(prelim$covariance)
    if (is.null(weights)) weights <- rep(1, n)
  }

  if (is.null(weights)) {
    es <- class_estep_cpp(data, as.numeric(mu), sigma, pat$obs, pat$rows, 0L, 1)
    weights <- .robust_weights(es$maha, pat$pobs, method, alpha, nu)
  }

  upd <- class_mstep_cpp(data, as.numeric(z), as.numeric(weights), as.numeric(mu), sigma,
                         pat$obs, pat$rows, lambda, identical(method, "t"))
  list(mean = upd$mean, covariance = upd$covariance, weights = weights)
}
