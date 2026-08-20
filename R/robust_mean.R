#' Calculate a Simple Robust (Trimmed) Mean
#'
#' Computes a one-step trimmed centroid: an observation is included in the
#' average only if its Euclidean distance to the coordinate-wise median of
#' \code{data} is below \code{threshold}. This is a quick, easy-to-reason-about
#' robust location estimate, not an iterative M-estimator; for the full
#' robust mixture-model estimation used elsewhere in this package, see
#' \code{\link{robust_lpa}}.
#'
#' @details
#' \code{threshold} is a distance \emph{from the coordinate-wise median of
#' \code{data}}, not from the origin -- so a sensible value depends on the
#' scale and spread of your variables. A reasonable starting point is a
#' small multiple of a typical per-variable standard deviation times
#' \code{sqrt(ncol(data))} (roughly the scale of a Euclidean distance across
#' all variables); \code{\link[stats]{mahalanobis}}-based thresholds (as used
#' internally by \code{\link{robust_lpa}}) account for correlation and scale
#' automatically and are preferable when variables are on very different
#' scales.
#'
#' @param data A matrix or data.frame of numeric observations.
#' @param threshold Maximum Euclidean distance to the coordinate-wise median
#'   for an observation to be included in the average. Default \code{10}.
#' @return A numeric vector representing the robust mean of the variables.
#'   If no observation falls within \code{threshold} of the median, returns a
#'   vector of zeros with a warning.
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' r_mean <- robust_mean(x, threshold = 3)
#'
#' # Print the calculated robust means
#' r_mean
#' @export
robust_mean <- function(data, threshold = 10) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a matrix or data.frame.")
  }
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "double"
  if (!is.numeric(data_mat)) stop("`data` must contain only numeric (or coercible-to-numeric) columns.")
  if (nrow(data_mat) < 1 || ncol(data_mat) < 1) stop("`data` must have at least one row and one column.")
  if (any(is.na(data_mat))) stop("`data` must not contain missing values; remove or impute them first.")
  if (!is.numeric(threshold) || length(threshold) != 1 || threshold <= 0) {
    stop("`threshold` must be a single positive number.")
  }
  
  # Calls the hidden C++ engine
  result <- robust_mean_cpp(X = data_mat, threshold = threshold)
  
  return(result)
}