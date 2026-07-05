#' Calculate the robust mean
#'
#' @param data A matrix or data.frame
#' @param threshold Maximum distance allowed to not be considered an outlier
#' @return A numeric vector representing the robust mean of the variables.
#' @examples
#' data(iris)
#' r_mean <- robust_mean(iris[1:30, 1:2])
#'
#' # Print the calculated robust means
#' r_mean
#' @export
robust_mean <- function(data, threshold = 10) {
  data_mat <- as.matrix(data)

  # Calls the hidden C++ engine
  result <- robust_mean_cpp(X = data_mat, threshold = threshold)

  return(result)
}
