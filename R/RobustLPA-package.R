#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib RobustLPA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qchisq rnorm runif
## usethis namespace: end
NULL

# Prevent CRAN notes for ggplot2 variables (Standard Non-Standard Evaluation fix)
utils::globalVariables(c("Variable", "Mean", "Class", "Ymin", "Ymax"))
