data(neuro_data)
x_neuro <- scale(as.matrix(neuro_data[, c("Memory", "Attention", "Executive_Functions",
                                          "RT_Stroop", "RT_TMT")]))

test_that("classical EM reaches the maximum-likelihood solution (reference: mclust)", {
  # Reference log-likelihoods from mclust::Mclust(x_neuro, G = 2, modelNames = ...)
  ref <- c(`1` = -1469.50, `2` = -1427.96, `3` = -1339.35, `6` = -1267.01)
  for (m in c(1, 2, 3, 6)) {
    set.seed(1)
    fit <- robust_lpa(x_neuro, G = 2, model = m, n_starts = 3, robust = FALSE)
    expect_lt(abs(fit$fit$LogLik - unname(ref[as.character(m)])), 0.05)
  }
})

test_that("fitted covariances satisfy the constraints of every model", {
  set.seed(2)
  xm <- x_neuro
  xm[sample(length(xm), 60)] <- NA
  offdiag <- function(S) S[upper.tri(S)]
  for (m in 1:6) {
    set.seed(3)
    fit <- robust_lpa(xm, G = 2, model = m, n_starts = 2, robust = FALSE)
    S1 <- fit$covariances[[1]]
    S2 <- fit$covariances[[2]]
    if (m %in% c(1, 2)) {
      expect_true(all(abs(offdiag(S1)) < 1e-12))
      expect_true(all(abs(offdiag(S2)) < 1e-12))
    }
    if (m %in% c(1, 3)) expect_equal(S1, S2)
    if (m == 4) expect_equal(stats::cov2cor(S1), stats::cov2cor(S2), tolerance = 1e-6)
    if (m == 5) expect_equal(diag(S1), diag(S2), tolerance = 1e-6)
    expect_true(all(abs(rowSums(fit$probabilities) - 1) < 1e-10))
  }
})

test_that("the likelihood-based EM is monotone for models 4 and 5", {
  for (m in 4:5) {
    set.seed(4)
    expect_no_warning(fit <- robust_lpa(x_neuro, G = 2, model = m, n_starts = 2, robust = FALSE))
    expect_true(fit$converged)
  }
})

test_that("missing values are handled by maximum likelihood under MAR", {
  set.seed(1)
  n <- 2000
  x1 <- stats::rnorm(n)
  x2 <- 0.8 * x1 + 0.6 * stats::rnorm(n)
  X <- cbind(x1 = x1, x2 = x2)
  X[x1 > 0.3 & stats::runif(n) < 0.8, "x2"] <- NA   # missingness depends on observed x1
  fit <- robust_lpa(X, G = 1, model = 6, n_starts = 1, robust = FALSE, tol = 1e-10, max_iter = 500)
  # Closed-form ML under MAR for this monotone pattern: regress x2 on x1 in
  # the complete cases and combine with the ML moments of x1 (all rows).
  ok <- !is.na(X[, "x2"])
  reg <- stats::lm(x2 ~ x1, data = as.data.frame(X), subset = ok)
  b <- stats::coef(reg)
  s2 <- mean(stats::resid(reg)^2)
  m1 <- mean(X[, "x1"])
  v1 <- mean((X[, "x1"] - m1)^2)
  expect_lt(abs(fit$means[[1]][["x2"]] - (b[[1]] + b[[2]] * m1)), 1e-4)
  expect_lt(abs(fit$covariances[[1]]["x1", "x2"] - b[[2]] * v1), 1e-4)
  expect_lt(abs(fit$covariances[[1]]["x2", "x2"] - (s2 + b[[2]]^2 * v1)), 1e-4)
  expect_lt(abs(unname(fit$means[[1]]["x2"])), 0.06)   # truth is 0
})

test_that("rows with no observed values are allowed", {
  xm <- x_neuro
  xm[3, ] <- NA
  set.seed(5)
  fit <- robust_lpa(xm, G = 2, model = 2, n_starts = 2)
  expect_equal(unname(fit$probabilities[3, ]), fit$proportions, tolerance = 1e-8)
})

test_that("posterior probabilities of extreme outliers are valid (no underflow)", {
  set.seed(3)
  Z <- rbind(matrix(stats::rnorm(1200), 300, 4), matrix(stats::rnorm(1200, 2), 300, 4), rep(60, 4))
  for (method in c("huber", "t")) {
    fit <- robust_lpa(Z, G = 2, model = 2, n_starts = 2, robust_method = method)
    expect_equal(sum(fit$probabilities[601, ]), 1)
    expect_true(all(is.finite(fit$probabilities)))
    expect_lt(fit$weights[601], 0.2)
  }
})

test_that("robust estimators resist gross outliers; the t mixture most of all", {
  set.seed(6)
  X <- matrix(stats::rnorm(1200), 400, 3)
  idx <- sample(400, 20)
  X[idx, ] <- X[idx, ] + matrix(stats::rnorm(60, 0, 15), 20, 3)
  classical <- robust_lpa(X, G = 1, model = 6, n_starts = 1, robust = FALSE)
  huber <- robust_lpa(X, G = 1, model = 6, n_starts = 1, robust_method = "huber")
  tfit <- robust_lpa(X, G = 1, model = 6, n_starts = 1, robust_method = "t")
  expect_lt(max(diag(huber$covariances[[1]])), min(diag(classical$covariances[[1]])))
  expect_lt(max(diag(tfit$covariances[[1]])), 1.3)
  expect_lt(tfit$nu, 10)
  expect_equal(tfit$fit$Parameters, 3 + 6 + 1)   # means + covariance + nu
})

test_that("results are reproducible with set.seed() for any number of cores", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(42)
  a <- robust_lpa(x_neuro, G = 3, model = 2, n_starts = 4, cores = 1, init = "random")
  set.seed(42)
  b <- robust_lpa(x_neuro, G = 3, model = 2, n_starts = 4, cores = 2, init = "random")
  expect_identical(a$means, b$means)
  expect_identical(a$probabilities, b$probabilities)
})

test_that("input validation", {
  expect_error(robust_lpa(x_neuro, G = 0), "`G`")
  expect_error(robust_lpa(x_neuro, G = 2, model = 7), "`model`")
  expect_error(robust_lpa(x_neuro, G = 2, nu = -1, robust_method = "t"), "`nu`")
  expect_error(robust_lpa(x_neuro, G = 2, robust_method = "tukey"))
})

test_that("means close to zero are counted as parameters unless the LASSO set them to zero", {
  # One profile on standardized data: every mean is ~0 but is still estimated.
  xs <- scale(as.matrix(neuro_data[, c("Memory", "Attention", "Executive_Functions", "RT_Stroop", "RT_TMT")]))
  p <- ncol(xs)
  f1 <- robust_lpa(xs, G = 1, model = 1, robust = FALSE, n_starts = 1)
  f6 <- robust_lpa(xs, G = 1, model = 6, robust = FALSE, n_starts = 1)
  expect_equal(f1$fit$Parameters, 2 * p)
  expect_equal(f6$fit$Parameters, p + p * (p + 1) / 2)
  expect_equal(f1$fit$BIC, -2 * f1$fit$LogLik + 2 * p * log(nrow(xs)))
  # With the EM LASSO, means shrunk exactly to zero are not free parameters.
  set.seed(1)
  fl <- suppressWarnings(robust_lpa(xs, G = 2, model = 1, robust = FALSE, n_starts = 2, lambda = 0.3))
  n_zero <- sum(vapply(fl$means, function(m) sum(m == 0), numeric(1)))
  expect_gt(n_zero, 0)
  expect_equal(fl$fit$Parameters, 2 * p + p + 1 - n_zero)
})
