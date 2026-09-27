test_that("MCMC and EM estimate the same model 1 (variable-specific shared variances)", {
  set.seed(2)
  cls <- rep(1:2, each = 200)
  Y <- cbind(a = stats::rnorm(400, ifelse(cls == 1, 0, 3), 1),
             b = stats::rnorm(400, ifelse(cls == 1, 0, 3), 4))
  em <- robust_lpa(Y, G = 2, model = 1, robust = FALSE)
  mc <- suppressWarnings(robust_lpa(Y, G = 2, model = 1, engine = "MCMC", mcmc_iter = 300,
                                    n_chains = 2, robust = FALSE))
  expect_equal(unname(diag(mc$covariances[[1]])), unname(diag(em$covariances[[1]])), tolerance = 0.15)
  expect_equal(mc$covariances[[1]], mc$covariances[[2]])
})

test_that("MCMC t mixture with missing data returns a complete, valid fit", {
  data(neuro_data)
  x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop", "RT_TMT")]))
  set.seed(1)
  x[sample(length(x), 40)] <- NA
  set.seed(3)
  fit <- suppressWarnings(robust_lpa(x, G = 2, model = 6, engine = "MCMC", mcmc_iter = 300,
                                     n_chains = 2, robust_method = "t"))
  expect_true(is.numeric(fit$nu) && fit$nu >= 1)
  expect_true("nu" %in% fit$mcmc_diagnostics$Parameter)
  expect_true(is.finite(fit$fit$WAIC))
  expect_true(all(abs(rowSums(fit$probabilities) - 1) < 1e-10))
  expect_s3_class(summary(fit), "summary.robust_lpa")
})

test_that("MCMC models 4 and 5 respect their constraints", {
  data(neuro_data)
  x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop", "RT_TMT")]))
  set.seed(4)
  f4 <- suppressWarnings(robust_lpa(x, G = 2, model = 4, engine = "MCMC", mcmc_iter = 200, n_chains = 1, robust = FALSE))
  set.seed(5)
  f5 <- suppressWarnings(robust_lpa(x, G = 2, model = 5, engine = "MCMC", mcmc_iter = 200, n_chains = 1, robust = FALSE))
  ch4 <- f4$mcmc_draws$chains[[1]]$sigma_chain[[200]]
  ch5 <- f5$mcmc_draws$chains[[1]]$sigma_chain[[200]]
  expect_equal(stats::cov2cor(ch4[[1]]), stats::cov2cor(ch4[[2]]), tolerance = 1e-10)
  expect_equal(diag(ch5[[1]]), diag(ch5[[2]]), tolerance = 1e-10)
})
