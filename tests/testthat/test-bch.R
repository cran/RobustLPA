simulate_bch <- function(n, props, sep, mu_y, seed) {
  set.seed(seed)
  cls <- sample(seq_along(props), n, TRUE, prob = props)
  x <- cbind(stats::rnorm(n, ifelse(cls == 2, sep, 0)), stats::rnorm(n, ifelse(cls == 2, sep, 0)))
  y <- stats::rnorm(n, mu_y[cls])
  fit <- robust_lpa(x, G = 2, model = 1, n_starts = 3, robust = FALSE)
  ord <- order(sapply(fit$means, `[`, 1))
  fit$probabilities <- fit$probabilities[, ord]
  fit$assignments <- match(fit$assignments, ord)
  fit$means <- fit$means[ord]
  fit$covariances <- fit$covariances[ord]
  fit$proportions <- fit$proportions[ord]
  list(fit = fit, y = y)
}

test_that("BCH classification error matrix has rows summing to 1", {
  sim <- simulate_bch(1000, c(0.7, 0.3), 1.5, c(0, 1), seed = 1)
  res <- bch_robust(sim$fit, sim$y)
  expect_equal(unname(rowSums(res$Classification_Matrix)), c(1, 1))
  expect_equal(unname(rowSums(res$Classification_Weights)), c(1, 1))
})

test_that("BCH recovers distal means with unequal, overlapping profiles", {
  # 80/20 profiles with substantial overlap (entropy ~ 0.55): the naive
  # modal-assignment means are strongly attenuated; BCH should not be.
  biases <- t(vapply(1:20, function(s) {
    sim <- simulate_bch(4000, c(0.8, 0.2), 1.5, c(0, 1), seed = 100 + s)
    unname(bch_robust(sim$fit, sim$y)$Profile_Means) - c(0, 1)
  }, numeric(2)))
  expect_lt(abs(mean(biases[, 1])), 0.03)
  expect_lt(abs(mean(biases[, 2])), 0.08)
})

test_that("BCH bootstrap correction runs and refits sequentially", {
  sim <- simulate_bch(400, c(0.6, 0.4), 2, c(0, 1), seed = 3)
  sim$fit$call_args$cores <- 4   # must be overridden inside the bootstrap refits
  set.seed(9)
  res <- bch_robust(sim$fit, sim$y, correction = "bootstrap", n_boot = 15)
  expect_false(is.null(res$Bootstrap_Correction))
  expect_true(all(res$Bootstrap_Correction$SE > 0))
  expect_lt(res$Bootstrap_Correction$Wald_p_value, 0.01)
})
