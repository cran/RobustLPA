data(neuro_long)

test_that("the person-level E-step equals the direct multivariate density", {
  set.seed(1)
  N <- 5; K <- 2
  rows <- do.call(rbind, lapply(1:N, function(i) {
    tt <- sort(sample(0:5, sample(2:5, 1)))
    expand.grid(t = tt, k = 0:(K - 1), i = i)
  }))
  rows <- rows[order(rows$i, rows$k, rows$t), ]
  rows <- rows[-c(2, 7), ]
  y <- stats::rnorm(nrow(rows))
  starts <- c(0, cumsum(table(factor(rows$i, levels = 1:N))))
  beta <- c(0.5, -0.2, 1, 0.1)
  D <- matrix(c(1, .3, .2, .1, .3, .5, .1, .05, .2, .1, .8, .2, .1, .05, .2, .3), 4)
  sig2 <- c(0.4, 0.9)
  for (dist in 0:1) {
    es <- RobustLPA:::gmm_class_estep_cpp(y, rows$t, as.integer(rows$k), as.integer(starts), beta, D,
                                          sig2, K, 1L, 2L, dist, 4)
    for (i in 1:N) {
      idx <- which(rows$i == i)
      X <- t(vapply(idx, function(j) { v <- numeric(4); v[rows$k[j] * 2 + 1:2] <- c(1, rows$t[j]); v }, numeric(4)))
      V <- X %*% D %*% t(X) + diag(sig2[rows$k[idx] + 1], length(idx))
      r <- y[idx] - X %*% beta
      d <- c(t(r) %*% solve(V, r))
      n <- length(idx)
      ld <- as.numeric(determinant(V)$modulus)
      direct <- if (dist == 0) -0.5 * (n * log(2 * pi) + ld + d) else
        lgamma((4 + n) / 2) - lgamma(2) - n / 2 * log(4 * pi) - 0.5 * ld - (4 + n) / 2 * log1p(d / 4)
      expect_equal(es$logdens[i], direct, tolerance = 1e-10)
    }
  }
})

test_that("a one-class growth model is the maximum-likelihood linear mixed model", {
  skip_if_not_installed("lme4")
  fit <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 1, robust = FALSE, n_starts = 1, tol = 1e-12,
                    max_iter = 2000)
  ref <- lme4::lmer(Memory ~ Year + (Year | ID), data = neuro_long, REML = FALSE)
  expect_equal(fit$fit$LogLik, as.numeric(stats::logLik(ref)), tolerance = 1e-6)
  expect_equal(unname(fit$coefficients[[1]][1, ]), unname(lme4::fixef(ref)), tolerance = 1e-4)
  expect_equal(unname(fit$residual_var[1, 1]), stats::sigma(ref)^2, tolerance = 1e-4)
})

test_that("the multivariate one-class model equals a multivariate mixed model fitted by nlme", {
  skip_on_cran()
  skip_if_not_installed("nlme")
  d <- neuro_long
  dl <- rbind(data.frame(ID = d$ID, Year = d$Year, y = d$Memory, outc = "Memory"),
              data.frame(ID = d$ID, Year = d$Year, y = d$Executive, outc = "Executive"))
  dl <- dl[!is.na(dl$y), ]
  dl$outc <- factor(dl$outc, levels = c("Memory", "Executive"))
  ref <- nlme::lme(y ~ 0 + outc + outc:Year, random = list(ID = nlme::pdSymm(~ 0 + outc + outc:Year)),
                   weights = nlme::varIdent(form = ~ 1 | outc), data = dl, method = "ML",
                   control = nlme::lmeControl(maxIter = 500, msMaxIter = 500))
  fit <- robust_gmm(d, "ID", "Year", c("Memory", "Executive"), G = 1, robust = FALSE, n_starts = 1,
                    tol = 1e-12, max_iter = 3000)
  expect_equal(fit$fit$LogLik, as.numeric(stats::logLik(ref)), tolerance = 1e-6)
})

test_that("the t growth mixture recovers the simulated classes and resists gross errors", {
  skip_on_cran()
  set.seed(2)
  outs <- c("Memory", "Executive")
  ft <- robust_gmm(neuro_long, "ID", "Year", outs, G = 3, robust_method = "t", n_starts = 3)
  truth <- as.integer(neuro_long$True_Class[!duplicated(neuro_long$ID)])[match(ft$ids, unique(neuro_long$ID))]
  tab <- table(truth, ft$assignments)
  perm <- as.integer(RobustLPA:::solve_lsap_cpp(-unclass(tab))) + 1L
  expect_gt(sum(tab[cbind(1:3, perm)]) / length(truth), 0.85)
  slopes <- sort(vapply(ft$coefficients, function(m) m["Memory", "Year"], numeric(1)))
  expect_lt(max(abs(slopes - c(-5, -2, 0))), 0.3)
  expect_gt(ft$fit$Min_Size, 0.1)
  expect_true(ft$nu > 3 && ft$nu < 50)
})

test_that("adaptive LASSO penalties select outcomes and stable classes", {
  skip_on_cran()
  N <- length(unique(neuro_long$ID))
  set.seed(3)
  f <- robust_gmm(neuro_long, "ID", "Year", c("Memory", "Speed"), G = 3, robust_method = "t", n_starts = 2,
                  lambda_diff = 9 / N, group_diff = TRUE, lambda_growth = 9 / N, relax = TRUE)
  expect_identical(unname(f$penalty$outcome_selected), c(TRUE, FALSE))
  speed <- vapply(f$coefficients, function(m) m["Speed", ], numeric(2))
  expect_lt(max(abs(speed - speed[, 1])), 1e-6)            # Speed identical across classes
  mem_slopes <- vapply(f$coefficients, function(m) m["Memory", "Year"], numeric(1))
  expect_equal(sum(mem_slopes == 0), 1)                     # the stable class
  expect_true(f$penalty$relaxed)
  expect_lt(f$fit$Parameters, 2 + 3 * 4 + 10 + 2 + 1)
})

test_that("the growth-mixture MCMC engine agrees with EM", {
  skip_on_cran()
  set.seed(4)
  em <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, robust_method = "t", n_starts = 2)
  mc <- suppressWarnings(robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, robust_method = "t",
                                    engine = "MCMC", mcmc_iter = 600, n_chains = 2, n_starts = 2))
  # The two fits are independent, so their class labels are arbitrary: compare the
  # classes after ordering them by their Memory slope (label-invariant comparison).
  by_slope <- function(f) {
    m <- do.call(rbind, lapply(f$coefficients, function(b) b["Memory", ]))
    unname(m[order(m[, "Year"]), , drop = FALSE])
  }
  expect_equal(by_slope(mc), by_slope(em), tolerance = 0.1)
  expect_equal(sort(unname(mc$proportions)), sort(unname(em$proportions)), tolerance = 0.1)
  expect_true(is.finite(mc$fit$WAIC))
  expect_lt(max(mc$mcmc_diagnostics$Rhat, na.rm = TRUE), 1.2)
  expect_s3_class(plot_mcmc_chains(mc, pars = c("beta[1,Memory:Year]", "nu")), "ggplot")
})

test_that("simulation, new-data likelihood, BCH and reproducibility work on growth-mixture fits", {
  set.seed(5)
  f <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, n_starts = 2)
  sim <- RobustLPA:::.gmm_simulate(f)
  expect_identical(is.na(sim$Memory), is.na(f$data$Memory))
  expect_equal(RobustLPA:::.gmm_loglik_newdata(f, neuro_long), f$fit$LogLik, tolerance = 1e-8)
  base <- neuro_long[!duplicated(neuro_long$ID), ]
  b <- bch_robust(f, base$Biomarker[match(f$ids, base$ID)])
  expect_equal(unname(rowSums(b$Classification_Matrix)), c(1, 1))
  expect_error(bch_robust(f, base$Biomarker[1:10]), "persons")
  expect_s3_class(plot_robust_gmm(f), "ggplot")
  expect_output(print(summary(f)), "Class trajectories")
  skip_on_cran()
  skip_on_os("windows")
  set.seed(6)
  a1 <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, n_starts = 3, cores = 1, init = "random")
  set.seed(6)
  a2 <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, n_starts = 3, cores = 2, init = "random")
  expect_identical(a1$coefficients, a2$coefficients)
})

test_that("robust_gmm validates its inputs", {
  expect_error(robust_gmm(neuro_long, "ID", "Year", "Nope", G = 2), "not found")
  expect_error(robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, degree = 0), "degree >= 1")
  expect_error(robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, lambda_growth = -1), "non-negative")
  expect_error(robust_gmm(as.matrix(neuro_long), "ID", "Year", "Memory", G = 2), "data.frame")
})

test_that("one outcome with a constant trajectory (one coefficient per class) works", {
  set.seed(7)
  for (rnd in c("intercept", "none")) {
    f <- robust_gmm(neuro_long, "ID", "Year", "Memory", G = 2, degree = 0, random = rnd, n_starts = 2)
    expect_equal(dim(f$coefficients[[1]]), c(1L, 1L))
    expect_equal(RobustLPA:::.gmm_loglik_newdata(f, neuro_long), f$fit$LogLik, tolerance = 1e-8)
    expect_identical(is.na(RobustLPA:::.gmm_simulate(f)$Memory), is.na(f$data$Memory))
    expect_output(print(summary(f)), "Class trajectories")
  }
  expect_equal(dim(f$coefficients[[1]]), c(1L, 1L))
  set.seed(8)
  b <- RobustLPA:::.gmm_refit_resample(f, sample(length(f$ids), replace = TRUE))
  expect_true(is.list(b))
})

test_that("random starts draw distinct class centres", {
  set.seed(9)
  Fs <- rbind(matrix(stats::rnorm(40, -3), 20), matrix(stats::rnorm(40, 3), 20))
  cl <- RobustLPA:::.gmm_random_centers(Fs, 3)
  expect_setequal(unique(cl), 1:3)
  expect_null(RobustLPA:::.gmm_random_centers(Fs[1:2, ], 3))
})
