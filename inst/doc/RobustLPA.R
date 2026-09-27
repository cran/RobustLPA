## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4.2,
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(RobustLPA)
set.seed(2026)  # every result below is reproducible, for any number of cores

## -----------------------------------------------------------------------------
data(neuro_data)
str(neuro_data)
table(neuro_data$True_Profile)

## -----------------------------------------------------------------------------
vars <- c("Memory", "Attention", "Executive_Functions", "RT_Stroop", "RT_TMT")
x <- scale(as.matrix(neuro_data[, vars]))
head(x)

## -----------------------------------------------------------------------------
fit_em <- robust_lpa(x, G = 2, model = 6, n_starts = 5)
fit_em

## -----------------------------------------------------------------------------
summary(fit_em)

## -----------------------------------------------------------------------------
table(True_Profile = neuro_data$True_Profile, Assigned = fit_em$assignments)

## -----------------------------------------------------------------------------
fit_t <- robust_lpa(x, G = 2, model = 6, n_starts = 5, robust_method = "t")
fit_classical <- robust_lpa(x, G = 2, model = 6, n_starts = 5, robust = FALSE)
rbind(
  huber     = sort(sapply(fit_em$means, `[`, "RT_Stroop")),
  t         = sort(sapply(fit_t$means, `[`, "RT_Stroop")),
  classical = sort(sapply(fit_classical$means, `[`, "RT_Stroop"))
)
fit_t$nu

## -----------------------------------------------------------------------------
set.seed(6)
contaminated <- matrix(rnorm(400 * 3), 400, 3)
idx <- sample(400, 20)
contaminated[idx, ] <- contaminated[idx, ] + matrix(rnorm(60, 0, 15), 20, 3)
sapply(list(
  classical = robust_lpa(contaminated, G = 1, model = 6, robust = FALSE),
  huber     = robust_lpa(contaminated, G = 1, model = 6),
  t         = robust_lpa(contaminated, G = 1, model = 6, robust_method = "t")
), function(f) round(diag(f$covariances[[1]]), 2))

## -----------------------------------------------------------------------------
head(order(fit_t$weights))
round(head(sort(fit_t$weights)), 3)

## -----------------------------------------------------------------------------
fit_lasso <- robust_lpa(x, G = 2, model = 6, n_starts = 3, lambda = 0.15)
summary(fit_lasso)

## -----------------------------------------------------------------------------
x_na <- x
set.seed(1)
na_idx <- cbind(
  sample(nrow(x_na), 15),
  sample(ncol(x_na), 15, replace = TRUE)
)
x_na[na_idx] <- NA
mean(is.na(x_na))

fit_fiml <- robust_lpa(x_na, G = 2, model = 6, n_starts = 5)
summary(fit_fiml)

## -----------------------------------------------------------------------------
grid <- estimate_profiles_robust(x, n_profiles = 1:3, models = 1:6, n_starts = 5)
grid$fit_table[order(grid$fit_table$BIC), ]

## -----------------------------------------------------------------------------
summary(grid$models[["model_6_profiles_2"]])

## ----fig.alt = "Profile plot of the best-fitting model"-----------------------
plot_robust_lpa(grid, title = "Best-fitting model (lowest BIC)")

## -----------------------------------------------------------------------------
grid_lasso <- estimate_profiles_robust(
  x, n_profiles = 2, models = 6, n_starts = 3,
  tune_lasso = TRUE, k_folds = 5, lambda_grid = c(0, 0.05, 0.1, 0.2)
)
grid_lasso$fit_table[, c("Model", "Profiles", "BIC", "Lambda")]

## -----------------------------------------------------------------------------
blrt_res <- blrt_robust(x, G = 2, model = 6, n_samples = 20, n_starts = 3)
blrt_res

## -----------------------------------------------------------------------------
fit_mcmc <- robust_lpa(x, G = 2, model = 6, engine = "MCMC", robust_method = "t",
                       mcmc_iter = 500, n_chains = 4, prior_laplace = 0.1)
summary(fit_mcmc)

## ----fig.alt = "MCMC trace plots for two profile means, one mixing proportion and the t degrees of freedom"----
plot_mcmc_chains(fit_mcmc, pars = c("mu[1,1]", "mu[2,1]", "pi[1]", "nu"))

## -----------------------------------------------------------------------------
x_reduced <- scale(as.matrix(neuro_data[, c("Memory", "Attention",
                                             "Executive_Functions", "RT_Stroop")]))
fit_reduced <- robust_lpa(x_reduced, G = 2, model = 6, n_starts = 5)

bch_res <- bch_robust(fit_reduced, neuro_data$RT_TMT)
bch_res$Profile_Means
bch_res$ANOVA_Table

## -----------------------------------------------------------------------------
bch_boot <- bch_robust(fit_reduced, neuro_data$RT_TMT,
                        correction = "bootstrap", n_boot = 30)
bch_boot$Bootstrap_Correction

## ----eval = FALSE-------------------------------------------------------------
# grid_parallel <- estimate_profiles_robust(x, n_profiles = 1:3, models = 1:6,
#                                            n_starts = 5, cores = 4)
# fit_mcmc_parallel <- robust_lpa(x, G = 2, model = 6, engine = "MCMC",
#                                  mcmc_iter = 2000, n_chains = 4, cores = 4)

