## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3.8,
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(RobustLPA)
set.seed(2026)

## -----------------------------------------------------------------------------
data(neuro_long)
head(neuro_long)
table(visits = table(neuro_long$ID))

## -----------------------------------------------------------------------------
fit <- robust_gmm(neuro_long, id = "ID", time = "Year",
                  outcomes = c("Memory", "Executive"), G = 3,
                  robust_method = "t", n_starts = 3)
fit

## -----------------------------------------------------------------------------
summary(fit)

## ----fig.alt = "Class mean trajectories over the individual trajectories"-----
plot_robust_gmm(fit)

## -----------------------------------------------------------------------------
head(sort(fit$weights))

## -----------------------------------------------------------------------------
classical <- estimate_gmm_robust(neuro_long, id = "ID", time = "Year",
                                 outcomes = c("Memory", "Executive"),
                                 n_classes = 2:4, robust = FALSE, n_starts = 2)
robust <- estimate_gmm_robust(neuro_long, id = "ID", time = "Year",
                              outcomes = c("Memory", "Executive"),
                              n_classes = 2:4, robust_method = "t", n_starts = 2)
classical$fit_table[, c("Model", "LogLik", "BIC", "Entropy", "Min_Size")]
robust$fit_table[, c("Model", "LogLik", "BIC", "Entropy", "Min_Size")]

## ----eval = FALSE-------------------------------------------------------------
# blrt_gmm_robust(neuro_long, id = "ID", time = "Year",
#                 outcomes = c("Memory", "Executive"), G = 3,
#                 robust_method = "t", n_samples = 200, cores = 4)

## -----------------------------------------------------------------------------
N <- length(unique(neuro_long$ID))
fit_l <- robust_gmm(neuro_long, id = "ID", time = "Year",
                    outcomes = c("Memory", "Executive", "Speed"), G = 3,
                    robust_method = "t", n_starts = 2,
                    lambda_growth = 9 / N, lambda_diff = 9 / N, group_diff = TRUE,
                    relax = TRUE)
summary(fit_l)

## -----------------------------------------------------------------------------
fit_u <- robust_gmm(neuro_long, id = "ID", time = "Year",
                    outcomes = c("Memory", "Executive", "Speed"), G = 3,
                    robust_method = "t", n_starts = 2)
rbind(unpenalized = fit_u$fit[, c("LogLik", "Parameters", "BIC")],
      lasso_relaxed = fit_l$fit[, c("LogLik", "Parameters", "BIC")])

## ----eval = FALSE-------------------------------------------------------------
# estimate_gmm_robust(neuro_long, id = "ID", time = "Year",
#                     outcomes = c("Memory", "Executive", "Speed"), n_classes = 3,
#                     tune_penalty = "bic", z_grid = c(1.5, 2, 2.5, 3, 4),
#                     robust_method = "t", group_diff = TRUE, relax = TRUE)

## -----------------------------------------------------------------------------
baseline <- neuro_long[!duplicated(neuro_long$ID), ]
baseline <- baseline[match(fit$ids, baseline$ID), ]
bch_biomarker <- bch_robust(fit, baseline$Biomarker)
round(bch_biomarker$Profile_Means)
bch_biomarker$ANOVA_Table

## -----------------------------------------------------------------------------
fit_b <- robust_gmm(neuro_long, id = "ID", time = "Year", outcomes = "Memory",
                    G = 3, robust_method = "t", engine = "MCMC",
                    mcmc_iter = 600, n_chains = 2, n_starts = 2)
fit_b

## ----fig.alt = "MCMC trace plots for the Memory slopes of the three classes"----
plot_mcmc_chains(fit_b, pars = c("beta[1,Memory:Year]", "beta[2,Memory:Year]",
                                 "beta[3,Memory:Year]", "nu"))

