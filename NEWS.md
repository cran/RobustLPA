# RobustLPA 1.0.0

Major release: a new Bayesian MCMC estimation engine alongside the existing
EM engine, an optional classical (non-robust) estimation mode, LASSO
regularization, a new BCH auxiliary-variable method, parallel computing
throughout, S3 print/summary methods, a new bundled dataset, and a new
package vignette.

* **New MCMC engine.** `robust_lpa(engine = "MCMC")` estimates the same six
  variance-covariance parameterizations via Gibbs sampling, with a Bayesian
  Lasso (Laplace-prior) option for the profile means (`prior_laplace`) and
  the same Huber-weighted robust/classical toggle as the EM engine
  (`robust`/`alpha`). Runs multiple independent chains by default
  (`n_chains = 4`), aligns each chain's arbitrarily-ordered profile labels
  to a common ordering before pooling, and reports classic Gelman-Rubin
  \eqn{\hat{R}} and effective sample size convergence diagnostics
  (`$mcmc_diagnostics`, via the new `coda` dependency) for every scalar
  parameter, warning if any exceeds the standard 1.1 R-hat threshold.
  `plot_mcmc_chains()` draws multi-chain trace plots via `bayesplot`.
* **Optional robustness.** `robust_lpa()` gains a `robust` argument (default
  `TRUE`) and `alpha`, so either engine can be run with Huber-weighted
  robust estimation (as in 0.1.0) or with classical (non-robust) maximum
  likelihood, using the same engine, missing-data handling, and variance-
  covariance parameterizations either way.
* **LASSO regularization.** The EM engine supports LASSO-type soft-
  thresholding shrinkage of the profile means (`lambda`), and
  `estimate_profiles_robust(tune_lasso = TRUE)` selects `lambda` by k-fold
  cross-validation (`k_folds`, `lambda_grid`).
* **Parallel computing.** A `cores` argument runs independent work in
  parallel throughout the package: EM random restarts or MCMC chains within
  a single `robust_lpa()` call; the model/profile grid (and cross-validation
  folds) in `estimate_profiles_robust()`; bootstrap replicates in
  `blrt_robust()`; and bootstrap correction replicates in `bch_robust()`.
  Uses `parallel::mclapply()` on macOS/Linux and a `parallel::makeCluster()`
  PSOCK backend on Windows.
* **BCH auxiliary-variable analysis.** New `bch_robust()` implements the
  Bolck-Croon-Hagenaars (2004) three-step method for relating fitted
  profiles to a continuous auxiliary/distal outcome variable, correcting for
  classification error. Optionally (`correction = "bootstrap"`) adds a
  nonparametric bootstrap correction for classification uncertainty in the
  step-1 model, reporting bootstrap standard errors, confidence intervals,
  and a Wald chi-square test for the profile means (`$Bootstrap_Correction`)
  -- a practical approximation to the Bakk, Oberski & Vermunt (2014)
  sandwich correction, preferable to the base `$ANOVA_Table` F-test (which
  treats the classification matrix as fixed/known) for publication-grade
  inference.
* **Enhanced bootstrapped likelihood ratio test.** `blrt_robust()` (present
  since 0.1.0) now supports both estimation engines (`engine = "EM"` or
  `"MCMC"`), replicates the observed data's FIML missingness pattern in
  every simulated bootstrap sample, and runs its bootstrap replicates in
  parallel via `cores`.
* **S3 classes for fitted models.** `robust_lpa()` now returns an object of
  class `"robust_lpa"` with `print()` and `summary()` methods: `print()` is
  a compact overview (engine, model, profile count, N, headline fit
  indices, mixing proportions), and `summary()` adds per-profile means,
  profile sizes, and, for the MCMC engine, the Gelman-Rubin/effective
  sample size range -- both focused on what's needed for a first read,
  rather than dumping the full fitted object.
* **Refit-ready fitted objects.** `robust_lpa()` fits now store `$data` (the
  numeric matrix used for estimation) and `$call_args` (every argument
  controlling the fit), so other functions can refit the identical
  specification on new or resampled data; used internally by
  `bch_robust()`'s bootstrap correction.
* **Regularization safety net for covariance estimation.** Variance-
  covariance models that allow off-diagonal terms (models 3-6) use a
  scale-relative eigenvalue floor to keep every estimated covariance matrix
  numerically positive-definite, without perturbing already well-behaved
  matrices. Fits (from either engine) whose smallest profile holds under 2%
  of observations emit a transparency warning identifying the model,
  profile count, and engine, since this pattern can indicate an unusually
  small genuine subgroup or an unconstrained-covariance model over-fitting
  a handful of near-coincident points -- worth a second look either way
  (e.g. a more constrained `model`, more `n_starts`, or fewer profiles).
* **New bundled dataset.** `neuro_data` (n = 250, two groups) is a new
  synthetic dataset with a genuine group-level difference in covariance
  structure (not only in means) and calibrated, variable-magnitude outlier
  contamination, so that model/profile selection, robust-vs-classical
  comparisons, and classification recovery all behave as a realistic
  worked example throughout the package's documentation. See `?neuro_data`
  for details and `data-raw/generate_neuro_data.R` for the full generative
  code.
* **Package vignette.** A new vignette (`vignette("RobustLPA")`) walks
  through the full workflow: choosing a variance-covariance model, EM vs.
  MCMC estimation, robust vs. classical estimation, LASSO regularization,
  model/profile selection, the bootstrapped likelihood ratio test, and BCH
  auxiliary-variable analysis, on `neuro_data`.
* Packaging: declared the `parallel` dependency under `Suggests`; reference
  DOIs throughout the documentation use the `\doi{}` Rd macro, as recommended
  by CRAN.

# RobustLPA 0.1.0

* Initial CRAN submission.
* Added robust Latent Profile Analysis (LPA) estimation using Huber
  weighting to handle multivariate outliers.
* Implemented Full Information Maximum Likelihood (FIML) via a
  high-performance C++ engine (`RcppArmadillo`) to natively handle missing
  data.
* Included multiple geometric variance-covariance model parameterizations
  (Models 1 to 6).
* Provided `blrt_robust()` for Bootstrapped Likelihood Ratio Tests.
* Provided `plot_robust_lpa()` for publication-ready visualizations using
  `ggplot2`.
