# RobustLPA 1.1.0

This release adds robust growth mixture models for longitudinal data and
is also a correctness release: every estimator in the package was
re-checked against a known reference (closed-form maximum likelihood,
`mclust`, `lme4`/`nlme`, or simulation with known truth), and the problems
found are fixed below. Several of them change numerical results relative
to 1.0.0; re-running analyses fitted with 1.0.0 is recommended.

## New: robust growth mixture models for longitudinal data

* **`robust_gmm()`** fits growth mixture models (finite mixtures of linear
  mixed-effects models; Verbeke & Lesaffre, 1996; Muthen & Shedden, 1999)
  and, with `random = "none"`, latent class growth analysis (Nagin, 1999),
  for one or several outcomes measured repeatedly on unbalanced occasions
  (long-format data; persons may have different visit times and missing
  outcomes). Each class has a polynomial mean trajectory per outcome
  (`degree`); within classes, persons have correlated random intercepts or
  intercepts and slopes (`random`, `re_structure`), with covariance and
  residual variances shared or class-specific (`re_cov`, `resid_var`).
  Every person contributes the exact likelihood of the values actually
  observed, so estimation is valid under missing-at-random drop-out.
  * Robust versions: multivariate-t classes (`robust_method = "t"`; the
    t linear mixed model of Pinheiro, Liu & Wu, 2001, with `nu` estimated),
    which down-weight whole outlying persons and keep a proper likelihood,
    or Huber weights.
  * EM engine: alternating ECM (Meng & van Dyk, 1997) with GLS updates of
    the trajectories and SQUAREM acceleration (Varadhan & Roland, 2008).
    The one-class model reproduces the maximum-likelihood linear mixed model
    of `lme4`/`nlme` exactly (checked in the tests).
  * MCMC engine: blocked Gibbs sampler (trajectories drawn with the random
    effects integrated out; persons allocated with random effects and latent
    scales integrated out), scaled inverse-Wishart prior on the random-effect
    covariances sampled by parameter expansion (Liu & Wu, 1999; Gelman et al.,
    2008), WAIC, Gelman-Rubin diagnostics and trace plots
    (`plot_mcmc_chains()`).
  * **LASSO for trajectories**, in both engines: `lambda_growth` shrinks the
    growth terms of every class toward zero (identifying classes that are
    stable on an outcome); `lambda_diff` shrinks the class trajectories toward
    each other, element-wise or by outcome (`group_diff = TRUE`, which selects
    the outcomes that differentiate the classes; Xie, Pan & Shen, 2008; the
    groups are standardized by the information of their coefficients, Simon &
    Tibshirani, 2012).
    Penalties are adaptive by default (Zou, 2006; Wang & Leng, 2008), so that
    `lambda = z^2 / N` sets to zero terms whose Wald statistic is below about
    `z`; the penalized step is solved exactly by ADMM, the BIC counts only the
    free coefficients, `relax = TRUE` refits the selected model without
    penalty (relaxed Lasso), and the MCMC engine uses the corresponding
    Bayesian (group, adaptive) Lasso priors.
* **`estimate_gmm_robust()`** fits a grid of numbers of classes (and
  random-effect specifications), optionally tuning the LASSO penalties by BIC
  or by K-fold cross-validation over persons.
* **`blrt_gmm_robust()`**: bootstrapped likelihood ratio test for the number
  of classes, simulating from the fitted null model on the observed visit
  schedule and missingness pattern.
* **`plot_robust_gmm()`** plots the class trajectories over the observed
  individual trajectories; `print()` and `summary()` methods.
* **`bch_robust()`** accepts `robust_gmm()` fits, to relate the trajectory
  classes to baseline characteristics or distal outcomes (one value per
  person, in the order of `$ids`); its bootstrap correction resamples whole
  persons.
* New example dataset **`neuro_long`** (400 simulated persons, up to six
  annual visits, three outcomes, drop-out and gross errors; its generating
  script is installed in `inst/scripts/`) and a second vignette,
  `vignette("robust-growth-mixture", package = "RobustLPA")`.

## Bug fixes

* **Parameter count, AIC, BIC and SABIC of `robust_lpa()`.** 1.0.0 left out
  of the count every mean with absolute value below 1e-5, even without a
  LASSO penalty. On standardized data with one profile all means are 0, so
  `p` parameters were missing (model 1 on the five `neuro_data` variables:
  5 instead of 10; model 6: 15 instead of 20) and the BIC of `G = 1` was
  about 27.6 too low, favouring one profile. Means are now left out only
  when the EM LASSO (`lambda > 0`) has set them exactly to zero.
* **Parallel workers on Windows.** With `cores > 1` the worker processes
  loaded RobustLPA from the default library, which could hold another
  version than the one used in the session (e.g. a new version installed in
  a separate library), silently mixing code from the two. They now load the
  session's copy, and a warning is given if they cannot.

* **BCH (`bch_robust()`) now estimates the classification error matrix
  correctly.** 1.0.0 computed `D[g, c]` as the average posterior of profile
  `g` among observations assigned to `c`, i.e. P(C = g | assigned c), the
  reverse of the P(assigned s | C = t) required by Bolck, Croon & Hagenaars
  (2004) / Vermunt (2010), and indexed the weight matrix `W = D^-1` with
  transposed subscripts. The two errors cancel only for symmetric
  classification error. With unequal, overlapping profiles (80/20, entropy
  about 0.55; true difference in the distal mean = 1) the corrected
  profile-2 mean had a median bias of +0.44 and RMSE of 1.9 in 1.0.0 (worse
  than the uncorrected comparison), versus -0.05 and 0.19 now. Rows of
  `$Classification_Matrix` now sum to 1.
* **Constrained models (1-5) could return an unconstrained solution.** The
  EM initialization did not apply the variance-covariance constraints, and
  the "keep the best iterate" logic could then return those unconstrained
  starting values. For example, `model = 1` on `neuro_data` returned
  non-diagonal, profile-specific covariances with a log-likelihood
  (-1395.4) *higher* than the model-1 maximum (-1469.5, as found by
  `mclust` EEI), which also invalidated BIC comparisons across models.
  Classical fits now reproduce the `mclust` maximum-likelihood solutions
  for models 1, 2, 3 and 6.
* **Missing data are now handled by exact maximum likelihood.** The E-step
  already used the observed-data (FIML) density, but the M-step used
  pairwise available-case moments, which are biased under MAR. The M-step
  now uses the exact EM for incomplete data (conditional expectations of
  the missing entries plus their conditional covariance; Ghahramani &
  Jordan, 1994). Example: with 30% MAR missingness the estimated mean of
  the incomplete variable was -0.37 (truth 0) in 1.0.0 and now matches the
  closed-form ML solution. The MCMC engine now uses the equivalent data
  augmentation step.
* **Numerical underflow for outlying observations.** Densities are now
  computed on the log scale throughout (log-sum-exp), so extreme
  observations receive valid posterior probabilities. In 1.0.0 an extreme
  outlier could get a posterior row of all zeros and a random profile.
* **MCMC models 1 and 5 now match the EM definitions.** 1.0.0's sampler
  used a single variance shared by all variables (a spherical model) for
  models 1 and 5; they now use variable-specific variances, as documented
  and as in the EM engine.
* **Models 4 and 5 are now fitted properly by both engines.** The EM engine
  used an approximate M-step that made the likelihood decrease and stopped
  at poor solutions; it now uses conditional-maximization steps that never
  decrease the likelihood. On `neuro_data` (classical fit, G = 2), model 4
  improves from log-likelihood -1358.7 (83% classification accuracy) to
  -1315.4 (94%), and model 5 from -1387.6 (52%) to -1318.8 (89%). The MCMC engine now samples
  the scales exactly (slice sampling) and the correlation matrices by
  Metropolis steps, instead of an approximate conditional that drifted
  away from the posterior.
* **Huber weights now accumulate across EM iterations.** In 1.0.0 each
  M-step recomputed the weights against a fresh *non-robust* estimate, so
  outliers could mask themselves; the weights are now computed against the
  current robust estimates (an iteratively reweighted M-estimator). Because
  the Huber fixed point is not a maximizer of the Gaussian likelihood, the
  final iterate is returned and no log-likelihood-decrease warnings are
  issued for Huber fits.
* `bch_robust(correction = "bootstrap")` refits now run with `cores = 1`
  inside each replicate (previously they inherited the original fit's
  `cores`, causing nested parallelism).
* Profile labels (MCMC draws, bootstrap refits) are now matched by an exact
  optimal assignment (Hungarian algorithm) on standardized means, instead
  of a greedy nearest-mean rule. MCMC draws are relabeled to a fixed pivot
  (the preliminary EM solution).
* `estimate_profiles_robust()` warns instead of silently ignoring
  `tune_lasso = TRUE` with the MCMC engine, validates its inputs, reports
  why a model failed, and uses the same parallel backend as the rest of the
  package.

## New features

* **Multivariate-t mixtures** (`robust_lpa(robust_method = "t")`, with `nu`
  estimated by default or fixed): a likelihood-based robust model (Peel &
  McLachlan, 2000) fitted by ECME in the EM engine and by an exact Gibbs
  sampler (latent scales; `nu` updated with the scales integrated out) in
  the MCMC engine. AIC/BIC, the BLRT and BCH are all used as intended with
  this model, and it is markedly more resistant than Huber weighting to
  gross outliers (in a test with 5% gross outliers on data with unit
  variances, the Huber variance estimates were 2.3-3.6 -- 3.5-6.1 in 1.0.0 --
  while the t fit estimated nu = 2.4 and a core scale of 0.7-0.8).
* `$weights`: every fit now reports each observation's robustness weight
  (small values flag outliers).
* The MCMC engine reports the **WAIC** (`$fit$WAIC`, `$waic`), starts its
  chains from a preliminary EM fit with dispersed perturbations, and adapts
  its Metropolis step sizes during burn-in.
* `blrt_robust()` simulates the null data from the fitted model family
  (multivariate t for `robust_method = "t"`).
* **Reproducibility:** `set.seed()` now gives identical results for any
  value of `cores` (one seed is drawn per work unit before dispatch).
* EM starts are initialized by k-means by default (`init = "kmeans"`;
  `init = "random"` restores random soft partitions).
* EM fits report `$converged` and `$iterations`.
* A `testthat` test suite checks every estimator against closed-form or
  reference solutions.

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
