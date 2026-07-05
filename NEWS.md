# RobustLPA 0.1.0

* Initial CRAN submission.
* Added robust Latent Profile Analysis (LPA) estimation using Huber weighting to handle multivariate outliers.
* Implemented Full Information Maximum Likelihood (FIML) via a high-performance C++ engine (`RcppArmadillo`) to natively handle missing data.
* Included multiple geometric variance-covariance model parameterizations (Models 1 to 6).
* Provided `blrt_robust()` for Bootstrapped Likelihood Ratio Tests.
* Provided `plot_robust_lpa()` for publication-ready visualizations using `ggplot2`.
