#' Robust Growth Mixture Model (Latent Class Growth Analysis)
#'
#' Fits a growth mixture model (GMM) -- a finite mixture of linear
#' mixed-effects models for longitudinal data (Verbeke & Lesaffre, 1996;
#' Muthen & Shedden, 1999) -- or, with \code{random = "none"}, a latent class
#' growth analysis (LCGA; Nagin, 1999), with one or several outcomes measured
#' repeatedly on the same persons. Each latent class ("profile of change")
#' has its own polynomial mean trajectory for every outcome; within a class,
#' persons deviate from it through correlated random effects and residual
#' errors. The model can be estimated by maximum likelihood (EM engine) or
#' by Gibbs sampling (MCMC engine), in a classical (Gaussian), Huber-weighted
#' or multivariate-t (robust) version, optionally with LASSO penalties on
#' the class trajectories.
#'
#' @section Model:
#' For person \eqn{i} in class \eqn{g}, the stacked vector \eqn{y_i} of all
#' his/her observed values (all outcomes, all occasions) is
#' \deqn{y_i = X_i \beta_g + Z_i b_i + e_i, \quad b_i \sim N(0, D_g / u_i),
#'   \quad e_i \sim N(0, R_{g} / u_i),}
#' where \eqn{X_i} contains, for each observation of outcome \eqn{k} at time
#' \eqn{t}, the polynomial basis \eqn{(1, t, \dots, t^{degree})} in the
#' block of outcome \eqn{k}; \eqn{Z_i} the first 1 (\code{random =
#' "intercept"}) or 2 (\code{random = "slope"}) of those terms; \eqn{R_g}
#' is diagonal with one residual variance per outcome; and \eqn{u_i = 1}
#' (Gaussian) or \eqn{u_i \sim \mathrm{Gamma}(\nu/2, \nu/2)} (multivariate
#' t: Pinheiro, Liu & Wu, 2001). Marginally, \eqn{y_i} is Gaussian or
#' multivariate t with mean \eqn{X_i \beta_g} and scale
#' \eqn{V_{ig} = Z_i D_g Z_i' + R_g}.
#'
#' Persons may be observed at different times, on different numbers of
#' occasions, and not on every outcome at every occasion: each person
#' contributes the exact likelihood of the values actually observed, so
#' the estimates are maximum likelihood under missing-at-random
#' missingness (including drop-out that depends on earlier observed
#' values). Every outcome is standardized internally (all reported
#' estimates are on the original scale); the time variable is used as
#' supplied, so center it where the intercept should be interpreted (e.g.
#' years since baseline).
#'
#' @section Robust estimation:
#' \describe{
#'   \item{\code{robust_method = "t"}}{A mixture of multivariate-t linear
#'     mixed models: the random effects and the errors of a person share
#'     the latent scale \eqn{u_i}, so a person whose trajectory is far from
#'     every class (or who has a few gross errors) is down-weighted as a
#'     whole through \eqn{E[u_i] = (\nu + n_i) / (\nu + d_i)}, with
#'     \eqn{d_i} the Mahalanobis distance of \eqn{y_i} and \eqn{n_i} its
#'     length. It is a proper likelihood model, fitted by ECM with an ECME
#'     update of \eqn{\nu}, so information criteria, the bootstrapped
#'     likelihood ratio test (\code{\link{blrt_gmm_robust}}) and the BCH
#'     method are used as intended. Recommended.}
#'   \item{\code{robust_method = "huber"}}{Huber weights computed from
#'     \eqn{d_i} with a \eqn{\chi^2_{n_i}} cutoff (\code{alpha}) down-weight
#'     outlying persons in the M-step; an estimating-equation approach
#'     whose reported log-likelihood is the Gaussian one evaluated at the
#'     robust estimates (a heuristic for information criteria).}
#' }
#'
#' @section LASSO penalties:
#' Two penalties on the class trajectories (fixed effects) are available,
#' separately or together; both act on the standardized outcome scale and
#' are expressed per person, i.e. the EM engine maximizes
#' \deqn{\ell(\theta)/N - \lambda_{growth} \sum_{g,k,l \ge 1} |\beta_{gkl}|
#'   - \lambda_{diff} \sum_{g,k,l} |\beta_{gkl} - \bar\beta_{kl}|,}
#' with \eqn{\bar\beta_{kl}} the (unweighted) mean over classes.
#' \describe{
#'   \item{\code{lambda_growth}}{Sparse trajectories: the growth terms
#'     (slope, quadratic, ...) of every class and outcome are shrunk toward
#'     0, intercepts are not penalized. A coefficient set exactly to 0 means
#'     that the class does not change on that outcome (e.g. a "stable"
#'     class). This is the Lasso of Tibshirani (1996) applied to the
#'     fixed effects of a mixture of regressions (Khalili & Chen, 2007;
#'     Du et al., 2013).}
#'   \item{\code{lambda_diff}}{Sparse class differences: every coefficient
#'     of every class is shrunk toward the across-class mean, so that
#'     coefficients (or, with \code{group_diff = TRUE}, whole outcomes) on
#'     which the classes do not differ are fused. With \code{group_diff =
#'     TRUE} the penalty is \eqn{\lambda_{diff} \sum_k w_k \| S_k (B_k -
#'     1\bar\beta_k') \|}, where \eqn{S_k} scales every coefficient by the
#'     square root of its information (so that intercepts and slopes are
#'     penalized on comparable scales; Simon & Tibshirani, 2012) and
#'     \eqn{w_k} is the adaptive weight; it removes outcome \eqn{k} from the
#'     class separation altogether when its group is set to zero -- the
#'     grouped variable-selection penalty of Xie, Pan & Shen (2008), which
#'     extends the penalized model-based clustering of Pan & Shen (2007)
#'     used by \code{\link{robust_lpa}}'s \code{lambda}.}
#' }
#' By default (\code{adaptive = TRUE}) the penalties are adaptive (Zou,
#' 2006; Wang & Leng, 2008 for the group version): every coefficient,
#' deviation or outcome group is weighted by the reciprocal of its
#' unpenalized estimate, which gives consistent selection and makes the
#' penalty level interpretable -- a term is set to zero approximately when
#' its Wald statistic is below \eqn{\sqrt{N\lambda}} in absolute value, so
#' \code{lambda = z^2 / N} corresponds to a threshold \code{z} (e.g.
#' \code{4 / N} for \eqn{|z| < 2}). The penalized fit is started from the
#' unpenalized maximum-likelihood solution (which also supplies the
#' weights). The penalized fixed-effects step of each iteration is solved
#' exactly (to numerical tolerance) by ADMM (Boyd et al., 2011); the number of
#' fixed-effect parameters used by AIC/BIC is the number of free
#' coefficients left by the zeros and fusions (the generalized-Lasso degrees
#' of freedom; Tibshirani & Taylor, 2011). Penalized estimates are shrunk
#' toward zero / toward each other: \code{relax = TRUE} refits the model
#' without penalty while keeping the selected zeros and fusions (the relaxed
#' Lasso; Meinshausen, 2007), which is preferable for reporting and
#' inference. Choose \code{lambda_growth} / \code{lambda_diff} with
#' \code{\link{estimate_gmm_robust}} (BIC or cross-validation over persons).
#' In the MCMC engine the same penalties become Bayesian-Lasso (Laplace)
#' priors with rate \eqn{N\lambda} times the adaptive weight of each term
#' (Park & Casella, 2008; Bayesian group Lasso of Kyung et al., 2010, when
#' \code{group_diff = TRUE}), whose posterior mode is the EM penalized
#' estimate; \code{relax} does not apply to the MCMC engine.
#'
#' @section Estimation:
#' The EM engine uses an alternating ECM algorithm (Meng & van Dyk, 1997):
#' each iteration updates the class proportions and the class trajectories
#' by (penalized) weighted generalized least squares with the random effects
#' integrated out, then, after a fresh E-step, the residual variances, the
#' random-effect covariances and (t model) \eqn{\nu}; every step increases
#' the (penalized) observed-data log-likelihood. The iterations are
#' accelerated by SQUAREM (Varadhan & Roland, 2008) with a monotonicity
#' safeguard. Several starts are run (k-means on person-level
#' least-squares trajectories by default) and the best is kept.
#' @section MCMC details:
#' The Gibbs sampler allocates persons with the random effects and latent
#' scales integrated out, then draws \eqn{\nu} (t model; random-walk
#' Metropolis with the latent scales integrated out, Gamma(2, 0.1) prior),
#' the latent scales, the class trajectories with the random effects
#' integrated out, the random effects, the random-effect covariance
#' matrices, the residual variances (inverse-gamma(1, 0.1) on the
#' standardized scale) and the mixing proportions (Dirichlet(1, ..., 1)).
#' The random-effect covariances have a scaled inverse-Wishart prior
#' (O'Malley & Zaslavsky, 2008), \eqn{D = \mathrm{diag}(a) \Psi
#' \mathrm{diag}(a)} with \eqn{\Psi \sim IW(q + 1, I)} per block and
#' \eqn{a_j \sim N(0, A_j^2)} (\eqn{A_j} = 5 standardized units for
#' intercepts, divided by the standard deviation of the time term for
#' slopes), sampled by parameter expansion (Liu & Wu, 1999; Gelman et al.,
#' 2008), which avoids the slow mixing of small variance components. The
#' fixed effects have vague normal priors, combined with the Bayesian-Lasso
#' priors when penalties are requested. The chains start from the EM fit
#' with perturbed trajectories, draws are relabeled to that EM solution by
#' an optimal assignment of the class mean trajectories, and the first half
#' of each chain is discarded. The Huber option is a heuristic in the MCMC
#' engine, as in \code{\link{robust_lpa}}.
#'
#' @param data A data.frame in long format: one row per person and
#'   occasion, with the person identifier, the time variable and the
#'   outcome columns (\code{NA} allowed in the outcomes).
#' @param id,time Names of the person-identifier and time columns.
#' @param outcomes Character vector with the names of one or more outcome
#'   columns (modelled jointly).
#' @param G Number of latent classes.
#' @param degree Degree of the polynomial trajectory in \code{time}
#'   (default \code{1}, linear).
#' @param random Random effects within classes: \code{"slope"} (random
#'   intercept and slope for every outcome; default), \code{"intercept"}
#'   (random intercept only) or \code{"none"} (no random effects: latent
#'   class growth analysis).
#' @param re_structure Covariance structure of the random effects:
#'   \code{"full"} (all correlated, default), \code{"block"} (correlated
#'   within an outcome, independent across outcomes) or \code{"diagonal"}.
#' @param re_cov Random-effect covariance \code{"equal"} across classes
#'   (default, the usual and more stable choice) or \code{"varying"}.
#' @param resid_var Residual variances \code{"equal"} across classes
#'   (default) or \code{"varying"}.
#' @param engine \code{"EM"} (default) or \code{"MCMC"}.
#' @param robust Logical; \code{FALSE} fits the classical Gaussian model.
#' @param robust_method \code{"huber"} (default) or \code{"t"}; see
#'   "Robust estimation".
#' @param nu Degrees of freedom of the t model: \code{NULL} (default)
#'   estimates it, a number fixes it.
#' @param alpha Huber significance level.
#' @param lambda_growth,lambda_diff Non-negative LASSO penalty levels (see
#'   "LASSO penalties"); default \code{0} (no penalty).
#' @param group_diff Logical; group-wise (by outcome) difference penalty.
#' @param adaptive Logical; adaptive-Lasso weights from the unpenalized fit
#'   (default \code{TRUE}; see "LASSO penalties").
#' @param relax Logical; refit without penalty keeping the selected zeros
#'   and fusions (EM engine only).
#' @param n_starts Number of EM starts (also used for the EM fit that
#'   initializes the MCMC engine).
#' @param max_iter Maximum number of EM iterations.
#' @param tol Relative convergence tolerance on the (penalized)
#'   log-likelihood.
#' @param init Initialization of the EM starts: \code{"kmeans"} (k-means on
#'   person-level least-squares trajectories; default) or \code{"random"}
#'   (\code{G} persons drawn at random as class centres, every person assigned
#'   to the nearest one).
#' @param mcmc_iter Iterations per MCMC chain (first half discarded).
#' @param n_chains Number of MCMC chains.
#' @param cores Number of cores for the EM starts / MCMC chains. Results are
#'   identical for any value given the same \code{set.seed()}.
#' @return An object of class \code{"robust_gmm"}, a list with:
#'   \describe{
#'     \item{engine, robust_method, nu}{Estimation settings (\code{nu}: the
#'       estimated or fixed t degrees of freedom, \code{NULL} otherwise).}
#'     \item{coefficients}{A list of \code{G} matrices (outcomes x
#'       polynomial terms): the class mean trajectories.}
#'     \item{random_cov}{A list of \code{G} random-effect covariance matrices
#'       (\code{NULL} for \code{random = "none"}).}
#'     \item{residual_var}{A \code{G x K} matrix of residual variances.}
#'     \item{proportions}{Class proportions.}
#'     \item{probabilities, assignments}{Posterior class probabilities
#'       (persons x classes) and modal classes, in the order of \code{ids}.}
#'     \item{ids}{The person identifiers, in the order used by every
#'       person-level output (use it to align auxiliary variables for
#'       \code{\link{bch_robust}}).}
#'     \item{weights}{Per-person robustness weights, averaged over classes
#'       with the posterior probabilities. With \code{robust_method = "huber"}
#'       they lie in (0, 1] (1 = full weight). With \code{robust_method = "t"}
#'       they are the E-step weights \eqn{(\nu + n_i)/(\nu + \delta_i)}, where
#'       \eqn{n_i} is the number of observed values and \eqn{\delta_i} the
#'       squared Mahalanobis distance of person \eqn{i}: they exceed 1 for
#'       persons closer to their class trajectory than expected and are small
#'       for outlying persons. All 1 when \code{robust = FALSE}.}
#'     \item{random_effects}{Posterior means of each person's random
#'       effects under his/her modal class.}
#'     \item{fit}{One-row data.frame: \code{Classes}, \code{LogLik},
#'       \code{Parameters}, \code{AIC}, \code{BIC} (with \eqn{\log N},
#'       \eqn{N} = persons), \code{SABIC}, \code{Entropy}, \code{Min_Size},
#'       \code{Max_Size}, the penalty levels, and \code{WAIC} (MCMC).}
#'     \item{penalty}{Penalty settings and the selected zeros / fusions.}
#'     \item{converged, iterations}{EM convergence information.}
#'     \item{mcmc_draws, mcmc_diagnostics, waic}{MCMC output (see
#'       \code{\link{plot_mcmc_chains}}).}
#'     \item{spec, internal, call_args, data}{Model specification,
#'       standardized-scale estimates, arguments and data used for refitting.}
#'   }
#' @references
#'   Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011).
#'   Distributed optimization and statistical learning via the alternating
#'   direction method of multipliers. \emph{Foundations and Trends in Machine
#'   Learning}, 3(1), 1-122. \doi{10.1561/2200000016}
#'
#'   Du, Y., Khalili, A., Neslehova, J. G., & Steele, R. J. (2013).
#'   Simultaneous fixed and random effects selection in finite mixture of
#'   linear mixed-effects models. \emph{Canadian Journal of Statistics},
#'   41(4), 596-616. \doi{10.1002/cjs.11192}
#'
#'   Gelman, A., van Dyk, D. A., Huang, Z., & Boscardin, W. J. (2008). Using
#'   redundant parameterizations to fit hierarchical models. \emph{Journal of
#'   Computational and Graphical Statistics}, 17(1), 95-122.
#'   \doi{10.1198/106186008X287337}
#'
#'   Khalili, A., & Chen, J. (2007). Variable selection in finite mixture of
#'   regression models. \emph{Journal of the American Statistical
#'   Association}, 102(479), 1025-1038. \doi{10.1198/016214507000000590}
#'
#'   Kyung, M., Gill, J., Ghosh, M., & Casella, G. (2010). Penalized
#'   regression, standard errors, and Bayesian lassos. \emph{Bayesian
#'   Analysis}, 5(2), 369-411. \doi{10.1214/10-BA607}
#'
#'   Liu, J. S., & Wu, Y. N. (1999). Parameter expansion for data
#'   augmentation. \emph{Journal of the American Statistical Association},
#'   94(448), 1264-1274. \doi{10.1080/01621459.1999.10473879}
#'
#'   Meinshausen, N. (2007). Relaxed Lasso. \emph{Computational Statistics &
#'   Data Analysis}, 52(1), 374-393. \doi{10.1016/j.csda.2006.12.019}
#'
#'   Meng, X.-L., & van Dyk, D. (1997). The EM algorithm -- an old folk-song
#'   sung to a fast new tune. \emph{Journal of the Royal Statistical Society:
#'   Series B}, 59(3), 511-567. \doi{10.1111/1467-9868.00082}
#'
#'   Muthen, B., & Shedden, K. (1999). Finite mixture modeling with mixture
#'   outcomes using the EM algorithm. \emph{Biometrics}, 55(2), 463-469.
#'   \doi{10.1111/j.0006-341X.1999.00463.x}
#'
#'   Nagin, D. S. (1999). Analyzing developmental trajectories: A
#'   semiparametric, group-based approach. \emph{Psychological Methods},
#'   4(2), 139-157. \doi{10.1037/1082-989X.4.2.139}
#'
#'   O'Malley, A. J., & Zaslavsky, A. M. (2008). Domain-level covariance
#'   analysis for multilevel survey data with structured nonresponse.
#'   \emph{Journal of the American Statistical Association}, 103(484),
#'   1405-1418. \doi{10.1198/016214508000000724}
#'
#'   Pan, W., & Shen, X. (2007). Penalized model-based clustering with
#'   application to variable selection. \emph{Journal of Machine Learning
#'   Research}, 8, 1145-1164.
#'
#'   Park, T., & Casella, G. (2008). The Bayesian Lasso. \emph{Journal of the
#'   American Statistical Association}, 103(482), 681-686.
#'   \doi{10.1198/016214508000000337}
#'
#'   Pinheiro, J. C., Liu, C., & Wu, Y. N. (2001). Efficient algorithms for
#'   robust estimation in linear mixed-effects models using the multivariate
#'   t distribution. \emph{Journal of Computational and Graphical
#'   Statistics}, 10(2), 249-276. \doi{10.1198/10618600152628059}
#'
#'   Simon, N., & Tibshirani, R. (2012). Standardization and the group lasso
#'   penalty. \emph{Statistica Sinica}, 22(3), 983-1001.
#'
#'   Tibshirani, R. (1996). Regression shrinkage and selection via the lasso.
#'   \emph{Journal of the Royal Statistical Society: Series B}, 58(1),
#'   267-288. \doi{10.1111/j.2517-6161.1996.tb02080.x}
#'
#'   Tibshirani, R. J., & Taylor, J. (2011). The solution path of the
#'   generalized lasso. \emph{The Annals of Statistics}, 39(3), 1335-1371.
#'   \doi{10.1214/11-AOS878}
#'
#'   Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods
#'   for accelerating the convergence of any EM algorithm. \emph{Scandinavian
#'   Journal of Statistics}, 35(2), 335-353.
#'   \doi{10.1111/j.1467-9469.2007.00585.x}
#'
#'   Verbeke, G., & Lesaffre, E. (1996). A linear mixed-effects model with
#'   heterogeneity in the random-effects population. \emph{Journal of the
#'   American Statistical Association}, 91(433), 217-221.
#'   \doi{10.1080/01621459.1996.10476679}
#'
#'   Wang, H., & Leng, C. (2008). A note on adaptive group lasso.
#'   \emph{Computational Statistics & Data Analysis}, 52(12), 5277-5286.
#'   \doi{10.1016/j.csda.2008.05.006}
#'
#'   Xie, B., Pan, W., & Shen, X. (2008). Variable selection in penalized
#'   model-based clustering via regularization on grouped parameters.
#'   \emph{Biometrics}, 64(3), 921-930. \doi{10.1111/j.1541-0420.2007.00955.x}
#'
#'   Zou, H. (2006). The adaptive lasso and its oracle properties.
#'   \emph{Journal of the American Statistical Association}, 101(476),
#'   1418-1429. \doi{10.1198/016214506000000735}
#' @seealso \code{\link{estimate_gmm_robust}} (number of classes and
#'   penalty selection), \code{\link{blrt_gmm_robust}},
#'   \code{\link{plot_robust_gmm}}, \code{\link{bch_robust}} (relating the
#'   classes to baseline or distal variables), \code{\link{plot_mcmc_chains}}.
#' @examples
#' data(neuro_long)
#' set.seed(1)
#' fit <- robust_gmm(neuro_long, id = "ID", time = "Year",
#'                   outcomes = c("Memory", "Executive"), G = 2,
#'                   robust_method = "t", n_starts = 2)
#' fit
#' summary(fit)
#' plot_robust_gmm(fit)
#'
#' \donttest{
#' # Sparse class differences (group LASSO by outcome), relaxed refit
#' # (adaptive weights: lambda = z^2 / N removes groups with |z| below ~z)
#' N <- length(unique(neuro_long$ID))
#' fit_l <- robust_gmm(neuro_long, id = "ID", time = "Year",
#'                     outcomes = c("Memory", "Executive", "Speed"), G = 3,
#'                     robust_method = "t", lambda_diff = 9 / N,
#'                     group_diff = TRUE, relax = TRUE, n_starts = 2)
#' fit_l$penalty$outcome_selected
#'
#' # MCMC engine
#' fit_b <- robust_gmm(neuro_long, id = "ID", time = "Year",
#'                     outcomes = "Memory", G = 2, engine = "MCMC",
#'                     robust_method = "t", mcmc_iter = 400, n_chains = 2)
#' summary(fit_b)
#' }
#' @export
robust_gmm <- function(data, id, time, outcomes, G, degree = 1,
                       random = c("slope", "intercept", "none"),
                       re_structure = c("full", "block", "diagonal"),
                       re_cov = c("equal", "varying"),
                       resid_var = c("equal", "varying"),
                       engine = "EM", robust = TRUE, robust_method = c("huber", "t"),
                       nu = NULL, alpha = 0.05,
                       lambda_growth = 0, lambda_diff = 0, group_diff = FALSE, adaptive = TRUE,
                       relax = FALSE,
                       n_starts = 5, max_iter = 500, tol = 1e-8, init = c("kmeans", "random"),
                       mcmc_iter = 2000, n_chains = 4, cores = 1) {
  random <- match.arg(random)
  re_structure <- match.arg(re_structure)
  re_cov <- match.arg(re_cov)
  resid_var <- match.arg(resid_var)
  robust_method <- match.arg(robust_method)
  init <- match.arg(init)
  if (!(length(engine) == 1 && engine %in% c("EM", "MCMC"))) stop("`engine` must be either 'EM' or 'MCMC'.")
  if (!is.numeric(G) || length(G) != 1 || G < 1 || G != round(G)) stop("`G` must be a single positive integer.")
  if (!is.numeric(degree) || length(degree) != 1 || degree < 0 || degree > 3 || degree != round(degree)) {
    stop("`degree` must be 0, 1, 2 or 3.")
  }
  if (!is.logical(robust) || length(robust) != 1 || is.na(robust)) stop("`robust` must be a single TRUE/FALSE value.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("`alpha` must be strictly between 0 and 1.")
  if (!is.null(nu) && (!is.numeric(nu) || length(nu) != 1 || !is.finite(nu) || nu <= 0)) {
    stop("`nu` must be NULL (estimate it) or a single positive number.")
  }
  for (nm in c("lambda_growth", "lambda_diff")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 0) stop("`", nm, "` must be a single non-negative number.")
  }
  if (!is.numeric(n_starts) || length(n_starts) != 1 || n_starts < 1) stop("`n_starts` must be at least 1.")
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) stop("`max_iter` must be a positive integer.")
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) stop("`tol` must be positive.")
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) stop("`cores` must be a single positive integer.")
  G <- as.integer(G)
  degree <- as.integer(degree)
  nrand <- switch(random, none = 0L, intercept = 1L, slope = 2L)
  if (nrand > degree + 1) stop("`random = \"slope\"` requires `degree >= 1`.")
  if (engine == "MCMC") {
    if (!is.numeric(mcmc_iter) || length(mcmc_iter) != 1 || mcmc_iter < 4) stop("`mcmc_iter` must be at least 4.")
    if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains < 1) stop("`n_chains` must be a positive integer.")
  }

  prep <- .gmm_prepare(data, id, time, outcomes)
  if (G > prep$N) stop("`G` cannot exceed the number of persons.")
  method <- .robust_method_name(robust, robust_method)

  call_args <- list(
    id = id, time = time, outcomes = outcomes, G = G, degree = degree, random = random,
    re_structure = re_structure, re_cov = re_cov, resid_var = resid_var, engine = engine,
    robust = robust, robust_method = robust_method, nu = nu, alpha = alpha,
    lambda_growth = lambda_growth, lambda_diff = lambda_diff, group_diff = group_diff,
    adaptive = adaptive, relax = relax,
    n_starts = n_starts, max_iter = max_iter, tol = tol, init = init,
    mcmc_iter = mcmc_iter, n_chains = n_chains, cores = cores
  )
  spec <- list(
    G = G, degree = degree, nrand = nrand, random = random, re_structure = re_structure,
    re_cov = re_cov, resid_var = resid_var, method = method, nu = nu, alpha = alpha,
    max_iter = max_iter, tol = tol, init = init, K = prep$K, outcomes = prep$outcomes,
    id = id, time = time, terms = .gmm_term_names(time, degree),
    nu_estimated = identical(method, "t") && is.null(nu)
  )
  idx <- .gmm_penalty_index(G, prep$K, degree)
  if (G == 1 && lambda_diff > 0) {
    warning("`lambda_diff` has no effect with a single class; ignoring it.", call. = FALSE)
    lambda_diff <- 0
  }
  penalized <- lambda_growth > 0 || lambda_diff > 0
  pen0 <- .gmm_penalty_spec(0, 0, FALSE, idx, G, prep$K, degree)

  # Penalized fits start from the unpenalized maximum-likelihood solution
  # (warm start), which also provides the adaptive-Lasso weights.
  res0 <- .gmm_em(prep, spec, pen0, idx, n_starts = n_starts, cores = cores)
  pen <- pen0
  if (penalized) {
    pen <- .gmm_penalty_spec(lambda_growth, lambda_diff, isTRUE(group_diff), idx, G, prep$K, degree,
                             beta_init = if (isTRUE(adaptive)) res0$beta else NULL,
                             info_scale = if (isTRUE(group_diff)) .gmm_info_scale(prep, spec, res0) else NULL)
  }
  warm <- function(r) list(beta = r$beta, D = r$D, sig2 = r$sig2, pi = r$pi, nu = r$nu)

  if (engine == "MCMC") {
    if (isTRUE(relax) && penalized) warning("`relax` applies to the EM engine only; ignoring it.", call. = FALSE)
    res_em <- if (penalized) suppressWarnings(.gmm_em(prep, spec, pen, idx, 1, 1, start = warm(res0))) else res0
    return(.robust_gmm_mcmc(prep, spec, pen, idx, res_em, mcmc_iter, n_chains, cores, call_args))
  }

  res <- if (penalized) .gmm_em(prep, spec, pen, idx, n_starts = 1, cores = 1, start = warm(res0)) else res0
  pen_used <- pen
  if (isTRUE(relax) && penalized) {
    pattern <- .gmm_pattern(res$beta, idx, pen)
    pen_used <- .gmm_penalty_spec(0, 0, FALSE, idx, G, prep$K, degree, constrain = pattern)
    res <- .gmm_em(prep, spec, pen_used, idx, n_starts = 1, cores = 1, start = warm(res))
    res$pattern <- pattern
  }
  .gmm_finalize(prep, spec, res, pen, pen_used, idx, engine = "EM", call_args = call_args, data = data,
                extra_penalty = list(unpenalized_loglik = res0$loglik - sum(log(prep$scale[prep$outc + 1L]))))
}

#' Zero / Fusion Pattern Selected by the Growth-Mixture Penalties
#' @keywords internal
#' @noRd
.gmm_pattern <- function(beta, idx, pen, tol = 1e-8) {
  b <- as.vector(t(beta))
  growth_zero <- if (pen$lambda_growth > 0) abs(b[idx$growth]) < tol else rep(FALSE, length(idx$growth))
  diff_zero <- rep(FALSE, length(b))
  if (pen$lambda_diff > 0) {
    dev <- as.vector(t(sweep(beta, 2, colMeans(beta), "-")))
    diff_zero <- abs(dev) < tol * pmax(1, abs(b))
  }
  list(growth_zero = growth_zero, diff_zero = diff_zero)
}

#' EM (ECM) Estimation of the Robust Growth Mixture Model
#'
#' Each iteration is an alternating ECM (AECM) iteration with two cycles
#' (see the comments in the code): (penalized) weighted GLS for the fixed
#' effects and the proportions, then a fresh E-step and the updates of the
#' residual variances, random-effect covariances, proportions and (t model)
#' an ECME step for nu.
#' For the likelihood-based fits (Gaussian and t) the iterations are
#' accelerated by the SQUAREM scheme of Varadhan & Roland (2008): two ECM
#' steps are extrapolated on an unconstrained parameterization (log
#' variances, log-Cholesky covariances, logit proportions) and the result is
#' accepted only if, after a stabilizing ECM step, the (penalized)
#' log-likelihood is at least that of the plain double ECM step, so the
#' algorithm stays monotone. The Huber fit is a plain fixed-point iteration.
#'
#' @param prep Result of \code{.gmm_prepare()}.
#' @param spec Model specification (see \code{robust_gmm()}).
#' @param pen Penalty specification; \code{idx} its index structure.
#' @param n_starts,cores Number of starts and cores.
#' @param start Optional list of starting values (single start).
#' @return A list with the standardized-scale estimates of the best start.
#' @keywords internal
#' @noRd
.gmm_em <- function(prep, spec, pen, idx, n_starts, cores, start = NULL) {
  G <- spec$G
  K <- prep$K
  degree <- spec$degree
  nrand <- spec$nrand
  Q <- K * nrand
  method <- spec$method
  alpha <- spec$alpha
  t_scatter <- identical(method, "t")
  monotone <- !identical(method, "huber")
  penalized <- pen$lambda_growth > 0 || pen$lambda_diff > 0
  ols <- if (is.null(start)) .gmm_person_ols(prep, degree) else NULL
  N <- prep$N
  tri <- if (Q > 0) lower.tri(diag(Q), diag = TRUE) else NULL

  estep <- function(th) {
    lf <- .gmm_logf(prep, th$beta, th$D, th$sig2, th$pi, method, th$nu, degree, nrand)
    post <- .posterior_from_logf(lf$logf)
    list(lf = lf, z = post$z, loglik = post$loglik,
         obj = post$loglik - N * .gmm_penalty_value(th$beta, pen, idx))
  }

  # Alternating ECM (AECM; Meng & van Dyk, 1997) with two cycles:
  # cycle 1 treats (class, latent scale) as missing and updates the class
  # proportions and the fixed effects by (penalized) weighted GLS given the
  # current covariance parameters; cycle 2 re-runs the E-step and treats the
  # random effects as missing too, updating the residual variances, the
  # random-effect covariances, the proportions and nu. Each cycle increases
  # the (penalized) observed-data log-likelihood.
  mstep <- function(th, e) {
    Dfun <- function(g, D) if (nrand > 0) D[[g]] else matrix(0, 1, 1)
    # ---- cycle 1: proportions and fixed effects ---------------------------
    z <- e$z
    W <- matrix(1, N, G)
    for (g in seq_len(G)) W[, g] <- .gmm_weights(e$lf$maha[, g], e$lf$nobs, method, alpha, th$nu %||% 4)
    pi_new <- pmax(colMeans(z), 1e-10)
    pi_new <- pi_new / sum(pi_new)
    A_list <- vector("list", G)
    c_list <- vector("list", G)
    for (g in seq_len(G)) {
      gl <- gmm_class_gls_cpp(prep$y, prep$time, prep$outc, prep$starts, Dfun(g, th$D),
                              as.numeric(th$sig2[g, ]), K, degree, nrand, z[, g], W[, g])
      A_list[[g]] <- gl$A
      c_list[[g]] <- as.numeric(gl$c)
    }
    beta_new <- th$beta
    if (penalized) {
      beta_new <- .gmm_admm(A_list, c_list, th$beta, pen, idx, N)
    } else {
      for (g in seq_len(G)) {
        if (sum(z[, g]) < 1e-8) next
        Ag <- A_list[[g]]
        ridge <- 1e-10 * max(mean(diag(Ag)), 1e-10)
        beta_new[g, ] <- solve(Ag + diag(ridge, nrow(Ag)), c_list[[g]])
      }
    }
    th_mid <- list(beta = beta_new, D = th$D, sig2 = th$sig2, pi = pi_new, nu = th$nu)

    # ---- cycle 2: variance components, proportions, nu ------------------------
    e2 <- estep(th_mid)
    z2 <- e2$z
    ss <- vector("list", G)
    for (g in seq_len(G)) {
      W[, g] <- .gmm_weights(e2$lf$maha[, g], e2$lf$nobs, method, alpha, th$nu %||% 4)
      ss[[g]] <- gmm_class_suffstats_cpp(prep$y, prep$time, prep$outc, prep$starts,
                                         as.numeric(beta_new[g, ]), Dfun(g, th$D), as.numeric(th$sig2[g, ]),
                                         K, degree, nrand, z2[, g], W[, g], t_scatter)
    }
    num <- matrix(0, G, K)
    den <- matrix(0, G, K)
    for (g in seq_len(G)) {
      rs <- gmm_class_resid_cpp(prep$y, prep$time, prep$outc, prep$starts, as.numeric(beta_new[g, ]),
                                ss[[g]]$bhat, ss[[g]]$zcz, K, degree, nrand, z2[, g], W[, g], t_scatter)
      num[g, ] <- rs$num
      den[g, ] <- rs$den
    }
    sig2 <- th$sig2
    if (spec$resid_var == "equal") {
      pooled <- colSums(num) / pmax(colSums(den), 1e-12)
      sig2 <- matrix(pmax(pooled, 1e-8), G, K, byrow = TRUE)
    } else {
      for (g in seq_len(G)) {
        ok <- den[g, ] > 1e-10
        sig2[g, ok] <- pmax(num[g, ok] / den[g, ok], 1e-8)
      }
    }
    D <- th$D
    if (nrand > 0) {
      D <- .gmm_update_D(lapply(ss, `[[`, "Dnum"), lapply(ss, `[[`, "Dden"), D, spec$re_cov,
                         spec$re_structure, K, nrand)
    }
    pi_new <- pmax(colMeans(z2), 1e-10)
    pi_new <- pi_new / sum(pi_new)
    nu_new <- th$nu
    if (spec$nu_estimated) {
      lf_new <- .gmm_logf(prep, beta_new, D, sig2, pi_new, "t", nu_new, degree, nrand)
      objf <- function(lv) .gmm_loglik_nu(exp(lv), lf_new$maha, lf_new$logdet, lf_new$nobs, pi_new)
      opt <- stats::optimize(objf, interval = log(c(1, 200)), maximum = TRUE, tol = 1e-4)
      if (opt$objective > objf(log(nu_new))) nu_new <- exp(opt$maximum)
    }
    list(beta = beta_new, D = D, sig2 = sig2, pi = pi_new, nu = nu_new)
  }

  # ---- unconstrained parameterization for SQUAREM ------------------------------
  to_vec <- function(th) {
    v <- c(as.vector(th$beta), log(as.vector(th$sig2)))
    if (nrand > 0) {
      for (g in seq_len(G)) {
        L <- t(chol(.force_pd(th$D[[g]], min_ratio = 1e-10)))
        diag(L) <- log(diag(L))
        v <- c(v, L[tri])
      }
    }
    v <- c(v, log(th$pi / th$pi[G]))
    if (spec$nu_estimated) v <- c(v, log(th$nu - 0.999))
    v
  }
  from_vec <- function(v) {
    pos <- 0
    take <- function(n) { out <- v[pos + seq_len(n)]; pos <<- pos + n; out }
    beta <- matrix(take(G * idx$P), G, idx$P)
    sig2 <- matrix(exp(take(G * K)), G, K)
    D <- NULL
    if (nrand > 0) {
      D <- lapply(seq_len(G), function(g) {
        L <- matrix(0, Q, Q)
        L[tri] <- take(sum(tri))
        diag(L) <- exp(diag(L))
        L %*% t(L)
      })
    }
    lp <- take(G)
    pi_g <- exp(lp - max(lp))
    pi_g <- pi_g / sum(pi_g)
    nu_v <- if (spec$nu_estimated) min(0.999 + exp(take(1)), 200) else spec$nu
    list(beta = beta, D = D, sig2 = pmax(sig2, 1e-8), pi = pi_g, nu = nu_v)
  }

  run_one_start <- function(start_id) {
    th <- if (!is.null(start)) start else .gmm_init(prep, ols, G, degree, nrand, spec$init, spec$re_structure)
    if (nrand == 0) th$D <- NULL
    th$nu <- if (t_scatter) (spec$nu %||% th$nu %||% 10) else NULL
    e <- estep(th)
    n_em <- 0L
    converged <- FALSE
    n_dec <- 0L
    step_max <- 1

    while (n_em < spec$max_iter) {
      obj_prev <- e$obj
      th1 <- mstep(th, e)
      e1 <- estep(th1)
      n_em <- n_em + 1L
      if (!monotone) {
        th <- th1
        e <- e1
        if (abs(e$obj - obj_prev) < spec$tol * (1 + abs(e$obj))) { converged <- TRUE; break }
        next
      }
      if (e1$obj < obj_prev - 1e-6 * (1 + abs(obj_prev))) n_dec <- n_dec + 1L
      if (abs(e1$obj - obj_prev) < spec$tol * (1 + abs(e1$obj))) {
        th <- th1
        e <- e1
        converged <- TRUE
        break
      }
      th2 <- mstep(th1, e1)
      e2 <- estep(th2)
      n_em <- n_em + 1L
      accepted <- FALSE
      v0 <- to_vec(th)
      r <- to_vec(th1) - v0
      vv <- to_vec(th2) - to_vec(th1) - r
      nr <- sqrt(sum(r^2))
      nv <- sqrt(sum(vv^2))
      if (is.finite(nr) && is.finite(nv) && nv > 1e-12) {
        a <- -nr / nv
        a <- max(min(a, -1), -step_max)
        th_ext <- tryCatch(from_vec(v0 - 2 * a * r + a^2 * vv), error = function(err) NULL)
        if (!is.null(th_ext)) {
          e_ext <- tryCatch(estep(th_ext), error = function(err) NULL)
          if (!is.null(e_ext) && is.finite(e_ext$obj)) {
            th3 <- mstep(th_ext, e_ext)
            e3 <- estep(th3)
            n_em <- n_em + 1L
            if (is.finite(e3$obj) && e3$obj >= e2$obj) {
              th <- th3
              e <- e3
              accepted <- TRUE
              if (a == -step_max) step_max <- step_max * 4
            }
          }
        }
      }
      if (!accepted) {
        th <- th2
        e <- e2
        step_max <- max(1, step_max / 4)
      }
      if (e$obj < obj_prev - 1e-6 * (1 + abs(obj_prev))) n_dec <- n_dec + 1L
      if (abs(e$obj - obj_prev) < spec$tol * (1 + abs(e$obj))) { converged <- TRUE; break }
    }

    warn <- character(0)
    if (monotone && n_dec > 0) {
      warn <- sprintf(
        "Start %d: the (penalized) log-likelihood decreased in %d step(s) (numerical safeguards on near-singular covariances). Consider more `n_starts`, a more constrained covariance (`re_cov = \"equal\"`, `re_structure = \"diagonal\"`), or fewer classes.",
        start_id, n_dec)
    }
    list(obj = e$obj, loglik = e$loglik, beta = th$beta, D = th$D, sig2 = th$sig2, pi = th$pi,
         nu = th$nu, z = e$z, converged = converged, iterations = n_em, warnings = warn)
  }

  results <- .run_parallel(cores, if (is.null(start)) n_starts else 1L, run_one_start,
                           export_vars = character(0), export_env = environment())
  objs <- vapply(results, function(r) r$obj, numeric(1))
  if (!any(is.finite(objs))) stop("All EM starts failed to produce a finite log-likelihood.")
  best <- results[[which.max(objs)]]
  for (w in best$warnings) warning(w, call. = FALSE)
  best
}

#' Assemble a robust_gmm Object
#' @keywords internal
#' @noRd
.gmm_finalize <- function(prep, spec, res, pen, pen_used, idx, engine, call_args, data = NULL,
                          extra = list(), extra_penalty = list()) {
  G <- spec$G
  K <- prep$K
  degree <- spec$degree
  nrand <- spec$nrand
  nb <- degree + 1
  Q <- K * nrand
  N <- prep$N
  method <- spec$method
  t_scatter <- identical(method, "t")

  lf <- .gmm_logf(prep, res$beta, res$D, res$sig2, res$pi, method, res$nu, degree, nrand)
  post <- .posterior_from_logf(lf$logf)
  z <- post$z
  jac <- sum(log(prep$scale[prep$outc + 1L]))
  loglik <- post$loglik - jac

  W <- matrix(1, N, G)
  bh_list <- vector("list", G)
  for (g in seq_len(G)) {
    W[, g] <- .gmm_weights(lf$maha[, g], lf$nobs, method, spec$alpha, res$nu %||% 4)
    if (nrand > 0) {
      bh_list[[g]] <- gmm_class_suffstats_cpp(prep$y, prep$time, prep$outc, prep$starts,
                                              as.numeric(res$beta[g, ]), res$D[[g]], as.numeric(res$sig2[g, ]),
                                              K, degree, nrand, z[, g], W[, g], t_scatter)$bhat
    }
  }
  obs_weights <- rowSums(z * W)
  assignments <- max.col(z, ties.method = "first")

  bt <- .gmm_backtransform(res$beta, res$D, res$sig2, prep, degree, nrand)
  terms <- spec$terms
  outcomes <- prep$outcomes
  re_names <- if (nrand > 0) paste(rep(outcomes, each = nrand), rep(terms[seq_len(nrand)], K), sep = ":") else NULL
  coefficients <- lapply(seq_len(G), function(g) {
    matrix(bt$beta[g, ], nrow = K, byrow = TRUE, dimnames = list(outcomes, terms))
  })
  names(coefficients) <- paste0("Class_", seq_len(G))
  random_cov <- NULL
  random_effects <- NULL
  if (nrand > 0) {
    random_cov <- lapply(bt$D, function(M) { dimnames(M) <- list(re_names, re_names); M })
    names(random_cov) <- paste0("Class_", seq_len(G))
    sq <- rep(prep$scale, each = nrand)
    random_effects <- matrix(vapply(seq_len(N), function(i) bh_list[[assignments[i]]][i, ] * sq, numeric(Q)),
                             N, Q, byrow = TRUE)
    dimnames(random_effects) <- list(NULL, re_names)
  }
  residual_var <- bt$sig2
  dimnames(residual_var) <- list(paste0("Class_", seq_len(G)), outcomes)

  # ---- parameter count --------------------------------------------------------
  use_growth <- pen_used$lambda_growth > 0
  use_diff <- pen_used$lambda_diff > 0
  df_beta <- .gmm_beta_df(res$beta, idx, use_growth, use_diff)
  n_D <- if (nrand == 0) 0 else switch(spec$re_structure,
    full = Q * (Q + 1) / 2, block = K * nrand * (nrand + 1) / 2, diagonal = Q)
  n_par <- (G - 1) + df_beta + n_D * (if (spec$re_cov == "varying") G else 1) +
    K * (if (spec$resid_var == "varying") G else 1) + as.integer(spec$nu_estimated)
  z_safe <- pmax(z, 1e-15)
  entropy <- if (G == 1) 1 else 1 - sum(-z * log(z_safe)) / (N * log(G))
  sizes <- tabulate(assignments, nbins = G)
  fit <- data.frame(
    Classes = G, LogLik = loglik, Parameters = n_par,
    AIC = -2 * loglik + 2 * n_par, BIC = -2 * loglik + n_par * log(N),
    SABIC = -2 * loglik + n_par * log((N + 2) / 24), Entropy = entropy,
    Min_Size = min(sizes) / N, Max_Size = max(sizes) / N,
    Lambda_Growth = pen$lambda_growth, Lambda_Diff = pen$lambda_diff
  )

  # ---- penalty summary ----------------------------------------------------------
  pattern <- res$pattern %||% if (pen$lambda_growth > 0 || pen$lambda_diff > 0) .gmm_pattern(res$beta, idx, pen) else NULL
  outcome_selected <- NULL
  if (!is.null(pattern) && pen$lambda_diff > 0) {
    outcome_selected <- vapply(seq_len(K), function(k) !all(pattern$diff_zero[idx$outcome == k]), logical(1))
    names(outcome_selected) <- outcomes
  }
  penalty <- c(list(lambda_growth = pen$lambda_growth, lambda_diff = pen$lambda_diff,
                    group_diff = pen$group, adaptive = isTRUE(pen$adaptive),
                    relaxed = isTRUE(pen_used$constraint), pattern = pattern,
                    outcome_selected = outcome_selected), extra_penalty)

  if (min(sizes) / N < 0.02) {
    warning(sprintf("The smallest class holds only %.1f%% of the persons; check for a spurious solution (more `n_starts`, fewer classes, or a more constrained covariance).",
                    100 * min(sizes) / N), call. = FALSE)
  }

  gram <- lapply(seq_len(K), function(k) {
    tt <- prep$time[prep$outc == k - 1L]
    B <- outer(tt, 0:degree, `^`)
    crossprod(B) / length(tt)
  })

  out <- c(
    list(
      engine = engine, robust_method = method, nu = if (t_scatter) res$nu else NULL,
      coefficients = coefficients, random_cov = random_cov, residual_var = residual_var,
      proportions = as.numeric(res$pi), probabilities = z, assignments = assignments,
      ids = prep$ids, weights = obs_weights, random_effects = random_effects, fit = fit,
      penalty = penalty, converged = res$converged, iterations = res$iterations
    ),
    extra,
    list(
      spec = spec,
      internal = list(beta = res$beta, D = res$D, sig2 = res$sig2, pi = res$pi, nu = res$nu,
                      gram = gram, center = prep$center, scale = prep$scale),
      call_args = call_args,
      data = prep$raw
    )
  )
  class(out) <- "robust_gmm"
  out
}

#' Optimal Relabeling of Growth-Mixture Classes
#'
#' Matches the classes of \code{beta_new} to those of \code{beta_ref}
#' (standardized coefficient matrices, classes x coefficients) by an exact
#' optimal assignment minimizing the mean squared distance between class
#' mean trajectories over the observed times of every outcome.
#' @return Integer vector: element g is the class of \code{beta_new} matched
#'   to reference class g.
#' @keywords internal
#' @noRd
.gmm_match_classes <- function(beta_ref, beta_new, gram, degree) {
  G <- nrow(beta_ref)
  if (G <= 1) return(seq_len(G))
  nb <- degree + 1
  K <- length(gram)
  cost <- matrix(0, G, G)
  for (g in seq_len(G)) for (h in seq_len(G)) {
    for (k in seq_len(K)) {
      sel <- (k - 1) * nb + seq_len(nb)
      d <- beta_new[h, sel] - beta_ref[g, sel]
      cost[g, h] <- cost[g, h] + as.numeric(t(d) %*% gram[[k]] %*% d)
    }
  }
  as.integer(solve_lsap_cpp(cost)) + 1L
}

#' MCMC Branch of robust_gmm()
#' @keywords internal
#' @noRd
.robust_gmm_mcmc <- function(prep, spec, pen, idx, em, mcmc_iter, n_chains, cores, call_args) {
  G <- spec$G
  K <- prep$K
  degree <- spec$degree
  nrand <- spec$nrand
  nb <- degree + 1
  P <- K * nb
  Q <- K * nrand
  N <- prep$N
  method <- spec$method

  sd_pow <- vapply(0:degree, function(l) if (l == 0) 1 else max(stats::sd(prep$time^l), 1e-8), numeric(1))
  beta_prior_sd <- rep(10 / sd_pow, K)
  re_prior_scale <- if (nrand > 0) rep(5 / sd_pow[seq_len(nrand)], K) else numeric(1)
  robust_type <- switch(method, none = 0L, huber = 1L, t = 2L)
  nu_start <- if (identical(method, "t")) (spec$nu %||% em$nu %||% 10) else 30
  D0 <- if (nrand > 0) lapply(em$D, unname) else rep(list(matrix(0, 1, 1)), G)
  structure_code <- switch(spec$re_structure, full = 0L, block = 1L, diagonal = 2L)
  rates <- .gmm_mcmc_rates(pen, idx, G, K, N)
  jitter_sd <- beta_prior_sd / 100

  run_one_chain <- function(chain_id) {
    beta0 <- em$beta + matrix(stats::rnorm(G * P), G, P) * matrix(jitter_sd, G, P, byrow = TRUE)
    gmm_mcmc_chain_cpp(
      prep$y, prep$time, prep$outc, prep$starts, K, degree, nrand, G, mcmc_iter,
      robust_type, spec$alpha, nu_start, spec$nu_estimated,
      spec$re_cov == "equal", spec$resid_var == "equal", structure_code,
      rates$growth, rates$diff, rates$group, rates$group_scale, isTRUE(pen$group),
      beta_prior_sd, re_prior_scale, beta0, D0, em$sig2, pmax(em$pi, 1e-3)
    )
  }
  chains <- .run_parallel(cores, n_chains, run_one_chain, export_vars = character(0),
                          export_env = environment())

  gram <- lapply(seq_len(K), function(k) {
    tt <- prep$time[prep$outc == k - 1L]
    B <- outer(tt, 0:degree, `^`)
    crossprod(B) / length(tt)
  })

  # ---- relabel every draw to the EM solution ------------------------------------
  for (ci in seq_along(chains)) {
    ch <- chains[[ci]]
    for (it in seq_len(mcmc_iter)) {
      Bm <- matrix(ch$beta_chain[it, ], G, P, byrow = TRUE)
      perm <- .gmm_match_classes(em$beta, Bm, gram, degree)
      if (!identical(perm, seq_len(G))) {
        ch$beta_chain[it, ] <- as.vector(t(Bm[perm, , drop = FALSE]))
        ch$sig2_chain[it, ] <- as.vector(t(matrix(ch$sig2_chain[it, ], G, K, byrow = TRUE)[perm, , drop = FALSE]))
        ch$pi_chain[it, ] <- ch$pi_chain[it, perm]
        if (Q > 0) {
          Dm <- matrix(ch$D_chain[it, ], Q * Q, G)
          ch$D_chain[it, ] <- as.vector(Dm[, perm, drop = FALSE])
        }
      }
    }
    chains[[ci]] <- ch
  }

  burnin <- floor(mcmc_iter / 2)
  keep <- (burnin + 1):mcmc_iter
  pool <- function(name) do.call(rbind, lapply(chains, function(ch) ch[[name]][keep, , drop = FALSE]))
  beta_hat <- matrix(colMeans(pool("beta_chain")), G, P, byrow = TRUE)
  sig2_hat <- matrix(colMeans(pool("sig2_chain")), G, K, byrow = TRUE)
  pi_hat <- colMeans(pool("pi_chain"))
  D_hat <- NULL
  if (Q > 0) {
    Dm <- matrix(colMeans(pool("D_chain")), Q * Q, G)
    D_hat <- lapply(seq_len(G), function(g) matrix(Dm[, g], Q, Q))
  }
  nu_hat <- if (identical(method, "t")) stats::median(unlist(lapply(chains, function(ch) ch$nu_chain[keep]))) else NULL

  res <- list(beta = beta_hat, D = D_hat, sig2 = sig2_hat, pi = pi_hat, nu = nu_hat,
              converged = NA, iterations = NA)

  # ---- WAIC -------------------------------------------------------------------------
  jac_person <- as.numeric(tapply(log(prep$scale[prep$outc + 1L]), factor(prep$pid, levels = seq_len(N)), sum))
  jac_person[is.na(jac_person)] <- 0
  draws <- expand.grid(it = keep, ch = seq_along(chains))
  if (nrow(draws) > 500) draws <- draws[unique(round(seq(1, nrow(draws), length.out = 500))), , drop = FALSE]
  ll <- matrix(NA_real_, nrow(draws), N)
  for (s in seq_len(nrow(draws))) {
    ch <- chains[[draws$ch[s]]]
    it <- draws$it[s]
    Bs <- matrix(ch$beta_chain[it, ], G, P, byrow = TRUE)
    Ss <- matrix(ch$sig2_chain[it, ], G, K, byrow = TRUE)
    Ds <- if (Q > 0) { Dm <- matrix(ch$D_chain[it, ], Q * Q, G); lapply(seq_len(G), function(g) matrix(Dm[, g], Q, Q)) } else NULL
    lf <- .gmm_logf(prep, Bs, Ds, Ss, ch$pi_chain[it, ], method, ch$nu_chain[it], degree, nrand)
    ll[s, ] <- .row_logsumexp(lf$logf) - jac_person
  }
  m <- apply(ll, 2, max)
  lppd <- sum(m + log(colMeans(exp(sweep(ll, 2, m, "-")))))
  p_waic <- sum(apply(ll, 2, stats::var))
  waic <- list(WAIC = -2 * (lppd - p_waic), lppd = lppd, p_waic = p_waic, n_draws = nrow(draws))

  # ---- draws on the original scale, for diagnostics and trace plots ---------------------
  draws_array <- .gmm_draws_array(chains, spec, prep, P, Q)

  out <- .gmm_finalize(prep, spec, res, pen, pen, idx, engine = "MCMC", call_args = call_args,
                       data = NULL)
  out$fit$WAIC <- waic$WAIC
  diagnostics <- .diagnostics_from_array(draws_array[keep, , , drop = FALSE])
  if (!is.null(diagnostics) && isTRUE(any(diagnostics$Rhat > 1.1, na.rm = TRUE))) {
    warning("Some MCMC parameters have a Gelman-Rubin R-hat > 1.1; consider increasing `mcmc_iter`. See `model$mcmc_diagnostics`.",
            call. = FALSE)
  }
  out$mcmc_draws <- list(draws = draws_array, n_chains = length(chains), mcmc_iter = mcmc_iter, burnin = burnin)
  out$mcmc_diagnostics <- diagnostics
  out$waic <- waic
  out <- out[c(setdiff(names(out), c("spec", "internal", "call_args", "data")), "spec", "internal", "call_args", "data")]
  class(out) <- "robust_gmm"
  out
}

#' Original-Scale Draws Array of a Growth-Mixture MCMC Fit
#'
#' Builds an \code{[iterations, chains, parameters]} array with the class
#' trajectories (\code{"beta[g,outcome:term]"}), random-effect variances
#' (\code{"D[g,outcome:term]"}), residual variances
#' (\code{"sigma2[g,outcome]"}), mixing proportions (\code{"pi[g]"}) and,
#' if estimated, \code{"nu"}. Shared variances appear once (class index
#' omitted).
#' @keywords internal
#' @noRd
.gmm_draws_array <- function(chains, spec, prep, P, Q) {
  G <- spec$G
  K <- prep$K
  nb <- spec$degree + 1
  nrand <- spec$nrand
  terms <- spec$terms
  outcomes <- prep$outcomes
  coef_names <- paste(rep(outcomes, each = nb), rep(terms, K), sep = ":")
  sc_beta <- rep(prep$scale, each = nb)
  shift <- rep(0, P)
  shift[(seq_len(K) - 1) * nb + 1] <- prep$center
  iters <- nrow(chains[[1]]$beta_chain)

  build <- function(ch) {
    cols <- list()
    for (g in seq_len(G)) {
      b <- ch$beta_chain[, (g - 1) * P + seq_len(P), drop = FALSE]
      b <- sweep(sweep(b, 2, sc_beta, "*"), 2, shift, "+")
      colnames(b) <- sprintf("beta[%d,%s]", g, coef_names)
      cols[[length(cols) + 1]] <- b
    }
    if (Q > 0) {
      re_names <- paste(rep(outcomes, each = nrand), rep(terms[seq_len(nrand)], K), sep = ":")
      sq <- rep(prep$scale, each = nrand)
      diag_pos <- (seq_len(Q) - 1) * Q + seq_len(Q)
      gs <- if (spec$re_cov == "equal") 1L else seq_len(G)
      for (g in gs) {
        d <- ch$D_chain[, (g - 1) * Q * Q + diag_pos, drop = FALSE]
        d <- sweep(d, 2, sq^2, "*")
        colnames(d) <- if (spec$re_cov == "equal") sprintf("D[%s]", re_names) else sprintf("D[%d,%s]", g, re_names)
        cols[[length(cols) + 1]] <- d
      }
    }
    gs <- if (spec$resid_var == "equal") 1L else seq_len(G)
    for (g in gs) {
      s <- ch$sig2_chain[, (g - 1) * K + seq_len(K), drop = FALSE]
      s <- sweep(s, 2, prep$scale^2, "*")
      colnames(s) <- if (spec$resid_var == "equal") sprintf("sigma2[%s]", outcomes) else sprintf("sigma2[%d,%s]", g, outcomes)
      cols[[length(cols) + 1]] <- s
    }
    pm <- ch$pi_chain
    colnames(pm) <- sprintf("pi[%d]", seq_len(G))
    cols[[length(cols) + 1]] <- pm
    if (isTRUE(spec$nu_estimated)) cols[[length(cols) + 1]] <- matrix(ch$nu_chain, ncol = 1, dimnames = list(NULL, "nu"))
    do.call(cbind, cols)
  }
  mats <- lapply(chains, build)
  arr <- array(NA_real_, dim = c(iters, length(chains), ncol(mats[[1]])))
  for (ci in seq_along(mats)) arr[, ci, ] <- mats[[ci]]
  dimnames(arr) <- list(NULL, paste0("chain:", seq_along(chains)), colnames(mats[[1]]))
  arr
}
