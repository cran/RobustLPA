#' Simulated Neuropsychological Dataset for Robust LPA
#'
#' A synthetic dataset of neuropsychological test scores and reaction times
#' for two latent groups, "Healthy" and "Pathological", designed as a
#' worked example for every estimation path in this package: the
#' Expectation-Maximization and MCMC engines, all six variance-covariance
#' parameterizations, robust vs. classical estimation, LASSO regularization,
#' model/profile selection (\code{\link{estimate_profiles_robust}}), the
#' bootstrapped likelihood ratio test (\code{\link{blrt_robust}}), and the
#' BCH auxiliary-variable method (\code{\link{bch_robust}}).
#'
#' @details
#' The two groups differ in more than location: \code{Attention} and
#' \code{Executive_Functions} have \emph{identical} means and variances in
#' both groups (deliberate noise variables, carrying no group signal), while
#' \code{Memory}, \code{RT_Stroop}, and \code{RT_TMT} differ in mean between
#' groups, and \code{RT_Stroop}/\code{RT_TMT} additionally differ in
#' variance and in their correlation with each other (0.35 in Healthy vs.
#' 0.90 in Pathological). This last feature is intentional: it is a genuine,
#' whole-group difference in covariance structure (not merely in means), so
#' that \code{model = 6} (a fully unconstrained covariance matrix per
#' profile) is the best-fitting parameterization for this dataset by BIC at
#' \code{G = 2} -- run \code{estimate_profiles_robust(scale(neuro_data[, 3:7]),
#' models = 1:6, n_profiles = 1:3)} and inspect \code{$fit_table} to see this
#' directly. A more parsimonious model (e.g. \code{model = 3}, a single
#' covariance matrix shared across profiles) fits these data measurably
#' worse, illustrating why the six parameterizations exist and how to choose
#' among them.
#'
#' On top of this, 6\% of the Pathological observations (chosen at random)
#' receive an additional, positive, randomly-sized shift on \code{RT_Stroop}
#' and \code{RT_TMT} (drawn from a Gamma distribution, so the contamination
#' varies in severity rather than landing on a single fixed value) --
#' measurement-error-like outliers on top of the two groups' otherwise
#' multivariate-normal structure. These are what \code{robust = TRUE} (the
#' default of \code{\link{robust_lpa}}) down-weights via Huber-type
#' estimation; compare \code{robust = TRUE} vs. \code{robust = FALSE} fits to
#' see their effect on the estimated Pathological-profile covariance.
#'
#' The contamination magnitude and rate were calibrated (by direct grid
#' search across the six variance-covariance models, replicated over
#' multiple random seeds) so that fitting \code{model = 6} with \code{G = 2}
#' reliably wins by BIC over both more-parsimonious models at \code{G = 2}
#' and less-parsimonious models at \code{G = 3}, and recovers
#' \code{True_Profile} with about 93\% accuracy (classical, Huber or t
#' estimation, \code{model = 6}, \code{G = 2}).
#'
#' @format A data frame with 250 rows and 7 variables:
#' \describe{
#'   \item{ID}{Unique identifier for each participant.}
#'   \item{True_Profile}{The true latent group, \code{"Healthy"} (n = 150) or
#'     \code{"Pathological"} (n = 100). Not used for estimation (LPA is
#'     unsupervised); included so that recovered profiles can be checked
#'     against ground truth, e.g. \code{table(neuro_data$True_Profile, fit$assignments)}.}
#'   \item{Memory}{Simulated memory test score. Differs in mean between groups.}
#'   \item{Attention}{Simulated attention test score. Identical distribution
#'     in both groups (no group signal); a noise variable.}
#'   \item{Executive_Functions}{Simulated executive functions score.
#'     Identical distribution in both groups (no group signal); a noise
#'     variable.}
#'   \item{RT_Stroop}{Reaction time in milliseconds. Differs in mean,
#'     variance, and correlation with \code{RT_TMT} between groups; a subset
#'     of Pathological observations carry an additional outlying shift.}
#'   \item{RT_TMT}{Reaction time in milliseconds. Differs in mean, variance,
#'     and correlation with \code{RT_Stroop} between groups; a subset of
#'     Pathological observations carry an additional outlying shift.}
#' }
#' @source Simulated data for testing and documentation purposes.
"neuro_data"

#' Simulated Longitudinal Neuropsychological Dataset for Robust Growth Mixture Models
#'
#' A synthetic longitudinal dataset (long format: one row per person and
#' annual visit) designed as a worked example for \code{\link{robust_gmm}}:
#' three latent classes of cognitive change, three outcomes, unbalanced
#' follow-up with missed visits and drop-out, and a few gross data-entry
#' errors.
#'
#' @details
#' 400 persons are assigned to three latent classes ("Stable", "Slow
#' decline", "Fast decline"; about 50/30/20\%) and followed for up to six
#' annual visits (\code{Year} 0 to 5). \code{Memory} and \code{Executive}
#' (T-score-like metric) decline at class-specific rates (Memory: 0, -2 and
#' -5 points per year; Executive: 0, -1.5 and -4 points per year), whereas
#' \code{Speed} declines by 0.5 points per year \emph{in every class} (an
#' outcome that does not differentiate the classes, useful to illustrate the
#' group LASSO of \code{robust_gmm(lambda_diff = , group_diff = TRUE)}; the
#' "Stable" class has exactly zero slopes on Memory and Executive, useful
#' to illustrate \code{lambda_growth}). Within classes, persons have
#' correlated random intercepts (SD 5) and slopes (SD 0.4) on every outcome,
#' and residual errors with SD 2.5. After each visit a person drops out
#' with a probability that increases as the last observed Memory score
#' decreases (missing at random); 8\% of the follow-up visits are missed and
#' 5\% of the single test scores are missing. For 4\% of the persons one
#' score is corrupted by a gross error of 25 to 40 points. \code{Age} and
#' \code{Biomarker} are baseline characteristics that differ between the
#' classes, for illustrating \code{\link{bch_robust}} on a
#' \code{robust_gmm()} fit. The generating code (with its fixed seed) is
#' installed with the package: run
#' \code{source(system.file("scripts", "generate_neuro_long.R", package = "RobustLPA"))}
#' to rebuild the data set.
#'
#' @format A data frame with one row per person and visit, and 8 variables:
#' \describe{
#'   \item{ID}{Person identifier.}
#'   \item{Year}{Years since baseline (0 to 5).}
#'   \item{Memory, Executive, Speed}{Simulated test scores (\code{NA} when
#'     not observed).}
#'   \item{True_Class}{The true latent class (not used for estimation).}
#'   \item{Age}{Age at baseline (constant within person).}
#'   \item{Biomarker}{A baseline biomarker level (constant within person).}
#' }
#' @source Simulated; see Details.
#' @examples
#' data(neuro_long)
#' head(neuro_long)
#' table(table(neuro_long$ID))   # number of visits per person
"neuro_long"
