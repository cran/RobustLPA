#' BCH Classification Matrix, Weights, and Bias-Corrected Profile Means
#'
#' Core BCH point-estimate computation, shared by \code{\link{bch_robust}}'s
#' main call and (on each resample) its \code{correction = "bootstrap"} loop.
#' See \code{\link{bch_robust}}'s Details for the underlying formulas.
#'
#' @param z An \code{n x G} matrix of posterior profile-membership probabilities.
#' @param assignments Integer vector of length \code{n}, modal (hard) profile assignments.
#' @param aux_var Numeric vector of length \code{n}, the auxiliary variable (may contain \code{NA}).
#' @param warn_zero_weight Logical; whether to \code{warning()} when a
#'   profile's BCH weights sum to (near) zero. Set \code{FALSE} inside the
#'   bootstrap loop to avoid emitting up to \code{n_boot} warnings.
#' @return A list with \code{D_matrix}, \code{W_matrix}, \code{weighted_means}
#'   (named numeric vector of length \code{G}), \code{assignments_valid},
#'   \code{aux_valid}, and \code{n_valid}.
#' @keywords internal
#' @noRd
.bch_point_estimate <- function(z, assignments, aux_var, warn_zero_weight = TRUE) {
  G <- ncol(z)

  # ---- Classification error matrix D and BCH weight matrix W = D^-1 --------
  # D[t, s] = P(hat_C = s | C = t), estimated as the share of profile t's
  # expected membership (sum_i p_it) that is modally assigned to s:
  #   D[t, s] = sum_{i: hat_C_i = s} p_it / sum_i p_it
  # (Bolck, Croon, & Hagenaars, 2004; Vermunt, 2010). Rows sum to 1.
  # NOTE: averaging posteriors *within* assigned class s (colMeans) would
  # instead estimate P(C = t | hat_C = s), the reverse conditional, which
  # leaves residual bias whenever profiles differ in size or separation.
  class_counts <- table(factor(assignments, levels = 1:G))
  if (any(class_counts == 0)) {
    stop("At least one profile has zero modally-assigned observations; the classification matrix D cannot be estimated. This can happen with highly overlapping or near-empty profiles.")
  }

  D_matrix <- matrix(0, nrow = G, ncol = G)
  z_tot <- colSums(z)
  for (s in 1:G) {
    D_matrix[, s] <- colSums(z[assignments == s, , drop = FALSE]) / z_tot
  }

  W_matrix <- tryCatch(solve(D_matrix), error = function(e) {
    stop("Classification error matrix D is singular (profiles are too overlapping to invert). ",
         "BCH weights cannot be computed for this model; consider a solution with fewer or less-overlapping profiles.")
  })

  # ---- Drop observations with missing aux_var (listwise) -------------------
  valid_idx <- which(!is.na(aux_var))
  n_valid <- length(valid_idx)
  if (n_valid < G) stop("Too few non-missing `aux_var` values to fit a BCH model with ", G, " profiles.")
  assignments_valid <- assignments[valid_idx]
  aux_valid <- aux_var[valid_idx]

  # ---- BCH-corrected profile means -----------------------------------------
  # mu_BCH_g = sum_i W[Chat_i, g] * Y_i / sum_i W[Chat_i, g], summed over ALL
  # n_valid observations -- see bch_robust()'s Details. W = D^-1 maps the
  # distribution over assigned classes back to the true classes, so its rows
  # are indexed by the assigned class and its columns by the true profile.
  weights_per_obs <- t(W_matrix[assignments_valid, , drop = FALSE])  # G x n_valid
  weighted_means <- numeric(G)
  for (g in 1:G) {
    w_g <- weights_per_obs[g, ]
    denom <- sum(w_g)
    if (abs(denom) < 1e-6) {
      if (warn_zero_weight) {
        warning(sprintf("BCH weights for profile %d sum to (near) zero; its corrected mean is unstable and set to NA. This can happen with poorly separated profiles.", g))
      }
      weighted_means[g] <- NA
    } else {
      weighted_means[g] <- sum(w_g * aux_valid) / denom
    }
  }
  names(weighted_means) <- paste0("Profile_", 1:G)

  list(
    D_matrix = D_matrix, W_matrix = W_matrix, weighted_means = weighted_means,
    assignments_valid = assignments_valid, aux_valid = aux_valid, n_valid = n_valid
  )
}

#' Weighted One-Way ANOVA on BCH-Reweighted Long Data
#'
#' Fits the significance test used by \code{\link{bch_robust}}'s main
#' (fixed-\strong{D}) F-test: a one-way weighted ANOVA on the "long" data set
#' (one row per observation per profile, weighted by \eqn{W_{\hat{C}_i,g}}),
#' via direct weighted normal equations rather than \code{stats::lm()}/
#' \code{stats::aov()}, because the BCH weights are frequently negative and
#' base R's weighted-least-squares machinery cannot handle that (it internally
#' scales rows by \code{sqrt(weights)}).
#'
#' @param W_matrix The \code{G x G} BCH weight matrix from \code{\link{.bch_point_estimate}}.
#' @param assignments_valid Integer vector, modal assignments after listwise deletion.
#' @param aux_valid Numeric vector, the auxiliary variable after listwise deletion.
#' @param G Integer, the number of profiles.
#' @param n_valid Integer, the number of observations retained.
#' @return A data.frame with \code{Df}, \code{Sum_Sq}, \code{Mean_Sq},
#'   \code{F_value}, and \code{p_value} for the "Class" and "Residuals" rows.
#' @keywords internal
#' @noRd
.bch_weighted_anova <- function(W_matrix, assignments_valid, aux_valid, G, n_valid) {
  long_Y <- rep(aux_valid, times = G)
  long_Class <- factor(rep(1:G, each = n_valid))
  long_w <- numeric(G * n_valid)
  for (g in 1:G) {
    idx_range <- ((g - 1) * n_valid + 1):(g * n_valid)
    long_w[idx_range] <- W_matrix[assignments_valid, g]
  }

  X_design <- stats::model.matrix(~long_Class)
  XtWX <- t(X_design) %*% (long_w * X_design)
  XtWY <- t(X_design) %*% (long_w * long_Y)

  beta <- tryCatch(solve(XtWX, XtWY), error = function(e) {
    stop("The BCH-weighted design matrix is singular; unable to fit the weighted ANOVA. This can happen with highly degenerate classification weights.")
  })
  fitted_vals <- as.vector(X_design %*% beta)
  resid_vals <- long_Y - fitted_vals

  grand_mean <- sum(long_w * long_Y) / sum(long_w)
  TSS_w <- sum(long_w * (long_Y - grand_mean)^2)
  RSS_w <- sum(long_w * resid_vals^2)
  SS_between <- TSS_w - RSS_w

  df_between <- G - 1
  df_within <- n_valid - G  # total weight sums to n_valid exactly

  if (RSS_w <= 0 || df_within <= 0) {
    warning("Residual sum of squares from the BCH-weighted ANOVA is non-positive (degenerate weighting); the F-test is not interpretable and is set to NA.")
    F_stat <- NA_real_
    p_val <- NA_real_
  } else {
    F_stat <- (SS_between / df_between) / (RSS_w / df_within)
    p_val <- stats::pf(F_stat, df_between, df_within, lower.tail = FALSE)
  }

  data.frame(
    Df = c(df_between, df_within),
    Sum_Sq = c(SS_between, RSS_w),
    Mean_Sq = c(SS_between / df_between, if (df_within > 0) RSS_w / df_within else NA_real_),
    F_value = c(F_stat, NA_real_),
    p_value = c(p_val, NA_real_),
    row.names = c("Class", "Residuals")
  )
}

#' BCH Method for Auxiliary Continuous Variables
#'
#' Applies the 3-step Bolck-Croon-Hagenaars (BCH) method to test the
#' relationship between robust latent profiles and a continuous auxiliary
#' (distal outcome) variable, adjusting for classification error in the
#' profile assignments. Optionally adds a bootstrap correction
#' (\code{correction = "bootstrap"}) for the F-test's main weakness: treating
#' the classification matrix \strong{D} as known/fixed (see Details).
#'
#' @details
#' Let \eqn{\hat{p}_{ig}} be the posterior probability that observation
#' \eqn{i} belongs to profile \eqn{g} (\code{model$probabilities}), and let
#' \eqn{\hat{C}_i} be its modal (hard) assignment (\code{model$assignments}).
#' The classification error matrix \strong{D} is estimated as
#' \deqn{D_{t,s} = P(\hat{C} = s \mid C = t) \approx \frac{\sum_{i:\, \hat{C}_i = s} \hat{p}_{it}}{\sum_{i=1}^{n} \hat{p}_{it}}}
#' (Bolck, Croon, & Hagenaars, 2004; Vermunt, 2010), so each row of
#' \strong{D} sums to 1. The BCH weight matrix is \eqn{W = D^{-1}}, and the
#' classification-error-corrected mean of the auxiliary variable \eqn{Y}
#' for profile \eqn{t} is
#' \deqn{\hat{\mu}^{BCH}_t = \frac{\sum_{i=1}^{n} W_{\hat{C}_i,t} Y_i}{\sum_{i=1}^{n} W_{\hat{C}_i,t}}}
#' summed over \emph{every} observation: each contributes to every profile's
#' mean with a (possibly negative) cross-class weight
#' \eqn{W_{\hat{C}_i,t}}, which is what removes the attenuation bias of a
#' naive comparison of modal-assignment groups. This is unbiased under the
#' BCH assumption that \eqn{Y} is independent of the modal assignment given
#' the true profile.
#'
#' The main (\code{$ANOVA_Table}) significance test is a one-way weighted
#' ANOVA on the equivalent "long" data set (one row per observation per
#' profile, weighted by \eqn{W_{\hat{C}_i,g}}), fit by direct weighted normal
#' equations rather than \code{stats::lm()}/\code{stats::aov()}, because the
#' BCH weights are frequently negative and base R's weighted-least-squares
#' machinery cannot handle that.
#'
#' \strong{The fixed-D caveat, and the bootstrap correction.} \code{$ANOVA_Table}'s
#' F-test treats \strong{D} as known/fixed. Bolck et al. (2004) and Vermunt
#' (2010) both note that this understates the true uncertainty, because
#' \strong{D} is itself estimated from the step-1 model; Vermunt (2010)
#' reports that naive (uncorrected) BCH p-values can be "much too small,"
#' particularly with poorly separated profiles or small samples. The
#' literature's analytic fix is a "sandwich" (pseudo-likelihood) variance
#' correction (Bakk, Oberski, & Vermunt, 2014), which requires the Fisher
#' information of the step-1 mixture log-likelihood -- intractable to derive
#' analytically here for the Huber-robust EM and Bayesian-Lasso MCMC engines.
#' Setting \code{correction = "bootstrap"} instead approximates that
#' correction \emph{nonparametrically}: it resamples observations with
#' replacement, refits the entire step-1 \code{\link{robust_lpa}} model
#' (using the exact same specification as \code{model}, via its stored
#' \code{$call_args}) and recomputes \strong{D}/\strong{W}/the profile means
#' on each resample, so the resulting bootstrap variability genuinely
#' reflects step-1 estimation uncertainty (unlike the fixed-D F-test); for
#' \code{robust_gmm()} fits whole persons (with all their occasions) are
#' resampled (a cluster bootstrap). This
#' yields bootstrap standard errors and percentile confidence intervals for
#' \code{Profile_Means}, plus a Wald chi-square test of "all profile means
#' equal" using the bootstrap covariance -- reported in
#' \code{$Bootstrap_Correction}, and preferable to \code{$ANOVA_Table}'s
#' p-value for publication-grade inference. It is \emph{not} the Bakk et al.
#' (2014) analytic formula; treat it as a practical approximation with the
#' same goal (each refit's arbitrary profile labels are first aligned to
#' \code{model}'s by an exact optimal assignment of standardized profile
#' means, to avoid mixing different real-world profiles together across
#' resamples). Each refit runs sequentially (\code{cores = 1}) inside its
#' bootstrap replicate, whatever \code{cores} the original model used, so
#' parallelism happens only across replicates. It is off (\code{"none"}) by default
#' because it requires \code{n_boot} additional full model refits and is
#' therefore substantially slower; use \code{cores > 1} to parallelize it.
#'
#' @param model A fitted model returned by \code{\link{robust_lpa}} (latent
#'   profiles) or \code{\link{robust_gmm}} (latent classes of trajectories).
#' @param aux_var A numeric vector of the continuous auxiliary (distal outcome)
#'   variable, of length \code{nrow(model$probabilities)}: one value per
#'   observation for \code{robust_lpa()} fits, one value per person, in the
#'   order of \code{model$ids}, for \code{robust_gmm()} fits (e.g. a
#'   baseline characteristic, or a distal outcome). \code{NA} values are
#'   dropped listwise before the BCH calculations.
#' @param correction String, either \code{"none"} (default; the fixed-D
#'   F-test only) or \code{"bootstrap"} (also compute the bootstrap
#'   classification-uncertainty correction described in Details). Ignored
#'   (with \code{$Bootstrap_Correction = NULL}) when \code{"none"}.
#' @param n_boot Integer, number of bootstrap resamples to use when
#'   \code{correction = "bootstrap"} (default \code{200}; at least \code{10}).
#'   Each resample refits the full \code{\link{robust_lpa}} model, so runtime
#'   scales linearly with \code{n_boot}. Ignored when \code{correction = "none"}.
#' @param cores Integer, number of CPU cores to use to run the \code{n_boot}
#'   bootstrap refits in parallel when \code{correction = "bootstrap"}
#'   (default \code{1}, sequential); uses the same backend as
#'   \code{\link{robust_lpa}}'s own \code{cores} argument. Ignored when
#'   \code{correction = "none"}.
#' @return A list containing:
#'   \describe{
#'     \item{Profile_Means}{Named numeric vector of BCH bias-corrected profile means of \code{aux_var}.}
#'     \item{ANOVA_Table}{A data.frame with \code{Df}, \code{Sum_Sq}, \code{Mean_Sq}, \code{F_value}, and \code{p_value} for the "Class" and "Residuals" rows (see the fixed-D caveat in Details).}
#'     \item{Classification_Matrix}{The \code{G x G} classification error matrix \strong{D}, with \code{D[t, s]} = P(assigned \code{s} | true profile \code{t}); rows sum to 1.}
#'     \item{Classification_Weights}{The \code{G x G} BCH weight matrix \eqn{W = D^{-1}} (rows: assigned profile; columns: true profile).}
#'     \item{N_Used}{Integer, the number of observations retained after removing missing \code{aux_var} values.}
#'     \item{Bootstrap_Correction}{\code{NULL} unless \code{correction = "bootstrap"}, in which case a list with \code{n_boot_used}, \code{n_boot_failed}, \code{SE} and \code{CI_lower}/\code{CI_upper} (per profile), and the overall \code{Wald_stat}/\code{Wald_df}/\code{Wald_p_value} test of equal profile means (see Details).}
#'   }
#' @references
#'   Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#'   models with categorical variables: One-step versus three-step
#'   estimators. \emph{Political Analysis}, 12(1), 3-27. \doi{10.1093/pan/mph001}
#'
#'   Vermunt, J. K. (2010). Latent class modeling with covariates: Two
#'   improved three-step approaches. \emph{Political Analysis}, 18(4), 450-469. \doi{10.1093/pan/mpq025}
#'
#'   Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#'   assignments to external variables: Standard errors for correct
#'   inference. \emph{Political Analysis}, 22(4), 520-540. \doi{10.1093/pan/mpu003}
#' @examples
#' data(neuro_data)
#' # Fit the model on Memory and RT_Stroop only
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' set.seed(1)
#' fit <- robust_lpa(data = x, G = 2, model = 1, n_starts = 3)
#' summary(fit)  # profile means for the two fitted variables
#' # Test RT_TMT (not used to fit the model) as an auxiliary outcome
#' bch_res <- bch_robust(fit, neuro_data$RT_TMT)
#' bch_res$Profile_Means
#' bch_res$ANOVA_Table
#'
#' \donttest{
#' # Add the bootstrap classification-uncertainty correction (slower: refits
#' # the model n_boot times). A small n_boot here is just for a fast demo --
#' # use several hundred for publication-grade inference.
#' bch_res_boot <- bch_robust(fit, neuro_data$RT_TMT, correction = "bootstrap", n_boot = 30)
#' bch_res_boot$Bootstrap_Correction
#' }
#' @export
bch_robust <- function(model, aux_var, correction = c("none", "bootstrap"), n_boot = 200, cores = 1) {
  correction <- match.arg(correction)

  if (is.null(model$probabilities) || is.null(model$assignments)) {
    stop("Invalid model object. Please provide a model fitted with robust_lpa() or robust_gmm().")
  }
  is_gmm <- inherits(model, "robust_gmm")
  if (length(aux_var) != nrow(model$probabilities)) {
    stop("The length of the auxiliary variable must match the number of ",
         if (is_gmm) "persons in the model (one value per element of `model$ids`)." else "observations in the model.")
  }
  if (!is.numeric(aux_var)) stop("`aux_var` must be numeric.")

  z <- model$probabilities
  assignments <- model$assignments
  G <- ncol(z)
  n <- nrow(z)

  if (G < 2) stop("BCH requires a model with at least 2 profiles.")

  pe <- .bch_point_estimate(z, assignments, aux_var)
  anova_table <- .bch_weighted_anova(pe$W_matrix, pe$assignments_valid, pe$aux_valid, G, pe$n_valid)

  boot_out <- NULL
  if (correction == "bootstrap") {
    if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 10 || n_boot != round(n_boot)) {
      stop("`n_boot` must be a single integer of at least 10.")
    }
    if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) {
      stop("`cores` must be a single positive integer.")
    }
    if (is.null(model$data) || is.null(model$call_args)) {
      stop(
        "`correction = \"bootstrap\"` requires `model$data` and `model$call_args`, which are only ",
        "present on fits from the current version of robust_lpa(); re-fit `model` and try again."
      )
    }

    orig_means <- model$means
    boot_data_full <- model$data
    call_args <- model$call_args
    # standardize variables when matching profile labels across refits
    match_scale <- if (is_gmm) NULL else sqrt(pmax(Reduce(`+`, Map(function(S, w) w * diag(S),
                                                                   model$covariances, model$proportions)), 1e-8))

    # One bootstrap replicate: resample observations (data row + aux_var
    # value together, preserving their pairing), refit the identical
    # robust_lpa() specification, align its (arbitrarily labeled) profiles to
    # `model`'s via .match_profile_labels(), then recompute the BCH point
    # estimate on the aligned refit. Wrapped in its own tryCatch so a
    # degenerate resample (e.g. a profile with zero modal assignments) drops
    # that replicate instead of aborting the whole bootstrap.
    run_one_boot <- function(b) {
      tryCatch({
        idx <- sample.int(n, size = n, replace = TRUE)
        boot_aux <- aux_var[idx]
        if (is_gmm) {
          # resample persons (with all their occasions) and refit robust_gmm()
          rf <- .gmm_refit_resample(model, idx)
          pe_b <- .bch_point_estimate(rf$z, rf$assignments, boot_aux, warn_zero_weight = FALSE)
          return(pe_b$weighted_means)
        }
        boot_data <- boot_data_full[idx, , drop = FALSE]

        boot_fit <- suppressWarnings(
          do.call(robust_lpa, modifyList(call_args, list(data = boot_data, cores = 1)))
        )

        perm <- .match_profile_labels(orig_means, boot_fit$means, scale = match_scale)
        inv_perm <- integer(G)
        inv_perm[perm] <- seq_len(G)
        z_aligned <- boot_fit$probabilities[, perm, drop = FALSE]
        assignments_aligned <- inv_perm[boot_fit$assignments]

        pe_b <- .bch_point_estimate(z_aligned, assignments_aligned, boot_aux, warn_zero_weight = FALSE)
        pe_b$weighted_means
      }, error = function(e) NULL)
    }

    boot_results <- .run_parallel(
      cores, n_boot, run_one_boot,
      export_vars = c("n", "G", "boot_data_full", "aux_var", "call_args", "orig_means", "match_scale", "is_gmm", "model"),
      export_env = environment()
    )
    boot_means_list <- Filter(Negate(is.null), boot_results)

    min_needed <- max(10, G + 2)
    if (length(boot_means_list) < min_needed) {
      warning(sprintf(
        "Only %d of %d bootstrap replicates succeeded (need at least %d); skipping the bootstrap correction (`$Bootstrap_Correction` is NULL). Increase `n_boot` or inspect why replicates are failing (e.g. near-empty profiles on resampled data).",
        length(boot_means_list), n_boot, min_needed
      ))
    } else {
      boot_means_mat <- do.call(rbind, boot_means_list)
      boot_means_mat <- boot_means_mat[stats::complete.cases(boot_means_mat), , drop = FALSE]

      if (nrow(boot_means_mat) < min_needed) {
        warning("Too many bootstrap replicates produced NA profile means (near-zero BCH weights); skipping the bootstrap correction (`$Bootstrap_Correction` is NULL).")
      } else {
        se <- apply(boot_means_mat, 2, stats::sd)
        ci <- apply(boot_means_mat, 2, stats::quantile, probs = c(0.025, 0.975))

        Sigma_hat <- stats::cov(boot_means_mat)
        Cmat <- cbind(-1, diag(G - 1))
        Ctheta <- as.numeric(Cmat %*% pe$weighted_means)
        CSigmaCt <- Cmat %*% Sigma_hat %*% t(Cmat)

        wald <- tryCatch({
          wald_stat <- as.numeric(t(Ctheta) %*% solve(CSigmaCt) %*% Ctheta)
          list(stat = wald_stat, p_value = stats::pchisq(wald_stat, df = G - 1, lower.tail = FALSE))
        }, error = function(e) {
          warning("The bootstrap covariance of the profile means is singular; the Wald test could not be computed (SEs/CIs are still reported).")
          list(stat = NA_real_, p_value = NA_real_)
        })

        boot_out <- list(
          n_boot_used = nrow(boot_means_mat),
          n_boot_failed = n_boot - nrow(boot_means_mat),
          SE = se,
          CI_lower = ci[1, ],
          CI_upper = ci[2, ],
          Wald_stat = wald$stat,
          Wald_df = G - 1,
          Wald_p_value = wald$p_value
        )
      }
    }
  }

  return(list(
    Profile_Means = pe$weighted_means,
    ANOVA_Table = anova_table,
    Classification_Matrix = pe$D_matrix,
    Classification_Weights = pe$W_matrix,
    N_Used = pe$n_valid,
    Bootstrap_Correction = boot_out
  ))
}
