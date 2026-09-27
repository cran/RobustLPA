#' Force a Matrix to be Symmetric Positive Definite
#'
#' Symmetrizes a matrix and, if its smallest eigenvalue is below a floor,
#' adds a ridge to the diagonal so that it becomes safely positive definite.
#' Used as a numerical safety net after every covariance update for models
#' 3-6 in \code{\link{robust_lpa}} (models 1-2 use a separate, simpler
#' per-variable floor -- see \code{\link{.apply_covariance_model}}).
#'
#' The floor is set \emph{relative to the matrix's own average eigenvalue}
#' (\code{min_ratio} times \code{mean(eigenvalues)}), not a fixed absolute
#' constant. An earlier version used a fixed constant (effectively ~1e-5),
#' which is a negligible floor relative to the variance scale of typical
#' (e.g. standardized, variance ~= 1) data -- too weak to stop a profile's
#' covariance from collapsing toward a near-singular, near-zero-variance
#' matrix that fits a tiny, near-coincident handful of points almost
#' perfectly. That is a classic spurious/degenerate local maximum of the
#' unconstrained Gaussian mixture likelihood (see e.g. McLachlan & Peel,
#' 2000, \emph{Finite Mixture Models}, Sec. 3.10): the resulting "profile"
#' can have both an implausibly high log-likelihood contribution and a
#' vanishingly small size, and is a common way multivariate EM for
#' unconstrained-covariance mixtures goes wrong in any implementation, not
#' specific to this package. A relative floor scales automatically with
#' whatever units \code{data} is on and meaningfully resists this collapse,
#' at the cost of (very slightly) regularizing genuinely tight-but-real
#' clusters toward less extreme covariances. This does not eliminate the
#' possibility of a spurious solution -- no finite floor can, since it is a
#' structural property of unconstrained-covariance finite mixture ML -- but
#' substantially raises how tight a cluster has to be before it can trigger
#' one. If you see an implausibly small profile size for model 3-6, compare
#' against a smaller \code{G}, more \code{n_starts}, or a more constrained
#' \code{model} (1-2, which do not go through this floor).
#'
#' @param mat A square numeric matrix.
#' @param min_ratio Numeric, the minimum eigenvalue as a fraction of
#'   \code{mat}'s own average eigenvalue. Default \code{1e-3}.
#' @return A symmetric positive-definite matrix of the same dimension.
#' @keywords internal
#' @noRd
.force_pd <- function(mat, min_ratio = 1e-3) {
  mat <- (mat + t(mat)) / 2
  ev <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  floor_ev <- max(mean(ev) * min_ratio, 1e-8)
  min_ev <- min(ev)
  if (min_ev < floor_ev) {
    mat <- mat + diag(floor_ev - min_ev, nrow(mat))
  }
  mat
}

#' Apply a Variance-Covariance Model to a Set of Raw Profile Covariances
#'
#' Given the (unconstrained) per-profile scatter matrices \eqn{S_g} produced
#' by the M-step (\code{\link{robust_m_step}}) and the mixing proportions
#' \eqn{\pi_g}, returns the constrained covariance matrices implied by
#' \code{model} (see \code{\link{robust_lpa}}) that minimize the
#' expected complete-data discrepancy
#' \deqn{F(\Sigma_1, \dots, \Sigma_G) = \sum_g \pi_g \left[\log|\Sigma_g| + \mathrm{tr}(\Sigma_g^{-1} S_g)\right],}
#' i.e. the constrained M-step.
#'
#' Models 1, 2, 3 and 6 have closed-form solutions. Models 4 and 5 do not,
#' and are solved by conditional maximization (see \code{.cm_model4()} and
#' \code{.cm_model5()}), started from the current covariances when supplied,
#' so that every call decreases \eqn{F} and the EM algorithm remains a
#' (generalized) EM with a monotone likelihood.
#'
#' @param raw_sigmas A list of length \code{G} with unconstrained covariance matrices.
#' @param pi_g Numeric vector of length \code{G} with the current mixing proportions.
#' @param model An integer between 1 and 6.
#' @param G Integer, the number of profiles.
#' @param p Integer, the number of variables.
#' @param current Optional list of the current (constrained) covariance
#'   matrices, used as the starting point for models 4 and 5.
#' @return A list of length \code{G} with the constrained covariance matrices.
#' @keywords internal
#' @noRd
.apply_covariance_model <- function(raw_sigmas, pi_g, model, G, p, current = NULL) {
  sigma <- vector("list", G)
  w <- pi_g / sum(pi_g)

  if (model == 1) {
    # Equal variances across profiles, covariances fixed to 0 (diagonal, shared
    # across profiles; not a single isotropic/spherical variance across variables).
    pooled_diag <- numeric(p)
    for (g in 1:G) pooled_diag <- pooled_diag + w[g] * diag(raw_sigmas[[g]])
    pooled_diag <- pmax(pooled_diag, 1e-8)
    for (g in 1:G) sigma[[g]] <- diag(pooled_diag, p)

  } else if (model == 2) {
    # Varying variances across profiles, covariances fixed to 0.
    for (g in 1:G) sigma[[g]] <- diag(pmax(diag(raw_sigmas[[g]]), 1e-8), p)

  } else if (model == 3) {
    # Equal variances and equal covariances across profiles (one shared,
    # full covariance matrix).
    pooled_sigma <- matrix(0, nrow = p, ncol = p)
    for (g in 1:G) pooled_sigma <- pooled_sigma + w[g] * raw_sigmas[[g]]
    pooled_sigma <- .force_pd(pooled_sigma)
    for (g in 1:G) sigma[[g]] <- pooled_sigma

  } else if (model == 4) {
    # Varying variances, shared correlation matrix: Sigma_g = D_g R D_g.
    S <- lapply(raw_sigmas, .force_pd)
    sigma <- lapply(.cm_model4(S, w, current), .force_pd)

  } else if (model == 5) {
    # Shared variances, profile-specific correlation matrices: Sigma_g = D R_g D.
    S <- lapply(raw_sigmas, .force_pd)
    sigma <- lapply(.cm_model5(S, w, current), .force_pd)

  } else if (model == 6) {
    # Fully unconstrained: each profile has its own variances and covariances.
    for (g in 1:G) sigma[[g]] <- .force_pd(raw_sigmas[[g]])

  } else {
    stop("`model` must be an integer between 1 and 6.")
  }

  sigma
}

#' Constrained-Covariance Discrepancy
#' @keywords internal
#' @noRd
.cov_discrepancy <- function(Sigma_list, S_list, w) {
  total <- 0
  for (g in seq_along(S_list)) {
    L <- tryCatch(chol(Sigma_list[[g]]), error = function(e) NULL)
    if (is.null(L)) return(Inf)
    Sinv <- chol2inv(L)
    total <- total + w[g] * (2 * sum(log(diag(L))) + sum(Sinv * S_list[[g]]))
  }
  total
}

#' Minimize -2 sum(log a) + a' M a over a > 0 (Coordinate Descent)
#'
#' The objective is strictly convex on the positive orthant; each coordinate
#' update is the positive root of \eqn{M_{jj} a_j^2 + b_j a_j - 1 = 0} with
#' \eqn{b_j = \sum_{k \ne j} M_{jk} a_k}.
#' @keywords internal
#' @noRd
.solve_inverse_scales <- function(M, a, sweeps = 50, tol = 1e-10) {
  p <- length(a)
  for (it in seq_len(sweeps)) {
    a_old <- a
    for (j in seq_len(p)) {
      b <- sum(M[j, -j] * a[-j])
      mjj <- max(M[j, j], 1e-12)
      a[j] <- (-b + sqrt(b * b + 4 * mjj)) / (2 * mjj)
    }
    if (max(abs(a - a_old) / pmax(abs(a_old), 1e-12)) < tol) break
  }
  a
}

#' Conditional-Maximization M-Step for Model 4 (Sigma_g = D_g R D_g)
#'
#' Alternates two exact conditional minimizations of the discrepancy
#' \eqn{F}: (i) given the profile scales \eqn{D_g}, the shared matrix
#' \eqn{\Psi = \sum_g \pi_g D_g^{-1} S_g D_g^{-1}}; (ii) given \eqn{\Psi},
#' each \eqn{D_g} by coordinate descent on its inverse scales. The family
#' \eqn{D_g \Psi D_g} with unconstrained \eqn{\Psi} is identical to
#' \eqn{D_g R D_g} with a correlation matrix \eqn{R} (rescale \eqn{\Psi} to
#' unit diagonal and absorb the scale into \eqn{D_g}), so no constraint is
#' lost.
#' @keywords internal
#' @noRd
.cm_model4 <- function(S, w, current = NULL, max_iter = 25, tol = 1e-9) {
  G <- length(S)
  p <- nrow(S[[1]])
  d <- lapply(seq_len(G), function(g) {
    src <- if (!is.null(current)) current[[g]] else S[[g]]
    sqrt(pmax(diag(src), 1e-8))
  })
  build <- function(d, Psi) lapply(seq_len(G), function(g) Psi * outer(d[[g]], d[[g]]))
  compute_psi <- function(d) {
    Psi <- matrix(0, p, p)
    for (g in seq_len(G)) Psi <- Psi + w[g] * (S[[g]] / outer(d[[g]], d[[g]]))
    .force_pd((Psi + t(Psi)) / 2)
  }
  Psi <- if (!is.null(current)) stats::cov2cor(current[[1]]) else compute_psi(d)
  F_old <- if (!is.null(current)) .cov_discrepancy(current, S, w) else Inf
  best <- if (!is.null(current)) current else build(d, Psi)
  F_best <- F_old

  for (it in seq_len(max_iter)) {
    Psi <- compute_psi(d)
    Psi_inv <- chol2inv(chol(Psi))
    for (g in seq_len(G)) {
      a <- .solve_inverse_scales(Psi_inv * S[[g]], 1 / d[[g]])
      d[[g]] <- 1 / a
    }
    cand <- build(d, Psi)
    F_new <- .cov_discrepancy(cand, S, w)
    if (F_new < F_best) {
      best <- cand
      F_best <- F_new
    }
    if (is.finite(F_old) && abs(F_old - F_new) < tol * max(1, abs(F_new))) break
    F_old <- F_new
  }
  best
}

#' Conditional-Maximization M-Step for Model 5 (Sigma_g = D R_g D)
#'
#' Alternates (i) an exact conditional minimization of the discrepancy
#' \eqn{F} over the shared scales \eqn{D} given the profile correlation
#' matrices (coordinate descent on the inverse scales), and (ii) for each
#' profile, a minimization over its correlation matrix \eqn{R_g} given the
#' scales (no closed form exists for a correlation matrix with known
#' variances): a step toward the correlation matrix of the standardized
#' scatter \eqn{D^{-1} S_g D^{-1}} followed by projected-gradient steps on
#' the off-diagonal entries, each with step halving. Every accepted step
#' decreases \eqn{F}, so the EM algorithm remains monotone.
#' @keywords internal
#' @noRd
.cm_model5 <- function(S, w, current = NULL, max_iter = 5, tol = 1e-9) {
  G <- length(S)
  p <- nrow(S[[1]])
  if (!is.null(current)) {
    d <- sqrt(pmax(diag(current[[1]]), 1e-8))
    R <- lapply(current, stats::cov2cor)
  } else {
    pooled <- numeric(p)
    for (g in seq_len(G)) pooled <- pooled + w[g] * diag(S[[g]])
    d <- sqrt(pmax(pooled, 1e-8))
    R <- lapply(S, stats::cov2cor)
  }
  build <- function(d, R) lapply(R, function(Rg) Rg * outer(d, d))
  F_cur <- .cov_discrepancy(build(d, R), S, w)

  for (it in seq_len(max_iter)) {
    F_start <- F_cur
    # (i) shared scales given the correlation matrices (exact)
    M <- matrix(0, p, p)
    for (g in seq_len(G)) M <- M + w[g] * (chol2inv(chol(.force_pd(R[[g]]))) * S[[g]])
    d_new <- 1 / .solve_inverse_scales(M, 1 / d)
    F_new <- .cov_discrepancy(build(d_new, R), S, w)
    if (F_new < F_cur) {
      d <- d_new
      F_cur <- F_new
    }
    # (ii) correlation matrices given the scales: move toward the
    # correlation matrix of the standardized scatter (step halving), then
    # refine with projected-gradient steps on the off-diagonal entries,
    # accepting only steps that decrease the objective.
    for (g in seq_len(G)) {
      Tg <- S[[g]] / outer(d, d)
      fg <- function(Rm) {
        L <- tryCatch(chol(Rm), error = function(e) NULL)
        if (is.null(L) || any(abs(Rm[upper.tri(Rm)]) >= 1)) return(Inf)
        2 * sum(log(diag(L))) + sum(chol2inv(L) * Tg)
      }
      f_cur <- fg(R[[g]])
      target <- stats::cov2cor(.force_pd(Tg))
      step <- 1
      while (step >= 1 / 64) {
        R_try <- (1 - step) * R[[g]] + step * target
        f_try <- fg(R_try)
        if (f_try < f_cur) {
          R[[g]] <- R_try
          f_cur <- f_try
          break
        }
        step <- step / 2
      }
      for (k in seq_len(5)) {
        Rinv <- chol2inv(chol(R[[g]]))
        grad <- Rinv - Rinv %*% Tg %*% Rinv
        diag(grad) <- 0
        if (max(abs(grad)) < 1e-8) break
        step <- 1
        improved <- FALSE
        while (step >= 1e-6) {
          R_try <- R[[g]] - step * grad
          diag(R_try) <- 1
          f_try <- fg(R_try)
          if (f_try < f_cur - 1e-12) {
            R[[g]] <- R_try
            f_cur <- f_try
            improved <- TRUE
            break
          }
          step <- step / 2
        }
        if (!improved) break
      }
    }
    F_cur <- .cov_discrepancy(build(d, R), S, w)
    if (abs(F_start - F_cur) < tol * max(1, abs(F_cur))) break
  }
  build(d, R)
}

#' Compute Fit Indices and Hard Assignments for a Fitted Mixture
#'
#' Shared helper used by both the EM and MCMC branches of
#' \code{\link{robust_lpa}} to compute AIC/BIC/SABIC, classification entropy,
#' profile sizes, and modal (hard) assignments from a fitted model's means
#' and posterior probabilities. Centralizing this logic keeps the two
#' engines' fit indices consistent by construction.
#'
#' @param mu A list of length \code{G} with the estimated profile means.
#' @param log_lik The observed-data log-likelihood of the fitted model.
#' @param z An \code{n x G} matrix of posterior profile-membership probabilities.
#' @param n Integer, the number of observations.
#' @param p Integer, the number of variables.
#' @param G Integer, the number of profiles.
#' @param model An integer between 1 and 6.
#' @param extra_params Integer, additional free parameters not implied by
#'   \code{model} (e.g. 1 for an estimated multivariate-t degrees of freedom).
#' @param penalized Logical; \code{TRUE} when the means were estimated with a
#'   LASSO penalty (EM engine, \code{lambda > 0}). Only then are the mean
#'   components shrunk exactly to zero excluded from the parameter count.
#' @return A list with \code{fit_indices} (a one-row data.frame) and
#'   \code{assignments} (an integer vector of length \code{n}).
#' @keywords internal
#' @noRd
.compute_fit_indices <- function(mu, log_lik, z, n, p, G, model, extra_params = 0, penalized = FALSE) {
  if (G == 1) {
    K_base <- switch(as.character(model),
                     "1" = p + p, "2" = p + p, "3" = p + (p * (p + 1) / 2),
                     "4" = p + (p * (p + 1) / 2), "5" = p + (p * (p + 1) / 2), "6" = p + (p * (p + 1) / 2))
  } else {
    K_base <- switch(as.character(model),
                     "1" = (G * p) + p + (G - 1), "2" = (G * p) + (G * p) + (G - 1),
                     "3" = (G * p) + (p * (p + 1) / 2) + (G - 1), "4" = (G * p) + (G * p) + (p * (p - 1) / 2) + (G - 1),
                     "5" = (G * p) + p + (G * p * (p - 1) / 2) + (G - 1), "6" = (G * p) + (G * (p * (p + 1) / 2)) + (G - 1))
  }
  
  # Means set exactly to zero by the LASSO are not free parameters. Without a
  # penalty a mean that happens to be ~0 (e.g. G = 1 on standardized data) is
  # still estimated, so it must be counted.
  zeroed_means <- if (penalized) sum(vapply(mu, function(m) sum(abs(m) < 1e-5), numeric(1))) else 0
  K_eff <- max(1, K_base - zeroed_means + extra_params)
  
  AIC_val <- -2 * log_lik + 2 * K_eff
  BIC_val <- -2 * log_lik + K_eff * log(n)
  SABIC_val <- -2 * log_lik + K_eff * log((n + 2) / 24)
  
  z_safe <- z
  z_safe[z_safe < 1e-15] <- 1e-15
  entropy <- if (G == 1) 1 else 1 - (sum(-z * log(z_safe)) / (n * log(G)))
  
  assignments <- max.col(z)
  sizes <- table(factor(assignments, levels = 1:G))
  
  list(
    fit_indices = data.frame(
      Model = model, Profiles = G, LogLik = log_lik, Parameters = K_eff,
      AIC = AIC_val, BIC = BIC_val, SABIC = SABIC_val, Entropy = entropy,
      Min_Size = min(sizes) / n, Max_Size = max(sizes) / n
    ),
    assignments = assignments
  )
}

#' Warn When the Smallest Fitted Profile Looks Implausibly Small
#'
#' A transparency check run on every \code{\link{robust_lpa}} fit: an
#' implausibly small smallest-profile size (well below what a genuine small
#' subgroup would usually look like) is the standard symptom of a
#' spurious/degenerate solution -- most often a profile's covariance having
#' collapsed toward a tiny, near-coincident handful of points (see
#' \code{\link{.force_pd}}'s documentation), which is structurally more
#' likely for models 3-6 (which allow off-diagonal/full covariance) than for
#' models 1-2 (diagonal only). This does not attempt to distinguish a
#' genuine small group from a spurious one -- that requires judgment (compare
#' across `model`, more `n_starts`, fewer profiles) -- it only flags the
#' pattern so it isn't missed, especially by code that automatically selects
#' a "best" fit by AIC/BIC (e.g. \code{\link{plot_robust_lpa}} on an
#' \code{\link{estimate_profiles_robust}} result), where a spurious solution
#' can otherwise look deceptively attractive (an artificially high
#' log-likelihood from the near-singular covariance).
#'
#' @param fit_indices The one-row fit-indices data.frame from \code{\link{.compute_fit_indices}}.
#' @param engine String, \code{"EM"} or \code{"MCMC"} (for the warning text only).
#' @return \code{NULL}, invisibly; called for its \code{warning()} side effect.
#' @keywords internal
#' @noRd
.warn_if_degenerate_profile <- function(fit_indices, engine) {
  if (isTRUE(fit_indices$Min_Size < 0.02)) {
    warning(sprintf(
      "The smallest fitted profile holds only %.1f%% of the observations (model %d, %d profile(s), %s engine). This can be a genuine small group, but is also the classic symptom of a spurious/degenerate solution (see ?robust_lpa, \".force_pd()\" in the source) -- especially for model 3-6 with few `n_starts`. Before trusting this fit, compare against a more constrained `model` (1-2), more `n_starts`, or fewer profiles.",
      100 * fit_indices$Min_Size, fit_indices$Model, fit_indices$Profiles, engine
    ), call. = FALSE)
  }
  invisible(NULL)
}