# Internal helpers for robust_gmm(): data preparation, initialization,
# constrained covariance updates, the penalized (LASSO) fixed-effects step,
# degrees of freedom, and back-transformation to the original scale.

#' Prepare Long-Format Longitudinal Data for the Growth-Mixture Engine
#'
#' Stacks the observed (time, outcome) pairs of every person, sorted by
#' person, outcome and time, standardizes every outcome (mean 0, SD 1 over
#' all its observations) and records where each person's rows start.
#'
#' @param data A data.frame in long format (one row per person and occasion).
#' @param id,time Column names of the person identifier and of the time variable.
#' @param outcomes Character vector with the column names of the outcomes.
#' @return A list with the standardized values \code{y}, \code{time},
#'   0-based outcome index \code{outc}, person row offsets \code{starts},
#'   the person identifiers \code{ids}, their number \code{N}, the number of
#'   outcomes \code{K}, the standardization constants \code{center} and
#'   \code{scale}, and the original data subset \code{raw}.
#' @keywords internal
#' @noRd
.gmm_prepare <- function(data, id, time, outcomes) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame in long format (one row per person and occasion).")
  if (!is.character(id) || length(id) != 1) stop("`id` must be a single column name.")
  if (!is.character(time) || length(time) != 1) stop("`time` must be a single column name.")
  if (!is.character(outcomes) || length(outcomes) < 1) stop("`outcomes` must be a character vector of column names.")
  if (anyDuplicated(outcomes)) stop("`outcomes` must not contain duplicates.")
  missing_cols <- setdiff(c(id, time, outcomes), names(data))
  if (length(missing_cols) > 0) stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "))
  if (!is.numeric(data[[time]])) stop("The `time` column must be numeric.")
  for (o in outcomes) if (!is.numeric(data[[o]])) stop("Outcome column '", o, "' must be numeric.")

  raw <- data[!is.na(data[[id]]), c(id, time, outcomes), drop = FALSE]
  ids_all <- unique(raw[[id]])
  K <- length(outcomes)

  pieces <- lapply(seq_len(K), function(k) {
    v <- raw[[outcomes[k]]]
    keep <- !is.na(v) & !is.na(raw[[time]])
    data.frame(pid = match(raw[[id]][keep], ids_all), t = raw[[time]][keep], k = k - 1L,
               y = v[keep])
  })
  long <- do.call(rbind, pieces)
  src <- do.call(c, lapply(seq_len(K), function(k) {
    v <- raw[[outcomes[k]]]
    which(!is.na(v) & !is.na(raw[[time]]))
  }))
  long$src_row <- src
  if (nrow(long) == 0) stop("`data` contains no observed outcome values.")

  present <- sort(unique(long$pid))
  if (length(present) < length(ids_all)) {
    warning(length(ids_all) - length(present), " person(s) with no observed outcome value were dropped.",
            call. = FALSE)
  }
  ids <- ids_all[present]
  long$pid <- match(long$pid, present)
  ord <- order(long$pid, long$k, long$t)
  long <- long[ord, , drop = FALSE]

  center <- vapply(seq_len(K), function(k) mean(long$y[long$k == k - 1L]), numeric(1))
  scale <- vapply(seq_len(K), function(k) stats::sd(long$y[long$k == k - 1L]), numeric(1))
  if (any(!is.finite(scale) | scale <= 0)) {
    stop("Every outcome needs at least two distinct observed values (outcome(s) with zero variance: ",
         paste(outcomes[!is.finite(scale) | scale <= 0], collapse = ", "), ").")
  }
  names(center) <- names(scale) <- outcomes
  N <- length(ids)
  counts <- tabulate(long$pid, nbins = N)

  list(
    y = (long$y - center[long$k + 1L]) / scale[long$k + 1L],
    y_raw = long$y, time = as.numeric(long$t), outc = as.integer(long$k),
    pid = long$pid, src_row = long$src_row, starts = as.integer(c(0L, cumsum(counts))), nobs = counts,
    ids = ids, N = N, K = K, outcomes = outcomes, center = center, scale = scale,
    id_name = id, time_name = time,
    raw = raw
  )
}

#' Names of the Polynomial Time Terms
#' @keywords internal
#' @noRd
.gmm_term_names <- function(time_name, degree) {
  c("(Intercept)", if (degree >= 1) time_name,
    if (degree >= 2) paste0(time_name, "^", 2:degree))
}

#' Person-Level Least-Squares Trajectories (Initialization Only)
#'
#' Fits each person's polynomial trajectory on every outcome by ordinary
#' least squares (standardized scale); coefficients that cannot be
#' estimated from too few occasions are imputed by the column median.
#' @keywords internal
#' @noRd
.gmm_person_ols <- function(prep, degree) {
  nb <- degree + 1
  P <- prep$K * nb
  feats <- matrix(NA_real_, prep$N, P)
  rss <- numeric(0)
  for (i in seq_len(prep$N)) {
    rows <- if (prep$starts[i + 1] > prep$starts[i]) (prep$starts[i] + 1):prep$starts[i + 1] else integer(0)
    for (k in seq_len(prep$K)) {
      rk <- rows[prep$outc[rows] == k - 1L]
      if (length(rk) == 0) next
      tt <- prep$time[rk]
      yy <- prep$y[rk]
      m <- min(nb, length(unique(tt)))
      Xp <- outer(tt, 0:(m - 1), `^`)
      fit <- tryCatch(stats::lm.fit(Xp, yy), error = function(e) NULL)
      if (is.null(fit) || anyNA(fit$coefficients)) next
      feats[i, (k - 1) * nb + seq_len(m)] <- fit$coefficients
      if (length(rk) > m) rss <- c(rss, sum(fit$residuals^2) / (length(rk) - m))
    }
  }
  for (j in seq_len(P)) {
    miss <- is.na(feats[, j])
    if (all(miss)) feats[, j] <- 0 else feats[miss, j] <- stats::median(feats[!miss, j])
  }
  list(features = feats, resid_var = if (length(rss) > 0) stats::median(rss) else 0.5)
}

#' Random-Centre Partition for a Growth-Mixture EM Start
#'
#' Draws \code{G} persons at random as class centres and assigns every person
#' to the nearest centre in the (standardized) space of person-level
#' trajectory features (Forgy initialization). Unlike random soft
#' memberships, which give \code{G} almost identical starting classes, the
#' starting classes differ from each other.
#' @param Fs Standardized person-level features (\code{N x P}).
#' @param G Number of classes.
#' @return An integer vector of class labels, or \code{NULL} if no draw out
#'   of ten gave \code{G} non-empty classes.
#' @keywords internal
#' @noRd
.gmm_random_centers <- function(Fs, G) {
  n <- nrow(Fs)
  if (n < G) return(NULL)
  for (attempt in 1:10) {
    ctr <- Fs[sample.int(n, G), , drop = FALSE]
    d2 <- matrix(vapply(seq_len(G), function(g) colSums((t(Fs) - ctr[g, ])^2), numeric(n)), n, G)
    cl <- max.col(-d2, ties.method = "first")
    if (length(unique(cl)) == G) return(cl)
  }
  NULL
}

#' Starting Values for One Growth-Mixture EM Start
#' @keywords internal
#' @noRd
.gmm_init <- function(prep, ols, G, degree, nrand, init, re_structure) {
  nb <- degree + 1
  K <- prep$K
  P <- K * nb
  Q <- K * nrand
  feats <- ols$features
  if (G == 1) {
    z0 <- matrix(1, prep$N, 1)
  } else {
    sds <- apply(feats, 2, stats::sd)
    sds[!is.finite(sds) | sds < 1e-10] <- 1
    Fs <- sweep(sweep(feats, 2, colMeans(feats), "-"), 2, sds, "/")
    cl <- NULL
    if (identical(init, "kmeans")) {
      km <- tryCatch(suppressWarnings(stats::kmeans(Fs, centers = G, nstart = 1, iter.max = 50)),
                     error = function(e) NULL)
      if (!is.null(km) && length(unique(km$cluster)) == G) cl <- km$cluster
    }
    if (is.null(cl)) cl <- .gmm_random_centers(Fs, G)
    if (!is.null(cl)) {
      z0 <- matrix(0.1 / (G - 1), prep$N, G)
      z0[cbind(seq_len(prep$N), cl)] <- 0.9
    } else {
      z0 <- matrix(stats::runif(prep$N * G), prep$N, G)
      z0 <- z0 / rowSums(z0)
    }
  }
  beta <- matrix(vapply(seq_len(G), function(g) colSums(feats * z0[, g]) / sum(z0[, g]), numeric(P)),
                 G, P, byrow = TRUE)
  sig2 <- matrix(max(ols$resid_var, 0.05), G, K)
  D <- NULL
  if (Q > 0) {
    re_cols <- as.vector(t(outer((seq_len(K) - 1) * nb, seq_len(nrand), `+`)))
    v <- apply(feats[, re_cols, drop = FALSE], 2, stats::var)
    v[!is.finite(v)] <- 0.1
    D0 <- diag(pmax(0.5 * v, 1e-3), Q)
    D <- rep(list(D0), G)
  }
  list(beta = beta, D = D, sig2 = sig2, pi = colMeans(z0))
}

#' Apply the Random-Effect Covariance Constraints
#'
#' Pools the per-profile sufficient statistics when the covariance is shared
#' across profiles, and zeroes the entries excluded by the structure
#' ("block": no correlation between outcomes; "diagonal": no correlation at
#' all). For a Gaussian random-effects distribution these are exactly the
#' constrained maximum-likelihood updates.
#' @keywords internal
#' @noRd
.gmm_update_D <- function(Dnum, Dden, D_old, re_cov, re_structure, K, nrand) {
  G <- length(Dnum)
  Q <- K * nrand
  mask <- switch(re_structure,
    full = matrix(1, Q, Q),
    block = kronecker(diag(K), matrix(1, nrand, nrand)),
    diagonal = diag(Q)
  )
  fin <- function(M) {
    M <- (M + t(M)) / 2 * mask
    .force_pd(M, min_ratio = 1e-8)
  }
  if (re_cov == "equal") {
    den <- sum(unlist(Dden))
    if (den < 1e-10) return(D_old)
    Dp <- fin(Reduce(`+`, Dnum) / den)
    return(rep(list(Dp), G))
  }
  lapply(seq_len(G), function(g) if (Dden[[g]] < 1e-10) D_old[[g]] else fin(Dnum[[g]] / Dden[[g]]))
}

#' Soft-Thresholding Operator
#' @keywords internal
#' @noRd
.soft <- function(x, thr) sign(x) * pmax(abs(x) - thr, 0)

#' Structure of the Growth-Mixture Penalties
#'
#' Index vectors used by the penalized fixed-effects step: the positions of
#' the growth terms (all polynomial terms except the intercept) in the
#' stacked coefficient vector \eqn{(\beta_1', \dots, \beta_G')'}, and the
#' outcome group of every coefficient.
#' @keywords internal
#' @noRd
.gmm_penalty_index <- function(G, K, degree) {
  nb <- degree + 1
  P <- K * nb
  term <- rep(rep(0:degree, K), G)
  outcome <- rep(rep(seq_len(K), each = nb), G)
  cls <- rep(seq_len(G), each = P)
  list(growth = which(term >= 1), outcome = outcome, class = cls, term = term, P = P, nb = nb)
}

#' Value of the Growth-Mixture Penalty
#'
#' \deqn{\lambda_{g} \sum_{g,k,l \ge 1} w^{u} |\beta_{gkl}| +
#'   \lambda_{d} \sum_{g,k,l} w^{v} |\beta_{gkl} - \bar\beta_{\cdot kl}|}
#' (element-wise), or with the difference part replaced by
#' \eqn{\lambda_d \sum_k w_k \| S_k (B_k - 1 \bar\beta_k') \|_F}
#' (group-wise by outcome, with \eqn{S_k} the information scaling of the
#' coefficients), where \eqn{\bar\beta_{\cdot kl}} is the unweighted mean
#' over profiles. Evaluated on the standardized scale.
#' @keywords internal
#' @noRd
.gmm_penalty_value <- function(beta, pen, idx) {
  if (is.null(pen) || (pen$lambda_growth <= 0 && pen$lambda_diff <= 0)) return(0)
  G <- nrow(beta)
  b <- as.vector(t(beta))
  val <- 0
  if (pen$lambda_growth > 0) val <- val + pen$lambda_growth * sum(pen$w_u * abs(b[idx$growth]))
  if (pen$lambda_diff > 0) {
    dev <- as.vector(t(sweep(beta, 2, colMeans(beta), "-")))
    if (pen$group) {
      for (k in unique(idx$outcome)) {
        sel <- idx$outcome == k
        val <- val + pen$lambda_diff * pen$w_group[k] * sqrt(sum((pen$s_v[sel] * dev[sel])^2))
      }
    } else {
      val <- val + pen$lambda_diff * sum(pen$w_v * abs(dev))
    }
  }
  val
}

#' Penalized Fixed-Effects Step by ADMM
#'
#' Minimizes, over the stacked fixed effects \eqn{\beta} of all profiles,
#' \deqn{\tfrac12 \beta' A \beta - c'\beta + N \cdot \mathrm{pen}(\beta)}
#' where \eqn{A = \mathrm{blockdiag}(A_g)} and \eqn{c} are the weighted
#' normal equations of the ECM step (log-likelihood scale) and pen is the
#' growth / difference penalty (see \code{.gmm_penalty_value()}). The
#' problem is a convex generalized Lasso, solved with the alternating
#' direction method of multipliers (Boyd et al., 2011) with splitting
#' variables \eqn{u = S\beta} (growth terms) and \eqn{v = \Delta\beta}
#' (deviations from the across-profile mean). At convergence the exact
#' zeros of \eqn{u} and \eqn{v} are imposed on \eqn{\beta}.
#'
#' @param A_list,c_list Per-profile normal equations.
#' @param beta0 Current coefficients (G x P), used as a warm start.
#' @param pen Penalty specification (see \code{.gmm_penalty_spec()}).
#' @param idx Result of \code{.gmm_penalty_index()}.
#' @param N Number of persons (the penalty is on the per-person scale).
#' @return A G x P matrix of updated coefficients.
#' @keywords internal
#' @noRd
.gmm_admm <- function(A_list, c_list, beta0, pen, idx, N, max_iter = 10000, tol = 1e-10) {
  G <- length(A_list)
  P <- idx$P
  n <- G * P
  A <- matrix(0, n, n)
  for (g in seq_len(G)) {
    sel <- (g - 1) * P + seq_len(P)
    A[sel, sel] <- A_list[[g]]
  }
  cvec <- unlist(c_list)
  use_u <- pen$lambda_growth > 0
  use_v <- pen$lambda_diff > 0

  S <- if (use_u) diag(n)[idx$growth, , drop = FALSE] else matrix(0, 0, n)
  Dm <- if (use_v) diag(n) - kronecker(matrix(1 / G, G, G), diag(P)) else matrix(0, 0, n)
  if (use_v && pen$group) Dm <- pen$s_v * Dm
  rho <- max(mean(diag(A)), 1e-6)
  H <- A + rho * (crossprod(S) + crossprod(Dm)) + diag(1e-10 * rho, n)
  Hc <- chol(H)

  beta <- as.vector(t(beta0))
  u <- as.vector(S %*% beta)
  v <- as.vector(Dm %*% beta)
  yu <- numeric(length(u))
  yv <- numeric(length(v))
  thr_u <- N * pen$lambda_growth * pen$w_u / rho

  prox_v <- function(x) {
    if (!use_v) return(x)
    if (pen$group) {
      out <- x
      for (k in unique(idx$outcome)) {
        sel <- idx$outcome == k
        nrm <- sqrt(sum(x[sel]^2))
        thr <- N * pen$lambda_diff * pen$w_group[k] / rho
        out[sel] <- if (nrm <= thr) 0 else (1 - thr / nrm) * x[sel]
      }
      out
    } else {
      .soft(x, N * pen$lambda_diff * pen$w_v / rho)
    }
  }

  for (it in seq_len(max_iter)) {
    rhs <- cvec
    if (use_u) rhs <- rhs + rho * as.vector(crossprod(S, u - yu))
    if (use_v) rhs <- rhs + rho * as.vector(crossprod(Dm, v - yv))
    beta <- backsolve(Hc, forwardsolve(t(Hc), rhs))
    Sb <- as.vector(S %*% beta)
    Db <- as.vector(Dm %*% beta)
    u_new <- if (use_u) .soft(Sb + yu, thr_u) else u
    v_new <- prox_v(Db + yv)
    r_pri <- sqrt(sum((Sb - u_new)^2) + sum((Db - v_new)^2))
    r_dual <- rho * sqrt(sum((u_new - u)^2) + sum((v_new - v)^2))
    yu <- yu + Sb - u_new
    yv <- yv + Db - v_new
    u <- u_new
    v <- v_new
    scale_ref <- max(1, sqrt(sum(beta^2)))
    if (r_pri < tol * scale_ref && r_dual < tol * scale_ref) break
  }

  B <- matrix(beta, G, P, byrow = TRUE)
  # impose the exact sparsity pattern found by the splitting variables:
  # profiles whose deviation is exactly zero sit at the across-profile mean,
  # which then equals the mean of the remaining profiles
  if (use_v) {
    Vm <- matrix(v, G, P, byrow = TRUE)
    for (p in seq_len(P)) {
      zero <- Vm[, p] == 0
      if (all(zero)) {
        B[, p] <- mean(B[, p])
      } else if (any(zero)) {
        B[zero, p] <- mean(B[!zero, p])
      }
    }
  }
  if (use_u) {
    Um <- rep(NA_real_, n)
    Um[idx$growth] <- u
    Um <- matrix(Um, G, P, byrow = TRUE)
    B[!is.na(Um) & Um == 0] <- 0
  }
  B
}

#' Penalty Specification for robust_gmm()
#'
#' @param lambda_growth,lambda_diff Non-negative penalty levels.
#' @param group Logical, group-wise difference penalty.
#' @param idx Result of \code{.gmm_penalty_index()}.
#' @param constrain Optional list with logical vectors \code{growth_zero}
#'   (over growth positions) and \code{diff_zero} (over all positions) used
#'   by the relaxed refit: those coefficients / deviations are fixed at 0,
#'   all the others are unpenalized.
#' @param beta_init Optional unpenalized estimates (G x P, standardized
#'   scale) defining adaptive-Lasso weights (Zou, 2006): the reciprocal of
#'   the absolute unpenalized coefficient (growth penalty), of the absolute
#'   unpenalized deviation from the across-class mean (element-wise
#'   difference penalty), or of the norm of an outcome's unpenalized
#'   deviations (group penalty; Wang & Leng, 2008).
#' @keywords internal
#' @noRd
.gmm_penalty_spec <- function(lambda_growth, lambda_diff, group, idx, G, K, degree,
                              constrain = NULL, beta_init = NULL, info_scale = NULL) {
  nb <- degree + 1
  if (!is.null(constrain)) {
    big <- 1e8
    return(list(
      lambda_growth = if (any(constrain$growth_zero)) big else 0,
      lambda_diff = if (any(constrain$diff_zero)) big else 0,
      group = FALSE, w_u = as.numeric(constrain$growth_zero), w_v = as.numeric(constrain$diff_zero),
      w_group = rep(1, K), s_v = rep(1, G * K * nb), constraint = TRUE, adaptive = FALSE
    ))
  }
  w_u <- rep(1, length(idx$growth))
  w_v <- rep(1, G * K * nb)
  # group penalty on information-scaled deviations (Simon & Tibshirani, 2012)
  s_v <- if (!is.null(info_scale)) as.vector(t(info_scale)) else rep(1, G * K * nb)
  s_v <- s_v / exp(mean(log(pmax(s_v, 1e-12))))
  w_group <- rep(sqrt(G * nb), K)
  if (!is.null(beta_init)) {
    b <- as.vector(t(beta_init))
    dev <- as.vector(t(sweep(beta_init, 2, colMeans(beta_init), "-")))
    w_u <- 1 / pmax(abs(b[idx$growth]), 1e-6)
    w_v <- 1 / pmax(abs(dev), 1e-6)
    w_group <- vapply(seq_len(K), function(k) {
      sel <- idx$outcome == k
      1 / max(sqrt(sum((s_v[sel] * dev[sel])^2)), 1e-6)
    }, numeric(1))
  }
  list(lambda_growth = lambda_growth, lambda_diff = lambda_diff, group = group,
       w_u = w_u, w_v = w_v, w_group = w_group, s_v = s_v, constraint = FALSE,
       adaptive = !is.null(beta_init))
}

#' Information Scale of the Fixed Effects
#'
#' Square roots of the diagonal of the marginal (GLS) information matrix of
#' every class's fixed effects at a fitted solution (G x P, standardized
#' scale); used to standardize the group difference penalty so that each
#' outcome group is penalized on the scale of its Wald statistic.
#' @keywords internal
#' @noRd
.gmm_info_scale <- function(prep, spec, res) {
  G <- spec$G
  K <- prep$K
  nrand <- spec$nrand
  lf <- .gmm_logf(prep, res$beta, res$D, res$sig2, res$pi, spec$method, res$nu, spec$degree, nrand)
  z <- .posterior_from_logf(lf$logf)$z
  out <- matrix(1, G, K * (spec$degree + 1))
  for (g in seq_len(G)) {
    w <- .gmm_weights(lf$maha[, g], lf$nobs, spec$method, spec$alpha, res$nu %||% 4)
    Dg <- if (nrand > 0) res$D[[g]] else matrix(0, 1, 1)
    gl <- gmm_class_gls_cpp(prep$y, prep$time, prep$outc, prep$starts, Dg, as.numeric(res$sig2[g, ]),
                            K, spec$degree, nrand, z[, g], w)
    out[g, ] <- sqrt(pmax(diag(gl$A), 1e-12))
  }
  out
}

#' Laplace Rates of the Bayesian (Adaptive) Lasso for the MCMC Engine
#' @keywords internal
#' @noRd
.gmm_mcmc_rates <- function(pen, idx, G, K, N) {
  P <- idx$P
  rg <- numeric(G * P)
  if (pen$lambda_growth > 0) rg[idx$growth] <- N * pen$lambda_growth * pen$w_u
  rd <- numeric(G * P)
  rk <- numeric(K)
  if (pen$lambda_diff > 0) {
    if (pen$group) rk <- N * pen$lambda_diff * pen$w_group else rd <- N * pen$lambda_diff * pen$w_v
  }
  list(growth = matrix(rg, G, P, byrow = TRUE), diff = matrix(rd, G, P, byrow = TRUE), group = rk,
       group_scale = matrix(if (isTRUE(pen$group)) pen$s_v else rep(1, G * P), G, P, byrow = TRUE))
}

#' Effective Number of Fixed-Effect Parameters
#'
#' For each coefficient position (outcome x term), the number of free
#' parameters left by the constraints that the penalties set exactly:
#' \eqn{G} minus the rank of the constraints "coefficient = 0" (growth
#' penalty) and "coefficient = across-profile mean" (difference penalty) --
#' the degrees of freedom of the generalized Lasso (Tibshirani & Taylor,
#' 2011). Without active penalties this is \code{G * K * (degree + 1)}.
#' @keywords internal
#' @noRd
.gmm_beta_df <- function(beta, idx, use_growth, use_diff, tol = 1e-8) {
  G <- nrow(beta)
  P <- ncol(beta)
  if (!use_growth && !use_diff) return(G * P)
  growth_pos <- unique(((idx$growth - 1) %% P) + 1)
  df <- 0
  for (p in seq_len(P)) {
    rows <- list()
    if (use_growth && p %in% growth_pos) {
      for (g in which(abs(beta[, p]) < tol)) { e <- numeric(G); e[g] <- 1; rows[[length(rows) + 1]] <- e }
    }
    if (use_diff && G > 1) {
      m <- mean(beta[, p])
      for (g in which(abs(beta[, p] - m) < tol * max(1, abs(m)))) {
        e <- rep(-1 / G, G); e[g] <- e[g] + 1; rows[[length(rows) + 1]] <- e
      }
    }
    rk <- if (length(rows) == 0) 0 else qr(do.call(rbind, rows), tol = 1e-7)$rank
    df <- df + G - rk
  }
  df
}

#' Back-Transform Standardized Growth-Mixture Parameters
#' @keywords internal
#' @noRd
.gmm_backtransform <- function(beta, D, sig2, prep, degree, nrand) {
  nb <- degree + 1
  K <- prep$K
  G <- nrow(beta)
  sc_beta <- rep(prep$scale, each = nb)
  shift <- rep(0, K * nb)
  shift[(seq_len(K) - 1) * nb + 1] <- prep$center
  beta_o <- sweep(sweep(beta, 2, sc_beta, "*"), 2, shift, "+")
  sig2_o <- sweep(sig2, 2, prep$scale^2, "*")
  D_o <- NULL
  if (nrand > 0) {
    sq <- rep(prep$scale, each = nrand)
    D_o <- lapply(D, function(M) M * outer(sq, sq))
  }
  list(beta = beta_o, D = D_o, sig2 = sig2_o)
}

#' Per-Profile Log-Densities of All Persons Under a Growth-Mixture Model
#' @keywords internal
#' @noRd
.gmm_logf <- function(prep, beta, D, sig2, proportions, method, nu, degree, nrand) {
  G <- nrow(beta)
  dist <- if (identical(method, "t")) 1L else 0L
  nu_val <- if (dist == 1L) nu else 1
  Q <- prep$K * nrand
  logf <- matrix(0, prep$N, G)
  maha <- matrix(0, prep$N, G)
  logdet <- matrix(0, prep$N, G)
  nobs <- NULL
  for (g in seq_len(G)) {
    Dg <- if (Q > 0) D[[g]] else matrix(0, 1, 1)
    es <- gmm_class_estep_cpp(prep$y, prep$time, prep$outc, prep$starts, as.numeric(beta[g, ]), Dg,
                              as.numeric(sig2[g, ]), prep$K, degree, nrand, dist, nu_val)
    logf[, g] <- es$logdens + log(max(proportions[g], 1e-300))
    maha[, g] <- es$maha
    logdet[, g] <- es$logdet
    nobs <- es$nobs
  }
  list(logf = logf, maha = maha, logdet = logdet, nobs = nobs)
}

#' Observed-Data Log-Likelihood of a t Growth Mixture as a Function of nu
#'
#' Uses the Mahalanobis distances and log-determinants (which do not depend
#' on nu) so that nu can be optimized at negligible cost (ECME step).
#' @keywords internal
#' @noRd
.gmm_loglik_nu <- function(nu, maha, logdet, nobs, proportions) {
  n <- matrix(nobs, nrow(maha), ncol(maha))
  lg <- lgamma((nu + n) / 2) - lgamma(nu / 2) - n / 2 * log(nu * pi) - 0.5 * logdet -
    (nu + n) / 2 * log1p(maha / nu)
  lg <- sweep(lg, 2, log(pmax(proportions, 1e-300)), "+")
  sum(.row_logsumexp(lg))
}

#' Robustness Weights for Persons (Growth Mixture)
#' @keywords internal
#' @noRd
.gmm_weights <- function(maha, nobs, method, alpha, nu) {
  .robust_weights(maha, nobs, method, alpha, nu)
}
