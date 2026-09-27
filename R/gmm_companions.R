#' Log-Likelihood of New Persons Under a Fitted Growth Mixture Model
#'
#' Observed-data log-likelihood (original scale) of the persons in
#' \code{newdata} under the parameters of \code{fit}. Used for
#' cross-validated penalty selection in \code{\link{estimate_gmm_robust}}.
#' @keywords internal
#' @noRd
.gmm_loglik_newdata <- function(fit, newdata) {
  sp <- fit$spec
  prep <- .gmm_prepare(newdata, sp$id, sp$time, sp$outcomes)
  nb <- sp$degree + 1
  K <- prep$K
  G <- sp$G
  beta_o <- matrix(vapply(fit$coefficients, function(m) as.vector(t(m)), numeric(K * nb)),
                   G, K * nb, byrow = TRUE)
  shift <- rep(0, K * nb)
  shift[(seq_len(K) - 1) * nb + 1] <- prep$center
  beta <- sweep(sweep(beta_o, 2, shift, "-"), 2, rep(prep$scale, each = nb), "/")
  sig2 <- sweep(fit$residual_var, 2, prep$scale^2, "/")
  D <- NULL
  if (sp$nrand > 0) {
    sq <- rep(prep$scale, each = sp$nrand)
    D <- lapply(fit$random_cov, function(M) unname(M) / outer(sq, sq))
  }
  lf <- .gmm_logf(prep, beta, D, sig2, fit$proportions, fit$robust_method, fit$nu, sp$degree, sp$nrand)
  .posterior_from_logf(lf$logf)$loglik - sum(log(prep$scale[prep$outc + 1L]))
}

#' Simulate Longitudinal Data From a Fitted Growth Mixture Model
#'
#' Simulates new outcome values for the persons, occasions and outcomes
#' observed in the data of \code{fit} (so the design, the visit schedule and
#' the missingness pattern are reproduced exactly), from the fitted
#' Gaussian or multivariate-t growth mixture.
#' @return A data.frame with the same columns as the fitted data.
#' @keywords internal
#' @noRd
.gmm_simulate <- function(fit) {
  sp <- fit$spec
  prep <- .gmm_prepare(fit$data, sp$id, sp$time, sp$outcomes)
  nb <- sp$degree + 1
  nrand <- sp$nrand
  K <- prep$K
  Q <- K * nrand
  it <- fit$internal
  G <- sp$G
  cls <- if (G == 1) rep(1L, prep$N) else sample.int(G, prep$N, replace = TRUE, prob = it$pi)
  u <- if (identical(fit$robust_method, "t")) stats::rchisq(prep$N, df = fit$nu) / fit$nu else rep(1, prep$N)
  B <- matrix(0, prep$N, max(Q, 1))
  if (Q > 0) {
    for (g in seq_len(G)) {
      idx <- which(cls == g)
      if (length(idx) == 0) next
      L <- chol(.force_pd(it$D[[g]], min_ratio = 1e-10))
      B[idx, ] <- (matrix(stats::rnorm(length(idx) * Q), length(idx), Q) %*% L) / sqrt(u[idx])
    }
  }
  i <- prep$pid
  k <- prep$outc + 1L
  g <- cls[i]
  tt <- prep$time
  n_obs <- length(tt)
  coef_mat <- matrix(0, n_obs, nb)
  for (l in seq_len(nb)) coef_mat[, l] <- it$beta[cbind(g, (k - 1) * nb + l)]
  fixed <- rowSums(outer(tt, 0:sp$degree, `^`) * coef_mat)
  rand <- 0
  if (Q > 0) {
    re_mat <- matrix(0, n_obs, nrand)
    for (l in seq_len(nrand)) re_mat[, l] <- B[cbind(i, (k - 1) * nrand + l)]
    rand <- rowSums(outer(tt, 0:(nrand - 1), `^`) * re_mat)
  }
  e <- stats::rnorm(length(tt), 0, sqrt(it$sig2[cbind(g, k)] / u[i]))
  y_std <- fixed + rand + e
  y_raw <- prep$center[k] + prep$scale[k] * y_std
  out <- fit$data
  for (kk in seq_len(K)) {
    sel <- k == kk
    out[prep$src_row[sel], sp$outcomes[kk]] <- y_raw[sel]
  }
  out
}

#' Refit a Growth Mixture Model on a Bootstrap Resample of Persons
#'
#' Persons are resampled with replacement (each copy receives a new
#' identifier), the model is refitted with the stored specification
#' (sequentially), and the classes are aligned to those of \code{fit} by an
#' optimal assignment of the class trajectories (on \code{fit}'s
#' standardized scale).
#' @return A list with the aligned posterior probabilities and modal classes
#'   of the resampled persons.
#' @keywords internal
#' @noRd
.gmm_refit_resample <- function(fit, idx) {
  sp <- fit$spec
  dat <- fit$data
  rows <- split(seq_len(nrow(dat)), factor(dat[[sp$id]], levels = fit$ids))
  new_rows <- unlist(rows[idx], use.names = FALSE)
  new_id <- rep(seq_along(idx), vapply(rows[idx], length, integer(1)))
  boot <- dat[new_rows, , drop = FALSE]
  boot[[sp$id]] <- new_id
  args <- utils::modifyList(fit$call_args, list(data = boot, cores = 1))
  refit <- suppressWarnings(do.call(robust_gmm, args))
  nb <- sp$degree + 1
  std_coef <- function(f) {
    b <- matrix(vapply(f$coefficients, function(m) as.vector(t(m)), numeric(length(sp$outcomes) * nb)),
                sp$G, length(sp$outcomes) * nb, byrow = TRUE)
    shift <- rep(0, ncol(b))
    shift[(seq_along(sp$outcomes) - 1) * nb + 1] <- fit$internal$center
    sweep(sweep(b, 2, shift, "-"), 2, rep(fit$internal$scale, each = nb), "/")
  }
  perm <- .gmm_match_classes(std_coef(fit), std_coef(refit), fit$internal$gram, sp$degree)
  inv_perm <- integer(sp$G)
  inv_perm[perm] <- seq_len(sp$G)
  ord <- match(seq_along(idx), refit$ids)
  list(z = refit$probabilities[ord, perm, drop = FALSE], assignments = inv_perm[refit$assignments[ord]])
}

#' Estimate Robust Growth Mixture Models Across Numbers of Classes
#'
#' Fits \code{\link{robust_gmm}} for every number of classes in
#' \code{n_classes} (and every random-effect specification in
#' \code{random}), optionally selecting the LASSO penalties of each model
#' by BIC or by K-fold cross-validation over persons, and collects the fit
#' indices in one table.
#'
#' @details
#' With \code{tune_penalty = "bic"} or \code{"cv"}, the penalty levels are
#' chosen from \code{z_grid}: each value \eqn{z} gives \eqn{\lambda = z^2/N}
#' (\eqn{N} = number of persons). With the default adaptive weights of
#' \code{\link{robust_gmm}}, a growth term (or class difference) is then set
#' to zero approximately when its Wald statistic is below \eqn{z} in
#' absolute value, which makes the grid directly interpretable (e.g.
#' \eqn{z = 2} or \eqn{z = \sqrt{\log N}}). \code{"bic"} fits every candidate
#' on the full data and keeps the lowest BIC (whose parameter count accounts
#' for the zeros and fusions); \code{"cv"} keeps the candidate with the
#' highest held-out log-likelihood over \code{k_folds} folds of persons and
#' refits it on the full data.
#'
#' @param data,id,time,outcomes See \code{\link{robust_gmm}}.
#' @param n_classes Integer vector of numbers of classes to fit.
#' @param random Character vector of random-effect specifications to fit
#'   (see \code{\link{robust_gmm}}).
#' @param cores Number of cores used across the grid of models (results are
#'   identical for any value given the same \code{set.seed()}).
#' @param tune_penalty \code{"none"} (default; penalties, if any, are taken
#'   from \code{...}), \code{"bic"} or \code{"cv"}.
#' @param penalty Which penalties to tune: \code{"both"} (default),
#'   \code{"growth"} or \code{"diff"}.
#' @param z_grid Candidate thresholds (see Details); \code{0} (no penalty)
#'   is always included.
#' @param k_folds Number of folds for \code{tune_penalty = "cv"}.
#' @param ... Further arguments passed to \code{\link{robust_gmm}} (e.g.
#'   \code{degree}, \code{robust_method}, \code{group_diff}, \code{relax},
#'   \code{n_starts}, \code{engine}).
#' @return A list with \code{fit_table} (one row per model), \code{models}
#'   (named \code{"G<g>_<random>"}) and, when penalties are tuned,
#'   \code{tuning} (the criterion of every candidate).
#' @seealso \code{\link{robust_gmm}}, \code{\link{blrt_gmm_robust}}
#' @examples
#' data(neuro_long)
#' set.seed(1)
#' res <- estimate_gmm_robust(neuro_long, id = "ID", time = "Year",
#'                            outcomes = "Memory", n_classes = 1:2, n_starts = 2)
#' res$fit_table
#' @export
estimate_gmm_robust <- function(data, id, time, outcomes, n_classes = 1:3, random = "slope",
                                cores = 1, tune_penalty = c("none", "bic", "cv"),
                                penalty = c("both", "growth", "diff"),
                                z_grid = c(1.5, 2, 2.5, 3, 4), k_folds = 5, ...) {
  tune_penalty <- match.arg(tune_penalty)
  penalty <- match.arg(penalty)
  if (!is.numeric(n_classes) || any(n_classes < 1) || any(n_classes != round(n_classes))) {
    stop("`n_classes` must contain positive integers.")
  }
  random <- unname(vapply(random, function(r) match.arg(r, c("slope", "intercept", "none")), character(1)))
  dots <- list(...)
  if (tune_penalty != "none" && any(c("lambda_growth", "lambda_diff") %in% names(dots))) {
    stop("Do not pass `lambda_growth` / `lambda_diff` when `tune_penalty` is used.")
  }
  ids <- unique(data[[id]][!is.na(data[[id]])])
  N <- length(ids)
  grid <- expand.grid(G = n_classes, random = random, stringsAsFactors = FALSE)
  zs <- sort(unique(c(0, z_grid)))
  cand <- data.frame(z = zs, lambda = zs^2 / N)

  # candidate penalty: lambda = z^2 / (number of persons in the data fitted)
  fit_one <- function(dat, G, rnd, z) {
    lam <- z^2 / length(unique(dat[[id]][!is.na(dat[[id]])]))
    lg <- if (penalty %in% c("both", "growth")) lam else 0
    ld <- if (penalty %in% c("both", "diff") && G > 1) lam else 0
    do.call(robust_gmm, c(list(data = dat, id = id, time = time, outcomes = outcomes, G = G, random = rnd,
                               lambda_growth = lg, lambda_diff = ld), dots))
  }

  worker <- function(row) {
    G <- grid$G[row]
    rnd <- grid$random[row]
    name <- paste0("G", G, "_", rnd)
    if (tune_penalty == "none") {
      fit <- tryCatch(do.call(robust_gmm, c(list(data = data, id = id, time = time, outcomes = outcomes,
                                                   G = G, random = rnd), dots)), error = function(e) e)
      if (inherits(fit, "error")) return(list(name = name, error = conditionMessage(fit)))
      return(list(name = name, fit = fit, tuning = NULL))
    }
    crit <- rep(NA_real_, nrow(cand))
    fits <- vector("list", nrow(cand))
    if (tune_penalty == "bic") {
      for (j in seq_len(nrow(cand))) {
        f <- tryCatch(suppressWarnings(fit_one(data, G, rnd, cand$z[j])), error = function(e) NULL)
        if (!is.null(f)) {
          fits[[j]] <- f
          crit[j] <- f$fit$BIC
        }
      }
      if (all(is.na(crit))) return(list(name = name, error = "all candidate fits failed"))
      best <- which.min(crit)
      fit <- fits[[best]]
    } else {
      folds <- sample(rep(seq_len(k_folds), length.out = N))
      fold_of <- folds[match(data[[id]], ids)]
      for (j in seq_len(nrow(cand))) {
        ll <- rep(NA_real_, k_folds)
        for (k in seq_len(k_folds)) {
          train <- data[!is.na(fold_of) & fold_of != k, , drop = FALSE]
          test <- data[!is.na(fold_of) & fold_of == k, , drop = FALSE]
          f <- tryCatch(suppressWarnings(fit_one(train, G, rnd, cand$z[j])), error = function(e) NULL)
          if (!is.null(f)) ll[k] <- tryCatch(.gmm_loglik_newdata(f, test), error = function(e) NA_real_)
        }
        crit[j] <- mean(ll)
      }
      if (all(is.na(crit))) return(list(name = name, error = "all candidate fits failed"))
      best <- which.max(crit)
      fit <- suppressWarnings(fit_one(data, G, rnd, cand$z[best]))
    }
    list(name = name, fit = fit,
         tuning = data.frame(G = G, random = rnd, z = cand$z, lambda = cand$lambda,
                             criterion = crit, selected = seq_len(nrow(cand)) == best))
  }

  results <- .run_parallel(cores, nrow(grid), worker, export_vars = character(0), export_env = environment())
  fit_table <- data.frame()
  models <- list()
  tuning <- data.frame()
  for (r in results) {
    if (!is.null(r$error)) {
      message(sprintf("Model %s failed: %s", r$name, r$error))
      next
    }
    row <- cbind(data.frame(Model = r$name, Random = r$fit$spec$random), r$fit$fit)
    fit_table <- rbind(fit_table, row)
    models[[r$name]] <- r$fit
    if (!is.null(r$tuning)) tuning <- rbind(tuning, r$tuning)
  }
  rownames(fit_table) <- NULL
  out <- list(fit_table = fit_table, models = models)
  if (tune_penalty != "none") {
    names(tuning)[names(tuning) == "criterion"] <- if (tune_penalty == "bic") "BIC" else "CV_LogLik"
    out$tuning <- tuning
  }
  out
}

#' Bootstrapped Likelihood Ratio Test for Robust Growth Mixture Models
#'
#' Tests a growth mixture model with \code{G} classes against the same model
#' with \code{G - 1} classes by parametric bootstrap (McLachlan, 1987;
#' Nylund, Asparouhov & Muthen, 2007): data are simulated from the fitted
#' \code{G - 1}-class model -- for the persons, occasions and outcomes
#' actually observed, so the visit schedule and the missingness pattern are
#' reproduced exactly -- and both models are refitted to every simulated
#' dataset. With \code{robust_method = "t"} the data are simulated from the
#' fitted multivariate-t model and the test compares two proper t-mixture
#' likelihoods; with the Huber estimator the statistic is a heuristic.
#'
#' @param data,id,time,outcomes See \code{\link{robust_gmm}}.
#' @param G Number of classes of the alternative model (at least 2).
#' @param n_samples Number of bootstrap samples (default 50; use 200 or more
#'   for publication).
#' @param n_starts Number of EM starts for every fit.
#' @param cores Number of cores across bootstrap samples.
#' @param ... Further arguments passed to \code{\link{robust_gmm}} (the
#'   engine is always EM and penalties are not allowed).
#' @return A list with \code{LRT_Observed}, \code{Bootstrap_LRTs},
#'   \code{p_value} (with the "+1" correction) and \code{Bootstrap_Failures}.
#' @references
#'   McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
#'   statistic for the number of components in a normal mixture.
#'   \emph{Journal of the Royal Statistical Society: Series C}, 36(3),
#'   318-324. \doi{10.2307/2347790}
#'
#'   Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
#'   number of classes in latent class analysis and growth mixture modeling:
#'   A Monte Carlo simulation study. \emph{Structural Equation Modeling},
#'   14(4), 535-569. \doi{10.1080/10705510701575396}
#' @seealso \code{\link{robust_gmm}}, \code{\link{estimate_gmm_robust}}
#' @examples
#' \donttest{
#' data(neuro_long)
#' set.seed(1)
#' blrt_gmm_robust(neuro_long, id = "ID", time = "Year", outcomes = "Memory",
#'                 G = 2, n_samples = 5, n_starts = 2)
#' }
#' @export
blrt_gmm_robust <- function(data, id, time, outcomes, G, n_samples = 50, n_starts = 3, cores = 1, ...) {
  if (!is.numeric(G) || length(G) != 1 || G < 2 || G != round(G)) stop("`G` must be a single integer of at least 2.")
  if (!is.numeric(n_samples) || n_samples < 1) stop("`n_samples` must be at least 1.")
  dots <- list(...)
  if (any(c("lambda_growth", "lambda_diff", "engine") %in% names(dots))) {
    stop("`blrt_gmm_robust()` uses the unpenalized EM engine; do not pass `lambda_*` or `engine`.")
  }
  fit_g <- function(dat, g) {
    do.call(robust_gmm, c(list(data = dat, id = id, time = time, outcomes = outcomes, G = g,
                               n_starts = n_starts, cores = 1), dots))
  }
  message(sprintf("Robust GMM BLRT: %d vs %d classes...", G - 1, G))
  mod_null <- fit_g(data, G - 1)
  mod_alt <- fit_g(data, G)
  lrt_obs <- max(0, 2 * (mod_alt$fit$LogLik - mod_null$fit$LogLik))

  one <- function(b) {
    tryCatch({
      sim <- .gmm_simulate(mod_null)
      f0 <- suppressWarnings(fit_g(sim, G - 1))
      f1 <- suppressWarnings(fit_g(sim, G))
      max(0, 2 * (f1$fit$LogLik - f0$fit$LogLik))
    }, error = function(e) NA_real_)
  }
  lrts <- unlist(.run_parallel(cores, n_samples, one, export_vars = character(0), export_env = environment()))
  n_failed <- sum(is.na(lrts))
  lrts <- lrts[!is.na(lrts)]
  if (length(lrts) == 0) stop("All bootstrap samples failed.")
  if (n_failed > 0) warning(sprintf("%d of %d bootstrap samples failed and were excluded.", n_failed, n_samples))
  list(LRT_Observed = lrt_obs, Bootstrap_LRTs = lrts,
       p_value = (sum(lrts >= lrt_obs) + 1) / (length(lrts) + 1), Bootstrap_Failures = n_failed)
}
