#' Print a Fitted Robust LPA Model
#'
#' A short, four-line-or-fewer overview of a model fitted by
#' \code{\link{robust_lpa}}: engine/model/profiles/N, the headline fit indices
#' (log-likelihood, AIC, BIC, entropy), the mixing proportions, and, for the
#' MCMC engine, the chain configuration. It deliberately omits profile means
#' and full diagnostics -- use \code{\link{summary.robust_lpa}} for those.
#'
#' @param x A \code{robust_lpa} object, as returned by \code{\link{robust_lpa}}.
#' @param ... Currently ignored (present for S3 consistency with the generic
#'   \code{\link[base]{print}}).
#' @return \code{x}, invisibly.
#' @seealso \code{\link{robust_lpa}}, \code{\link{summary.robust_lpa}}
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' fit <- suppressWarnings(robust_lpa(x, G = 2, model = 2, n_starts = 3, max_iter = 20))
#' print(fit)
#' @export
print.robust_lpa <- function(x, ...) {
  cat(sprintf(
    "<robust_lpa> %s | model %d | G = %d | N = %d\n",
    x$engine, x$fit$Model, x$fit$Profiles, nrow(x$probabilities)
  ))
  cat(sprintf(
    "LogLik = %.1f | AIC = %.1f | BIC = %.1f | Entropy = %.3f\n",
    x$fit$LogLik, x$fit$AIC, x$fit$BIC, x$fit$Entropy
  ))
  props <- round(as.numeric(x$proportions), 3)
  cat("Proportions: ", paste(sprintf("P%d=%.2f", seq_along(props), props), collapse = ", "), "\n", sep = "")
  if (identical(x$engine, "MCMC") && !is.null(x$mcmc_draws)) {
    cat(sprintf(
      "MCMC: %d chain(s) x %d iter (see summary() for Rhat/ESS)\n",
      x$mcmc_draws$n_chains, x$mcmc_draws$mcmc_iter
    ))
  }
  invisible(x)
}

#' Summarize a Fitted Robust LPA Model
#'
#' Builds a compact summary of a model fitted by \code{\link{robust_lpa}},
#' limited to the information most people actually need to interpret a fit:
#' per-profile means, profile sizes/mixing proportions, the headline fit
#' indices, and, for the MCMC engine, the Gelman-Rubin \eqn{\hat{R}} /
#' effective sample size convergence range. Returns an object of class
#' \code{"summary.robust_lpa"} with its own \code{print} method
#' (\code{\link{print.summary.robust_lpa}}), following the usual
#' \code{summary()}/\code{print(summary())} convention used throughout R (e.g.
#' \code{summary.lm}). For the full per-parameter Rhat/ESS table, use
#' \code{object$mcmc_diagnostics} directly.
#'
#' @param object A \code{robust_lpa} object, as returned by \code{\link{robust_lpa}}.
#' @param ... Currently ignored (present for S3 consistency with the generic
#'   \code{\link[base]{summary}}).
#' @return An object of class \code{"summary.robust_lpa"}, a list with
#'   \code{engine}, \code{model}, \code{G}, \code{n}, \code{means} (a
#'   variables x profiles matrix), \code{sizes} (a data.frame of profile
#'   sizes/mixing proportions), \code{fit} (the one-row fit-indices
#'   data.frame, restricted to the headline columns), and, for the MCMC
#'   engine only, \code{mcmc_info} (chain configuration and convergence
#'   diagnostics).
#' @seealso \code{\link{robust_lpa}}, \code{\link{print.robust_lpa}}
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' fit <- suppressWarnings(robust_lpa(x, G = 2, model = 2, n_starts = 3, max_iter = 20))
#' summary(fit)
#' @export
summary.robust_lpa <- function(object, ...) {
  G <- length(object$means)
  p <- length(object$means[[1]])
  var_names <- names(object$means[[1]])
  if (is.null(var_names)) var_names <- paste0("V", seq_len(p))
  profile_names <- paste0("P", seq_len(G))

  means_mat <- matrix(NA_real_, nrow = p, ncol = G, dimnames = list(var_names, profile_names))
  for (g in seq_len(G)) {
    means_mat[, g] <- as.numeric(object$means[[g]])
  }

  sizes_tab <- table(factor(object$assignments, levels = seq_len(G)))
  size_df <- data.frame(
    Profile = profile_names,
    N = as.integer(sizes_tab),
    Proportion = round(as.numeric(object$proportions), 3),
    row.names = NULL
  )

  mcmc_info <- NULL
  if (identical(object$engine, "MCMC")) {
    mcmc_info <- list(
      n_chains = if (!is.null(object$mcmc_draws)) object$mcmc_draws$n_chains else NA,
      mcmc_iter = if (!is.null(object$mcmc_draws)) object$mcmc_draws$mcmc_iter else NA,
      diagnostics = object$mcmc_diagnostics
    )
  }

  out <- list(
    engine = object$engine,
    model = object$fit$Model,
    G = G,
    n = nrow(object$probabilities),
    means = means_mat,
    sizes = size_df,
    fit = object$fit[, c("LogLik", "AIC", "BIC", "Entropy"), drop = FALSE],
    mcmc_info = mcmc_info
  )
  class(out) <- "summary.robust_lpa"
  out
}

#' Print a Summarized Robust LPA Model
#'
#' Prints the object returned by \code{\link{summary.robust_lpa}}: profile
#' means, profile sizes/mixing proportions, the headline fit indices, and, for
#' the MCMC engine, the chain configuration and Gelman-Rubin \eqn{\hat{R}} /
#' effective sample size ranges.
#'
#' @param x An object of class \code{"summary.robust_lpa"}, as returned by
#'   \code{\link{summary.robust_lpa}}.
#' @param digits Integer, number of decimal places to display. Default \code{2}.
#' @param ... Currently ignored (present for S3 consistency with the generic
#'   \code{\link[base]{print}}).
#' @return \code{x}, invisibly.
#' @keywords internal
#' @export
print.summary.robust_lpa <- function(x, digits = 2, ...) {
  cat(sprintf("robust_lpa summary -- %s | model %d | G = %d | N = %d\n", x$engine, x$model, x$G, x$n))

  cat("\nProfile means:\n")
  print(round(x$means, digits))

  cat("\nProfile sizes:\n")
  print(x$sizes, row.names = FALSE)

  cat(sprintf(
    "\nFit: LogLik = %.1f | AIC = %.1f | BIC = %.1f | Entropy = %.3f\n",
    x$fit$LogLik, x$fit$AIC, x$fit$BIC, x$fit$Entropy
  ))

  if (!is.null(x$mcmc_info)) {
    diag <- x$mcmc_info$diagnostics
    if (!is.null(diag)) {
      rhat_range <- suppressWarnings(range(diag$Rhat, na.rm = TRUE))
      ess_range <- suppressWarnings(range(diag$ESS, na.rm = TRUE))
      cat(sprintf(
        "MCMC: %s chains x %s iter | Rhat [%.2f, %.2f] | ESS [%.0f, %.0f]",
        x$mcmc_info$n_chains, x$mcmc_info$mcmc_iter,
        rhat_range[1], rhat_range[2], ess_range[1], ess_range[2]
      ))
      if (isTRUE(any(diag$Rhat > 1.1, na.rm = TRUE))) cat("  [!] some Rhat > 1.1")
      cat("\n")
    } else {
      cat("MCMC: convergence diagnostics not available (see ?robust_lpa).\n")
    }
  }

  invisible(x)
}
