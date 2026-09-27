#' Print a Fitted Robust Growth Mixture Model
#'
#' Compact overview of a model fitted by \code{\link{robust_gmm}}: engine,
#' number of classes and persons, outcomes, random-effect and robustness
#' settings, headline fit indices, class proportions and, if used, the
#' LASSO penalties.
#'
#' @param x A \code{robust_gmm} object.
#' @param ... Currently ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{robust_gmm}}, \code{\link{summary.robust_gmm}}
#' @export
print.robust_gmm <- function(x, ...) {
  sp <- x$spec
  cat(sprintf("<robust_gmm> %s | G = %d | persons = %d | outcomes: %s\n",
              x$engine, sp$G, length(x$ids), paste(sp$outcomes, collapse = ", ")))
  cat(sprintf("Trajectory: degree %d | random: %s%s | %s\n", sp$degree, sp$random,
              if (sp$nrand > 0) sprintf(" (%s, %s across classes)", sp$re_structure, sp$re_cov) else "",
              .describe_method(x)))
  cat(sprintf("LogLik = %.1f | BIC = %.1f | Entropy = %.3f%s\n", x$fit$LogLik, x$fit$BIC, x$fit$Entropy,
              if (!is.null(x$fit$WAIC)) sprintf(" | WAIC = %.1f", x$fit$WAIC) else ""))
  props <- round(x$proportions, 3)
  cat("Proportions: ", paste(sprintf("C%d=%.2f", seq_along(props), props), collapse = ", "), "\n", sep = "")
  if (x$penalty$lambda_growth > 0 || x$penalty$lambda_diff > 0) {
    cat(sprintf("LASSO: growth = %g, difference = %g%s%s%s\n", x$penalty$lambda_growth, x$penalty$lambda_diff,
                if (isTRUE(x$penalty$group_diff)) " (group)" else "",
                if (isTRUE(x$penalty$adaptive)) ", adaptive" else "",
                if (isTRUE(x$penalty$relaxed)) ", relaxed refit" else ""))
  }
  if (isFALSE(x$converged)) cat(sprintf("EM did not converge within %d iterations.\n", x$iterations))
  invisible(x)
}

#' Summarize a Fitted Robust Growth Mixture Model
#'
#' @param object A \code{robust_gmm} object.
#' @param ... Currently ignored.
#' @return An object of class \code{"summary.robust_gmm"}: a list with the
#'   class trajectories (a data.frame with one row per class and outcome),
#'   class sizes, variance components, fit indices, penalty information and,
#'   for the MCMC engine, the range of the convergence diagnostics.
#' @seealso \code{\link{robust_gmm}}
#' @examples
#' data(neuro_long)
#' set.seed(1)
#' fit <- robust_gmm(neuro_long, id = "ID", time = "Year", outcomes = "Memory",
#'                   G = 2, n_starts = 2)
#' summary(fit)
#' @export
summary.robust_gmm <- function(object, ...) {
  sp <- object$spec
  G <- sp$G
  traj <- do.call(rbind, lapply(seq_len(G), function(g) {
    m <- object$coefficients[[g]]
    data.frame(Class = g, Outcome = rownames(m), m, check.names = FALSE, row.names = NULL)
  }))
  sizes <- data.frame(Class = seq_len(G), N = tabulate(object$assignments, nbins = G),
                      Proportion = round(object$proportions, 3))
  out <- list(
    engine = object$engine, method = .describe_method(object), spec = sp,
    n = length(object$ids), trajectories = traj, sizes = sizes,
    random_cov = object$random_cov, residual_var = object$residual_var,
    fit = object$fit, penalty = object$penalty, converged = object$converged,
    iterations = object$iterations, diagnostics = object$mcmc_diagnostics
  )
  class(out) <- "summary.robust_gmm"
  out
}

#' Print a Summarized Robust Growth Mixture Model
#'
#' @param x An object of class \code{"summary.robust_gmm"}.
#' @param digits Number of decimal places. Default \code{3}.
#' @param ... Currently ignored.
#' @return \code{x}, invisibly.
#' @keywords internal
#' @export
print.summary.robust_gmm <- function(x, digits = 3, ...) {
  sp <- x$spec
  cat(sprintf("robust_gmm summary -- %s | G = %d | persons = %d\n", x$engine, sp$G, x$n))
  cat(sprintf("Estimation: %s | degree %d | random: %s\n", x$method, sp$degree, sp$random))
  cat("\nClass trajectories (original scale):\n")
  tr <- x$trajectories
  num <- vapply(tr, is.numeric, logical(1)) & names(tr) != "Class"
  tr[num] <- lapply(tr[num], round, digits)
  print(tr, row.names = FALSE)
  cat("\nClass sizes (modal assignment):\n")
  print(x$sizes, row.names = FALSE)
  if (!is.null(x$random_cov)) {
    cat(sprintf("\nRandom-effect covariance (%s across classes; class 1 shown):\n", sp$re_cov))
    print(round(x$random_cov[[1]], digits))
  }
  cat("\nResidual variances:\n")
  print(round(x$residual_var, digits))
  f <- x$fit
  cat(sprintf("\nFit: LogLik = %.2f | Parameters = %d | AIC = %.1f | BIC = %.1f | SABIC = %.1f | Entropy = %.3f%s\n",
              f$LogLik, as.integer(f$Parameters), f$AIC, f$BIC, f$SABIC, f$Entropy,
              if (!is.null(f$WAIC)) sprintf(" | WAIC = %.1f", f$WAIC) else ""))
  pen <- x$penalty
  if (pen$lambda_growth > 0 || pen$lambda_diff > 0) {
    cat(sprintf("\nLASSO: growth = %g, difference = %g%s%s%s\n", pen$lambda_growth, pen$lambda_diff,
                if (isTRUE(pen$group_diff)) " (group by outcome)" else "",
                if (isTRUE(pen$adaptive)) ", adaptive weights" else "",
                if (isTRUE(pen$relaxed)) ", estimates from the relaxed (unpenalized) refit" else ""))
    if (!is.null(pen$outcome_selected)) {
      cat("Outcomes differentiating the classes: ",
          paste(names(pen$outcome_selected)[pen$outcome_selected], collapse = ", "),
          if (any(!pen$outcome_selected)) paste0(" (removed: ", paste(names(pen$outcome_selected)[!pen$outcome_selected], collapse = ", "), ")") else "",
          "\n", sep = "")
    }
    if (!is.null(pen$pattern) && any(pen$pattern$growth_zero)) {
      cat(sprintf("Growth terms set to zero: %d of %d\n", sum(pen$pattern$growth_zero), length(pen$pattern$growth_zero)))
    }
  }
  if (isFALSE(x$converged)) cat(sprintf("Note: EM did not converge within %d iterations.\n", x$iterations))
  if (!is.null(x$diagnostics)) {
    d <- x$diagnostics
    cat(sprintf("MCMC: Rhat [%.2f, %.2f] | ESS [%.0f, %.0f]%s\n",
                min(d$Rhat, na.rm = TRUE), max(d$Rhat, na.rm = TRUE), min(d$ESS, na.rm = TRUE),
                max(d$ESS, na.rm = TRUE), if (isTRUE(any(d$Rhat > 1.1, na.rm = TRUE))) "  [!] some Rhat > 1.1" else ""))
  }
  invisible(x)
}

#' Plot the Class Trajectories of a Robust Growth Mixture Model
#'
#' Draws, for every outcome, the estimated mean trajectory of each latent
#' class over the observed time range, optionally over the observed
#' individual trajectories coloured by modal class.
#'
#' @param model A \code{robust_gmm} object.
#' @param outcomes Optional subset of outcomes to plot (default: all).
#' @param individual Logical; draw the observed individual trajectories
#'   (default \code{TRUE}).
#' @param max_individuals Maximum number of persons whose trajectories are
#'   drawn (a random subset if there are more). Default \code{300}.
#' @param title Plot title.
#' @return A ggplot object.
#' @seealso \code{\link{robust_gmm}}
#' @examples
#' data(neuro_long)
#' set.seed(1)
#' fit <- robust_gmm(neuro_long, id = "ID", time = "Year", outcomes = "Memory",
#'                   G = 2, n_starts = 2)
#' plot_robust_gmm(fit)
#' @export
plot_robust_gmm <- function(model, outcomes = NULL, individual = TRUE, max_individuals = 300,
                            title = "Latent class trajectories") {
  if (!inherits(model, "robust_gmm")) stop("`model` must be a robust_gmm() fit.")
  sp <- model$spec
  outs <- outcomes %||% sp$outcomes
  if (!all(outs %in% sp$outcomes)) stop("Unknown outcome(s): ", paste(setdiff(outs, sp$outcomes), collapse = ", "))
  dat <- model$data
  tt <- dat[[sp$time]]
  grid <- seq(min(tt, na.rm = TRUE), max(tt, na.rm = TRUE), length.out = 60)
  curves <- do.call(rbind, lapply(seq_len(sp$G), function(g) {
    do.call(rbind, lapply(outs, function(o) {
      b <- model$coefficients[[g]][o, ]
      data.frame(Time = grid, Value = as.vector(outer(grid, 0:sp$degree, `^`) %*% b),
                 Outcome = factor(o, levels = outs), Class = factor(g, levels = seq_len(sp$G)))
    }))
  }))
  p <- ggplot2::ggplot()
  if (individual) {
    cls <- model$assignments[match(dat[[sp$id]], model$ids)]
    keep_ids <- model$ids
    if (length(keep_ids) > max_individuals) keep_ids <- sample(keep_ids, max_individuals)
    ind <- do.call(rbind, lapply(outs, function(o) {
      sel <- dat[[sp$id]] %in% keep_ids & !is.na(dat[[o]]) & !is.na(cls)
      data.frame(Time = dat[[sp$time]][sel], Value = dat[[o]][sel], Outcome = factor(o, levels = outs),
                 Person = factor(dat[[sp$id]][sel]), Class = factor(cls[sel], levels = seq_len(sp$G)))
    }))
    ind <- ind[order(ind$Person, ind$Time), ]
    p <- p + ggplot2::geom_line(data = ind, ggplot2::aes(x = Time, y = Value, group = Person, colour = Class),
                                alpha = 0.15, linewidth = 0.3)
  }
  p +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = Time, y = Value, colour = Class), linewidth = 1.3) +
    ggplot2::facet_wrap(~Outcome, scales = "free_y") +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = sp$time, y = "Value", colour = "Class") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
}
