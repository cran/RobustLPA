#' Estimate Robust Latent Profile Models Across Profiles and Models
#'
#' Fits every combination of \code{n_profiles} and \code{models} with
#' \code{\link{robust_lpa}} and collects their fit indices in one table.
#'
#' @param data A matrix or data.frame.
#' @param n_profiles A vector of integers specifying the number of profiles to run.
#' @param models A vector of LPA models to run.
#' @param engine String. Either "EM" or "MCMC".
#' @param cores Integer. Number of CPU cores to use for parallel processing
#'   \emph{across} the requested \code{n_profiles} x \code{models} grid. This
#'   is a different axis of parallelism from \code{\link{robust_lpa}}'s own
#'   \code{cores} argument; if you also pass \code{cores} through \code{...},
#'   keep the product of the two at or below your machine's core count.
#'   Results are identical for any value of \code{cores} given the same
#'   \code{set.seed()}.
#' @param n_starts Number of initializations per model.
#' @param lambda Fixed penalty for LASSO (EM engine).
#' @param tune_lasso Logical. If TRUE, selects \code{lambda} from
#'   \code{lambda_grid} by k-fold cross-validation of the held-out
#'   observed-data log-likelihood (EM engine only; ignored with a warning
#'   for \code{engine = "MCMC"}, whose shrinkage is governed by
#'   \code{prior_laplace}).
#' @param k_folds Number of folds for cross-validation.
#' @param lambda_grid Vector of penalty values to test.
#' @param ... Additional arguments passed on to every internal
#'   \code{\link{robust_lpa}} call (e.g. \code{max_iter}, \code{tol},
#'   \code{robust}, \code{robust_method}, \code{nu}, \code{alpha},
#'   \code{init}, and, when \code{engine = "MCMC"},
#'   \code{mcmc_iter}/\code{n_chains}/\code{prior_laplace}).
#' @return A list containing \code{fit_table} (one row per successfully
#'   fitted model, with a \code{Lambda} column) and \code{models} (the
#'   fitted \code{robust_lpa} objects, named \code{"model_<m>_profiles_<G>"}).
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' set.seed(1)
#' res <- estimate_profiles_robust(x, n_profiles = 1:2, models = 1, n_starts = 3)
#' res$fit_table
#' summary(res$models[[1]])
#' @export
estimate_profiles_robust <- function(data, n_profiles = 1:3, models = c(1, 2, 3, 4, 5, 6),
                                     engine = "EM", cores = 1,
                                     n_starts = 5, lambda = 0, tune_lasso = FALSE,
                                     k_folds = 5, lambda_grid = c(0.01, 0.05, 0.1, 0.2),
                                     ...) {
  if (!(length(engine) == 1 && engine %in% c("EM", "MCMC"))) {
    stop("`engine` must be either 'EM' or 'MCMC'.")
  }
  if (!is.numeric(n_profiles) || any(n_profiles < 1) || any(n_profiles != round(n_profiles))) {
    stop("`n_profiles` must contain positive integers.")
  }
  if (!all(models %in% 1:6)) stop("`models` must contain integers between 1 and 6.")
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive integer.")
  }
  if (tune_lasso) {
    if (engine == "MCMC") {
      warning("`tune_lasso = TRUE` applies to the EM engine only (the MCMC engine uses `prior_laplace`); ignoring it.")
      tune_lasso <- FALSE
    } else if (!is.numeric(k_folds) || length(k_folds) != 1 || k_folds < 2) {
      stop("`k_folds` must be a single integer of at least 2.")
    } else if (!is.numeric(lambda_grid) || length(lambda_grid) < 1 || any(lambda_grid < 0)) {
      stop("`lambda_grid` must be a vector of non-negative numbers.")
    }
  }

  X <- as.matrix(data)
  storage.mode(X) <- "double"
  n <- nrow(X)
  grid <- expand.grid(G = n_profiles, Model = models)

  worker_fun <- function(row_idx) {
    g <- grid$G[row_idx]
    m <- grid$Model[row_idx]
    current_lambda <- lambda

    if (tune_lasso && g > 1) {
      folds <- sample(rep(seq_len(k_folds), length.out = n))
      best_cv_loglik <- -Inf
      best_l <- lambda_grid[1]

      for (l_test in lambda_grid) {
        fold_logliks <- rep(NA_real_, k_folds)
        for (k in seq_len(k_folds)) {
          train_data <- X[folds != k, , drop = FALSE]
          test_data <- X[folds == k, , drop = FALSE]
          fit_train <- tryCatch(
            suppressWarnings(robust_lpa(data = train_data, G = g, model = m, engine = engine,
                                        n_starts = max(1, n_starts - 2), lambda = l_test, ...)),
            error = function(e) NULL
          )
          if (!is.null(fit_train)) fold_logliks[k] <- .loglik_newdata(fit_train, test_data)
        }
        mean_ll <- mean(fold_logliks, na.rm = TRUE)
        if (is.finite(mean_ll) && mean_ll > best_cv_loglik) {
          best_cv_loglik <- mean_ll
          best_l <- l_test
        }
      }
      current_lambda <- best_l
    }

    fit_out <- tryCatch(
      robust_lpa(data = X, G = g, model = m, engine = engine, n_starts = n_starts,
                 lambda = current_lambda, ...),
      error = function(e) e
    )
    if (inherits(fit_out, "error")) {
      return(list(model = m, profiles = g, success = FALSE, error = conditionMessage(fit_out)))
    }
    fit_out$fit$Lambda <- current_lambda
    list(name = paste0("model_", m, "_profiles_", g), model = m, profiles = g, success = TRUE,
         fit_row = fit_out$fit, model_data = fit_out)
  }

  results <- .run_parallel(
    cores, nrow(grid), worker_fun,
    export_vars = c("X", "n", "grid", "engine", "n_starts", "lambda", "tune_lasso", "k_folds", "lambda_grid"),
    export_env = environment()
  )

  master_fit_table <- data.frame()
  model_list <- list()
  for (res in results) {
    if (isTRUE(res$success)) {
      message(sprintf("Model %d with %d profile(s) estimated.", res$model, res$profiles))
      master_fit_table <- rbind(master_fit_table, res$fit_row)
      model_list[[res$name]] <- res$model_data
    } else {
      message(sprintf("Model %d with %d profile(s) failed: %s", res$model, res$profiles, res$error))
    }
  }
  rownames(master_fit_table) <- NULL
  list(fit_table = master_fit_table, models = model_list)
}
