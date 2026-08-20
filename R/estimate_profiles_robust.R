#' Estimate Robust Latent Profile Models Across Profiles and Models
#'
#' @param data A matrix or data.frame.
#' @param n_profiles A vector of integers specifying the number of profiles to run.
#' @param models A vector of LPA models to run.
#' @param engine String. Either "EM" or "MCMC".
#' @param cores Integer. Number of CPU cores to use for parallel processing
#'   \emph{across} the requested \code{n_profiles} x \code{models} grid (and
#'   across \code{k_folds} within \code{tune_lasso}). This is a different
#'   axis of parallelism from \code{\link{robust_lpa}}'s own \code{cores}
#'   argument (which parallelizes across EM restarts / MCMC chains
#'   \emph{within} a single model fit); if you also pass \code{cores} through
#'   \code{...} to \code{\link{robust_lpa}}, keep the product of the two
#'   roughly at or below your machine's core count to avoid oversubscription.
#' @param n_starts Number of initializations per model.
#' @param lambda Fixed penalty for LASSO.
#' @param tune_lasso Logical. If TRUE, finds optimal lambda via cross-validation.
#' @param k_folds Number of folds for cross-validation.
#' @param lambda_grid Vector of penalty values to test.
#' @param ... Additional arguments passed on to every internal
#'   \code{\link{robust_lpa}} call (e.g. \code{max_iter}, \code{tol},
#'   \code{robust}, \code{alpha}, \code{cores}, and, when
#'   \code{engine = "MCMC"}, \code{mcmc_iter}/\code{n_chains}/\code{prior_laplace})
#'   -- for instance, pass \code{robust = FALSE} here to compare
#'   model/profile combinations using classical (non-robust) estimation
#'   throughout, for either engine.
#' @return A list containing the fit comparison table and the estimated models.
#' @examples
#' # Quick evaluation of multiple profiles
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' res <- suppressWarnings(estimate_profiles_robust(x, n_profiles = 1:2, models = 1, n_starts = 3))
#' res$fit_table
#' # Each element of `res$models` is a `robust_lpa` object with print/summary methods
#' summary(res$models[[1]])
#' @export
estimate_profiles_robust <- function(data, n_profiles = 1:3, models = c(1, 2, 3, 4, 5, 6), 
                                     engine = "EM", cores = 1,
                                     n_starts = 5, lambda = 0, tune_lasso = FALSE, 
                                     k_folds = 5, lambda_grid = c(0.01, 0.05, 0.1, 0.2),
                                     ...) {
  
  X <- as.matrix(data)
  n <- nrow(X)
  has_na <- any(is.na(X))
  grid <- expand.grid(G = n_profiles, Model = models)
  
  worker_fun <- function(row_idx) {
    g <- grid$G[row_idx]
    m <- grid$Model[row_idx]
    current_lambda <- lambda
    
    if (tune_lasso && g > 1 && engine == "EM") {
      folds <- sample(rep(1:k_folds, length.out = n))
      best_cv_loglik <- -Inf
      best_l <- lambda_grid[1]
      
      for (l_test in lambda_grid) {
        fold_logliks <- numeric(k_folds)
        
        for (k in 1:k_folds) {
          train_data <- X[folds != k, , drop = FALSE]
          test_data <- X[folds == k, , drop = FALSE]
          
          fit_train <- tryCatch({
            robust_lpa(data = train_data, G = g, model = m, engine = engine, n_starts = max(1, n_starts - 2), lambda = l_test, ...)
          }, error = function(e) NULL)
          
          if (is.null(fit_train)) {
            fold_logliks[k] <- -Inf
            next
          }
          
          test_n <- nrow(test_data)
          density_matrix <- matrix(0, nrow = test_n, ncol = g)
          
          for (c_g in 1:g) {
            if (!has_na) {
              density_matrix[, c_g] <- dmvnorm_cpp(test_data, fit_train$means[[c_g]], fit_train$covariances[[c_g]]) * fit_train$proportions[c_g]
            } else {
              density_matrix[, c_g] <- dmvnorm_fiml_cpp(test_data, fit_train$means[[c_g]], fit_train$covariances[[c_g]]) * fit_train$proportions[c_g]
            }
          }
          
          row_sums <- rowSums(density_matrix)
          row_sums[row_sums < 1e-300] <- 1e-300
          fold_logliks[k] <- sum(log(row_sums))
        }
        
        mean_ll <- mean(fold_logliks[fold_logliks > -Inf])
        if (!is.na(mean_ll) && mean_ll > best_cv_loglik) {
          best_cv_loglik <- mean_ll
          best_l <- l_test
        }
      }
      current_lambda <- best_l
    }
    
    fit_out <- tryCatch({
      robust_lpa(data = X, G = g, model = m, engine = engine, n_starts = n_starts, lambda = current_lambda, ...)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(fit_out)) {
      fit_out$fit$Lambda <- current_lambda
      if (cores == 1) {
        message(sprintf("Model %d with %d profile(s) estimated.", m, g))
      }
      return(list(name = paste0("model_", m, "_profiles_", g), 
                  model = m, profiles = g, success = TRUE,
                  fit_row = fit_out$fit, model_data = fit_out))
    }
    
    if (cores == 1) {
      message(sprintf("Model %d with %d profile(s) failed.", m, g))
    }
    return(list(model = m, profiles = g, success = FALSE))
  }
  
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("RobustLPA", quietly = TRUE)) {
          library(RobustLPA)
        }
      })
      
      parallel::clusterExport(
        cl,
        varlist = c("X", "grid", "engine", "n_starts", "lambda", "has_na", 
                    "tune_lasso", "k_folds", "lambda_grid", "robust_lpa", "robust_m_step"),
        envir = environment()
      )
      
      results <- parallel::parLapply(cl, 1:nrow(grid), worker_fun)
      parallel::stopCluster(cl)
    } else {
      results <- parallel::mclapply(1:nrow(grid), worker_fun, mc.cores = cores)
    }
  } else {
    results <- lapply(1:nrow(grid), worker_fun)
  }
  
  master_fit_table <- data.frame()
  model_list <- list()
  
  for (res in results) {
    if (!is.null(res)) {
      if (cores > 1) {
        if (isTRUE(res$success)) {
          message(sprintf("Model %d with %d profile(s) estimated.", res$model, res$profiles))
        } else {
          message(sprintf("Model %d with %d profile(s) failed.", res$model, res$profiles))
        }
      }
      
      if (isTRUE(res$success)) {
        master_fit_table <- rbind(master_fit_table, res$fit_row)
        model_list[[res$name]] <- res$model_data
      }
    }
  }
  
  return(list(fit_table = master_fit_table, models = model_list))
}