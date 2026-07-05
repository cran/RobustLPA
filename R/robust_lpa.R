#' Fit a Single Robust Latent Profile Analysis Model
#'
#' This function estimates a single robust Latent Profile Analysis (LPA) model
#' for a specified number of profiles and model structure. It automatically
#' handles missing data via robust FIML if present.
#'
#' @param data A matrix or data.frame of observations.
#' @param G The number of latent profiles to extract.
#' @param model An integer (1 to 6) specifying the model parameterization.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Tolerance for convergence.
#' @param n_starts Number of random initializations.
#' @return A list containing parameters, fit indices, and assignments.
#' @examples
#' # Quick example using a small subset of the iris dataset
#' data(iris)
#' fit <- robust_lpa(data = iris[1:30, 1:2], G = 2, model = 1, n_starts = 1)
#' # Print the fit indices
#' fit$fit
#' @export
robust_lpa <- function(data, G, model = 6, max_iter = 100, tol = 1e-6, n_starts = 5) {
  X <- as.matrix(data)
  n <- nrow(X)
  p <- ncol(X)

  best_log_lik <- -Inf
  best_model <- NULL

  # Helper function to mathematically guarantee Positive-Definiteness
  force_pd <- function(mat) {
    mat <- (mat + t(mat)) / 2 # Ensure perfect symmetry
    ev <- eigen(mat, symmetric = TRUE)$values
    min_ev <- min(ev)
    if (min_ev < 1e-5) {
      mat <- mat + diag(abs(min_ev) + 1e-4, nrow(mat))
    }
    return(mat)
  }

  # Preventive check for NA presence
  has_na <- any(is.na(X))

  for(start in 1:n_starts) {
    # INITIALIZATION
    pi_g <- rep(1 / G, G)
    mu <- list()
    sigma <- list()

    z_init <- matrix(runif(n * G), nrow = n, ncol = G)
    z_init <- z_init / rowSums(z_init)

    for(g in 1:G) {
      init_results <- robust_m_step(data = X, z = z_init[, g], alpha = 0.05)
      mu[[g]] <- init_results$mean
      sigma[[g]] <- force_pd(init_results$covariance)
    }

    log_lik <- -Inf

    # EM LOOP
    for(iter in 1:max_iter) {

      # --- E-STEP (WITH FULL FIML C++ OPTIMIZATION) ---
      density_matrix <- matrix(0, nrow = n, ncol = G)

      if (!has_na) {
        # Fast path without Missing Data
        for(g in 1:G) {
          density_matrix[, g] <- dmvnorm_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
        }
      } else {
        # C++ FIML path (Extremely Fast)
        for(g in 1:G) {
          density_matrix[, g] <- dmvnorm_fiml_cpp(X, mu[[g]], sigma[[g]]) * pi_g[g]
        }
      }

      row_sums <- rowSums(density_matrix)
      row_sums[row_sums < 1e-300] <- 1e-300

      new_log_lik <- sum(log(row_sums))
      z <- density_matrix / row_sums

      if(abs(new_log_lik - log_lik) < tol) break
      log_lik <- new_log_lik

      # --- M-STEP ---
      raw_sigmas <- list()
      for(g in 1:G) {
        m_step_results <- robust_m_step(data = X, z = z[, g], alpha = 0.05)
        mu[[g]] <- m_step_results$mean
        raw_sigmas[[g]] <- m_step_results$covariance
        pi_g[g] <- mean(z[, g])
      }

      # Apply tidyLPA geometric constraints (Models 1 to 6)
      if(model == 1) {
        pooled_diag <- numeric(p)
        for(g in 1:G) { pooled_diag <- pooled_diag + pi_g[g] * diag(raw_sigmas[[g]]) }
        for(g in 1:G) { sigma[[g]] <- diag(pooled_diag, p) }

      } else if(model == 2) {
        for(g in 1:G) { sigma[[g]] <- diag(diag(raw_sigmas[[g]]), p) }

      } else if(model == 3) {
        pooled_sigma <- matrix(0, nrow = p, ncol = p)
        for(g in 1:G) { pooled_sigma <- pooled_sigma + pi_g[g] * raw_sigmas[[g]] }
        pooled_sigma <- force_pd(pooled_sigma) # Ensure PD
        for(g in 1:G) { sigma[[g]] <- pooled_sigma }

      } else if(model == 4) {
        pooled_cov <- matrix(0, nrow = p, ncol = p)
        for(g in 1:G) { pooled_cov <- pooled_cov + pi_g[g] * raw_sigmas[[g]] }
        for(g in 1:G) {
          sigma[[g]] <- pooled_cov
          diag(sigma[[g]]) <- diag(raw_sigmas[[g]])
          sigma[[g]] <- force_pd(sigma[[g]]) # Ensure Hybrid Matrix is PD
        }

      } else if(model == 5) {
        pooled_diag <- numeric(p)
        for(g in 1:G) { pooled_diag <- pooled_diag + pi_g[g] * diag(raw_sigmas[[g]]) }
        for(g in 1:G) {
          sigma[[g]] <- raw_sigmas[[g]]
          diag(sigma[[g]]) <- pooled_diag
          sigma[[g]] <- force_pd(sigma[[g]]) # Ensure Hybrid Matrix is PD
        }

      } else if(model == 6) {
        for(g in 1:G) { sigma[[g]] <- force_pd(raw_sigmas[[g]]) }
      }
    }

    # SAVE BEST START
    if(log_lik > best_log_lik) {
      best_log_lik <- log_lik

      # Calculate K (Number of estimated parameters)
      if(G == 1) {
        K <- switch(as.character(model),
                    "1" = p + p,
                    "2" = p + p,
                    "3" = p + (p * (p + 1)/2),
                    "4" = p + (p * (p + 1)/2),
                    "5" = p + (p * (p + 1)/2),
                    "6" = p + (p * (p + 1)/2))
      } else {
        K <- switch(as.character(model),
                    "1" = (G * p) + p + (G - 1),
                    "2" = (G * p) + (G * p) + (G - 1),
                    "3" = (G * p) + (p * (p + 1) / 2) + (G - 1),
                    "4" = (G * p) + (G * p) + (p * (p - 1) / 2) + (G - 1),
                    "5" = (G * p) + p + (G * p * (p - 1) / 2) + (G - 1),
                    "6" = (G * p) + (G * (p * (p + 1) / 2)) + (G - 1))
      }

      # Fit Indices Calculation
      AIC_val <- -2 * best_log_lik + 2 * K
      BIC_val <- -2 * best_log_lik + K * log(n)
      SABIC_val <- -2 * best_log_lik + K * log((n + 2) / 24)

      z_safe <- z; z_safe[z_safe < 1e-15] <- 1e-15
      entropy <- 1 - (sum(-z * log(z_safe)) / (n * log(G)))
      if(G == 1) entropy <- 1

      assignments <- max.col(z)
      sizes <- table(factor(assignments, levels = 1:G))

      fit_indices <- data.frame(
        Model = model,
        Profiles = G,
        LogLik = best_log_lik,
        Parameters = K,
        AIC = AIC_val,
        BIC = BIC_val,
        SABIC = SABIC_val,
        Entropy = entropy,
        Min_Size = min(sizes) / n,
        Max_Size = max(sizes) / n
      )

      best_model <- list(means = mu, covariances = sigma, proportions = pi_g,
                         probabilities = z, fit = fit_indices, assignments = assignments)
    }
  }
  return(best_model)
}
