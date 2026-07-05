#' Estimate Robust Latent Profile Models
#'
#' This function provides a user-friendly interface to fit multiple robust LPA
#' models simultaneously across different numbers of profiles and model structures,
#' returning a comprehensive fit summary table.
#'
#' @param data A matrix or data.frame.
#' @param n_profiles A vector of integers specifying the number of profiles to run (e.g., 1:3).
#' @param models A vector of LPA models to run (e.g., c(1, 2, 3, 4, 5, 6)). Default is c(1, 2, 3, 4, 5, 6).
#' @param n_starts Number of initializations per model.
#' @return A list containing the fit comparison table and the estimated models.
#' @examples
#' # Quick evaluation of multiple profiles
#' data(iris)
#' res <- estimate_profiles_robust(iris[1:30, 1:2], n_profiles = 1:2, models = 1, n_starts = 1)
#' res$fit_table
#' @export
estimate_profiles_robust <- function(data, n_profiles = 1:3, models = c(1, 2, 3, 4, 5, 6), n_starts = 5) {
  master_fit_table <- data.frame()
  model_list <- list()

  message("Running Robust LPA Models...")

  for(m in models) {
    for(g in n_profiles) {
      message(paste0("Fitting Model ", m, " with ", g, " Profiles... "))

      fit_out <- tryCatch({
        robust_lpa(data = data, G = g, model = m, n_starts = n_starts)
      }, error = function(e) {
        message(paste("[FATAL ERROR]:", e$message, "- "))
        NULL
      })

      if(!is.null(fit_out)) {
        master_fit_table <- rbind(master_fit_table, fit_out$fit)
        model_name <- paste0("model_", m, "_profiles_", g)
        model_list[[model_name]] <- fit_out
        message("Done.")
      } else {
        message("Failed.")
      }
    }
  }

  message("All models finished.")
  return(list(fit_table = master_fit_table, models = model_list))
}
