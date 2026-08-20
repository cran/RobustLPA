#' Plot Robust Latent Profiles
#'
#' Automatically generates a professional profile plot using ggplot2 from an
#' estimated robust LPA model.
#'
#' @param model Either a single fitted model object returned by
#'   \code{\link{robust_lpa}}, or the list returned by
#'   \code{\link{estimate_profiles_robust}}. In the latter case, the model
#'   with the lowest BIC in \code{model$fit_table} is selected automatically
#'   (a message reports which one), unless \code{which_model} is given.
#' @param which_model Optional string, the name of a specific model to plot
#'   when \code{model} is an \code{\link{estimate_profiles_robust}} result
#'   (one of \code{names(model$models)}, e.g. \code{"model_6_profiles_3"}).
#'   Ignored when \code{model} is already a single fitted model.
#' @param title The title of the plot. Default is "Robust Latent Profiles".
#' @param xlab The x-axis label. Default is "Variables".
#' @param ylab The y-axis label. Default is "Value".
#' @param var_labels A character vector to manually rename the variables on the X axis. Default is NULL (auto-detect).
#' @param legend_title The title of the legend. Default is "Class".
#' @return A ggplot object.
#' @examples
#' data(neuro_data)
#' x <- scale(as.matrix(neuro_data[, c("Memory", "RT_Stroop")]))
#' fit <- suppressWarnings(robust_lpa(data = x, G = 2, model = 1, n_starts = 3))
#' print(fit)  # concise overview (print.robust_lpa())
#' plot_robust_lpa(fit)
#' @export
plot_robust_lpa <- function(model,
                            which_model = NULL,
                            title = "Robust Latent Profiles",
                            xlab = "Variables",
                            ylab = "Value",
                            var_labels = NULL,
                            legend_title = "Class") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2' to use this function.")
  }
  
  # ---- Accept either a single robust_lpa() fit, or an estimate_profiles_robust()
  # result (a list with $fit_table and $models); auto-select the best-BIC
  # model in the latter case. Without this, the documented claim that this
  # function accepts estimate_profiles_robust() output directly would not
  # actually work, since that object has no top-level $means.
  if (is.null(model$means)) {
    if (!is.null(model$models) && !is.null(model$fit_table)) {
      if (length(model$models) == 0) {
        stop("`model` is an estimate_profiles_robust() result with no successfully fitted models to plot.")
      }
      if (!is.null(which_model)) {
        if (!which_model %in% names(model$models)) {
          stop("`which_model` = '", which_model, "' not found. Available: ", paste(names(model$models), collapse = ", "))
        }
        selected_name <- which_model
      } else {
        best_idx <- which.min(model$fit_table$BIC)
        selected_name <- paste0("model_", model$fit_table$Model[best_idx], "_profiles_", model$fit_table$Profiles[best_idx])
        if (!selected_name %in% names(model$models)) {
          # Fall back to matching by row position if name construction ever
          # drifts from estimate_profiles_robust()'s own naming convention.
          selected_name <- names(model$models)[best_idx]
        }
        message("`model` looks like an estimate_profiles_robust() result; plotting the lowest-BIC fit ('", selected_name, "'). Pass `which_model` to choose a different one.")
      }
      model <- model$models[[selected_name]]
    } else {
      stop("`model` must be a fitted model from robust_lpa() (with a `$means` component), or the list returned by estimate_profiles_robust().")
    }
  }
  
  if (is.null(model$means) || is.null(model$covariances)) {
    stop("`model` does not look like a valid robust_lpa() fit (missing `$means` or `$covariances`).")
  }
  
  G <- length(model$means)
  p <- length(model$means[[1]])
  
  auto_names <- names(model$means[[1]])
  if (is.null(auto_names)) {
    auto_names <- paste0("Var", 1:p)
  }
  
  if (!is.null(var_labels)) {
    if (length(var_labels) == p) {
      axis_names <- var_labels
    } else {
      warning("Length of 'var_labels' does not match number of variables. Using automatic names.")
      axis_names <- auto_names
    }
  } else {
    axis_names <- auto_names
  }
  
  df_list <- list()
  for (g in 1:G) {
    means_g <- as.numeric(model$means[[g]])
    sds_g <- sqrt(diag(model$covariances[[g]]))
    
    df_list[[g]] <- data.frame(
      Class = as.factor(g),
      Variable = auto_names,
      Mean = means_g,
      SD = sds_g
    )
  }
  plot_data <- do.call(rbind, df_list)
  
  plot_data$Variable <- factor(plot_data$Variable, levels = auto_names, labels = axis_names)
  
  plot_data$Ymin <- plot_data$Mean - plot_data$SD
  plot_data$Ymax <- plot_data$Mean + plot_data$SD
  
  # Points use pch 21-25 (recycled if G > 5) so that the "fill = white, colour
  # = Class" styling below actually has a visible effect: the default shape
  # palette ggplot2 assigns to a mapped `shape` aesthetic (pch 15-18-ish,
  # solid glyphs) has no separate fill channel, which silently made the
  # `fill = "white"` in geom_point() a no-op in the previous version of this
  # function.
  shape_values <- rep(21:25, length.out = G)
  
  p_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Variable, y = Mean, group = Class, color = Class, shape = Class, fill = Class)) +
    
    ggplot2::geom_crossbar(ggplot2::aes(ymin = Ymin, ymax = Ymax),
                           width = 0.25,
                           alpha = 0.12,
                           linewidth = 0.2,
                           fatten = 0) +
    
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Ymin, ymax = Ymax),
                           width = 0.10,
                           linewidth = 0.4,
                           alpha = 0.6) +
    
    ggplot2::geom_line(linewidth = 0.7) +
    
    ggplot2::geom_point(size = 2.5, fill = "white", stroke = 1.0) +
    
    ggplot2::scale_shape_manual(values = shape_values) +
    
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title,
                  x = xlab,
                  y = ylab,
                  color = legend_title,
                  fill = legend_title,
                  shape = legend_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
      axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(size = 11, margin = ggplot2::margin(t = 10)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    
    ggplot2::scale_color_brewer(palette = "Set1") +
    # geom_point()'s literal fill = "white" above overrides this per-layer
    # (rendering white-centered pch 21-25 points), but geom_crossbar() has
    # no such override and needs this scale so its box fill matches the
    # same "Set1" palette used for colour/lines elsewhere.
    ggplot2::scale_fill_brewer(palette = "Set1")
  
  return(p_plot)
}