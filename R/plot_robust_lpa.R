#' Plot Robust Latent Profiles
#'
#' Automatically generates a professional profile plot using ggplot2 from an estimated
#' robust LPA model.
#'
#' @param model A model object returned by robust_lpa or estimate_profiles_robust.
#' @param title The title of the plot. Default is "Robust Latent Profiles".
#' @param xlab The x-axis label. Default is "Variables".
#' @param ylab The y-axis label. Default is "Value".
#' @param var_labels A character vector to manually rename the variables on the X axis. Default is NULL (auto-detect).
#' @param legend_title The title of the legend. Default is "Class".
#' @return A ggplot object.
#' @examples
#' data(iris)
#' fit <- robust_lpa(data = iris[1:30, 1:2], G = 2, model = 1, n_starts = 1)
#' plot_robust_lpa(fit)
#' @export
plot_robust_lpa <- function(model,
                            title = "Robust Latent Profiles",
                            xlab = "Variables",
                            ylab = "Value",
                            var_labels = NULL,
                            legend_title = "Class") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2' to use this function.")
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
    ggplot2::scale_fill_brewer(palette = "Set1")

  return(p_plot)
}
