
#' Plot Hawkes Process Points
#'
#' Creates a scatterplot of points from a spatio-temporal Hawkes process, colored either by time or by whether they are background events.
#'
#' @param hawkes A `hawkes` object
#' @param color A character string. Either `"time"` to color by time, or `"background"` to color by whether an event is a background event (only works for simulated processes).
#' @param ... Additional arguments passed to `ggplot2::geom_point()`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' plot_hawkes(hawkes, color = "time")
plot_hawkes <- function(hawkes, color = "time", ...) {
  .unpack_hawkes(hawkes)

  if (color == "time") {
    plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = spatial_region) +
      ggplot2::geom_sf(data = hawkes, ggplot2::aes(color = t)) +
      labs(color = "Time")
  }
  if (color == "background") {
    plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = spatial_region) +
      ggplot2::geom_sf(data = hawkes |>
                         dplyr::mutate(background = factor(gen == 0, levels = c("TRUE", "FALSE"))),
                       aes(color = background)) +
      ggplot2::labs(color = "Background Event")
  }
  plot
}


#' Plot Conditional Intensity of Hawkes Process
#'
#' Visualizes the conditional intensity function at a fixed time over a spatial grid, based on a fitted Hawkes model.
#'
#' @param time A numeric value giving the time at which to evaluate the conditional intensity.
#' @param hawkes A `hawkes` object
#' @param est A list of fitted parameter values for the Hawkes process.
#' @param step A numeric value for the spatial resolution of the evaluation grid.
#'
#' @return A `ggplot2` object showing the spatial intensity at the given time, with observed events overlaid.
#' @export
#'
#' @examples
#' # example code
#' set.seed(123)
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.25),temporal = list(rate = 1), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 5))
#' plot_intensity(hawkes, est, time = 40)
#' plot_intensity(hawkes, est, coordinates = c(5,5))
#' plot_intensity(hawkes, est, coordinates = c(5,5), time = 40)
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -3.5, X1 = 1, X2 = .75, X3 = .5),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 1), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, region, spatial_family = "Gaussian", temporal_family = "Exponential", cov_map = example_background_covariates)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' plot_hawkes(hawkes)
#' plot_intensity(hawkes, est, time = 40, coordinates = c(4.5, 5))
plot_intensity <- function(hawkes, est, time = NULL, coordinates = NULL, step = 0.1) {
  region <- attr(hawkes, "region")
  x_min <- region$x[1]
  x_max <- region$x[2]
  y_min <- region$y[1]
  y_max <- region$y[2]
  t_min <- region$t[1]
  t_max <- region$t[2]

  plots <- list()

  if (!is.null(time)) {
    spatial <- spatial_conditional_intensity(hawkes, est, time, step = step)

    plots$spatial <- spatial |>
      ggplot2::ggplot() +
      ggplot2::geom_raster(ggplot2::aes(x, y, fill = intensity)) +
      ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, NA)) +
      ggplot2::geom_point(data = dplyr::filter(hawkes, t < time), ggplot2::aes(x, y)) +
      ggplot2::labs(x = "X", y = "Y", fill = "Intensity",
                    title = paste("Spatial Intensity at t =", time))
  }

  if (!is.null(coordinates)) {
    temporal <- temporal_conditional_intensity(hawkes, est, coordinates, step = step)

    plots$temporal <- temporal |>
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(t, intensity)) +
      ggplot2::labs(x = "Time", y = "Conditional Intensity",
                    title = paste0("Temporal Intensity at (", coordinates[1], ", ", coordinates[2], ")"))
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}
