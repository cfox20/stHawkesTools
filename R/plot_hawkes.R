
#' Plot Hawkes Process Points
#'
#' Creates a scatterplot of points from a spatio-temporal Hawkes process, colored either by time or by whether they are background events.
#'
#' @param hawkes A `hawkes` object
#' @param color A string, either `"time"` to color by time, or `"background"` to color by whether an event is a background event (only works for simulated processes).
#' @param ... Additional arguments passed to `ggplot2::geom_sf(data = hawkes, ...)`. For full customization directly use ggplot() + geom_sf()
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region, spatial_burnin = 1)
#'
#' plot_hawkes(hawkes, color = "time")
#'
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#'
#' plot_hawkes(hawkes, color = "time")
plot_hawkes <- function(hawkes, color = "time",...) {
  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  time_window <- attrs$time_window
  spatial_region <- attrs$spatial_region
  covariate_columns    <- attrs$covariate_columns
  spatial_family    <- attrs$spatial_family
  temporal_family    <- attrs$temporal_family
  spatial_sampler    <- attrs$spatial_sampler
  temporal_sampler    <- attrs$temporal_sampler
  spatial_pdf  <- attrs$spatial_pdf
  temporal_pdf <- attrs$temporal_pdf
  spatial_cdf  <- attrs$spatial_cdf
  temporal_cdf <- attrs$temporal_cdf
  spatial_is_separable <- isTRUE(attrs$spatial_is_separable)


  if (color == "time") {
    plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = spatial_region) +
      ggplot2::geom_sf(data = hawkes, ggplot2::aes(color = .data$t), ...) +
      ggplot2::labs(color = "Time")
  }
  if (color == "background") {
    plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = spatial_region) +
      ggplot2::geom_sf(data = hawkes |>
                         dplyr::mutate(background = factor(.data$gen == 0, levels = c("TRUE", "FALSE"))),
                       ggplot2::aes(color = .data$background), ...) +
      ggplot2::labs(color = "Background Event")
  }
  plot
}


#' Plot Conditional Intensity of Hawkes Process
#'
#' Visualizes the conditional intensity function at a fixed time over a spatial grid, based on a fitted Hawkes model.
#'
#' @param hawkes A `hawkes` object
#' @param est A list of fitted parameter values for the Hawkes process.
#' @param stepsize A numeric value for the spatial resolution of the evaluation grid.
#' @param time A numeric value giving the time at which to evaluate the conditional intensity.
#' @param coordinates A numeric vector of length 2 specifying the coordinates where to evalute the conditional intensity over time.
#'
#' @return A list of `ggplot2` objects showing the spatial intensity and temporal intensity at a given time or location.
#' @export
#'
#' @examples
#' # example code
#' set.seed(123)
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' plot_intensity(hawkes, est, stepsize = .25, time = 40)
#' plot_intensity(hawkes, est, stepsize = .25, coordinates = c(5,5))
#' plot_intensity(hawkes, est, stepsize = .25, coordinates = c(5,5), time = 40)
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' plot_hawkes(hawkes)
#' plot_intensity(hawkes, est, stepsize = .1, time = 40, coordinates = c(4.5, 5))
plot_intensity <- function(hawkes, est, stepsize, time = NULL, coordinates = NULL) {
  plots <- list()

  if (!is.null(time)) {
    spatial <- spatial_conditional_intensity(hawkes, est, time, stepsize)

    plots$spatial <- spatial |>
      ggplot2::ggplot() +
      ggplot2::geom_raster(ggplot2::aes(.data$x, .data$y, fill = .data$intensity)) +
      ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, NA)) +
      ggplot2::geom_point(data = dplyr::filter(hawkes, .data$t < time), ggplot2::aes(.data$x, .data$y)) +
      ggplot2::labs(x = "X", y = "Y", fill = "Intensity",
                    title = paste("Spatial Intensity at t =", time))
  }

  if (!is.null(coordinates)) {
    temporal <- temporal_conditional_intensity(hawkes, est, coordinates, stepsize)

    plots$temporal <- temporal |>
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(.data$t, .data$intensity)) +
      ggplot2::labs(x = "Time", y = "Conditional Intensity",
                    title = paste0("Temporal Intensity at (", coordinates[1], ", ", coordinates[2], ")"))
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}
