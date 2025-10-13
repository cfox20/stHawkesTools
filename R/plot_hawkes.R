
#' Plot Hawkes Process Points
#'
#' Creates a scatter plot of Hawkes events colored by time or background status.
#'
#' @param hawkes A `hawkes` object
#' @param color Either `"time"` to color by event time or `"background"` to show
#'   simulated background events.
#' @param ... Additional arguments forwarded to `ggplot2::geom_sf()`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.75,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(
#'   background_process = ~ 1,
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   spatial_burnin = 1
#' )
#'
#' plot_hawkes(hawkes, color = "time")
#'
#'
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = example_background_covariates,
#'   background_process = ~ X1 + X2,
#'   spatial_burnin = 1
#' )
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


#' Plot conditional intensity of a Hawkes process
#'
#' Visualizes the conditional intensity over space or time given a fitted Hawkes model.
#'
#' @param hawkes A `hawkes` object
#' @param est A list of fitted parameter values for the Hawkes process.
#' @param stepsize A numeric value for the spatial resolution of the evaluation grid.
#' @param time A numeric value giving the time at which to evaluate the conditional intensity.
#' @param coordinates Numeric vector of length two giving the evaluation location.
#' @param interpolate If TRUE interpolate linearly, if FALSE (the default) don't interpolate.
#' @param intensity_type Name of type of intensity to plot past as un unquoted name. Defaults to intensity for full intensity. Other options include triggering for triggering intensity and background for background intensity.
#' @param point_size Numeric to set the size of points for the spatial plot
#' @param spatial_zoom Optional numeric vector of length four giving xmin, xmax, ymin, and ymax
#'   bounds for a spatial zoom window.
#' @param temporal_zoom Optional numeric vector of length two giving the start and end times for a
#'   temporal zoom window.
#' @param zoom_to_most_recent Logical; if `TRUE`, overrides the manual zoom settings and centers the
#'   view on the most recent event observed prior to `time`.
#' @param recent_spatial_radius Optional numeric vector of length one or two giving the half-width of
#'   the spatial window used when `zoom_to_most_recent = TRUE`.
#' @param recent_time_window Optional numeric vector of length one or two giving the amount of time
#'   before and after the most recent event to include when `zoom_to_most_recent = TRUE`.
#'
#' @return A list of `ggplot2` objects showing spatial and temporal intensity views.
#' @export
#'
#' @examples
#' # example code
#' set.seed(123)
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.75,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' plot_intensity(hawkes, est, stepsize = 0.25, time = 40)
#' plot_intensity(hawkes, est, stepsize = 0.25, coordinates = c(5, 5))
#' plot_intensity(hawkes, est, stepsize = 0.25, coordinates = c(5, 5), time = 40)
#'
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.15,
#'   spatial = list(mean = 0, sd = 0.25),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = example_background_covariates,
#'   background_process = ~ X1 + X2,
#'   spatial_burnin = 1
#' )
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' plot_hawkes(hawkes)
#' plot_intensity(hawkes, est, stepsize = .05, time = 50, coordinates = c(4.5, 5))
#' plot_intensity(hawkes, est, stepsize = .025, time = 32, interpolate = TRUE, intensity_type = triggering)
plot_intensity <- function(hawkes, est, stepsize, time = NULL, coordinates = NULL,
                           interpolate = FALSE, intensity_type = intensity,
                           point_size = 1.5, spatial_zoom = NULL,
                           temporal_zoom = NULL, zoom_to_most_recent = FALSE,
                           recent_spatial_radius = NULL, recent_time_window = NULL) {
  if (is.null(time) && is.null(coordinates)) {
    stop("At least 1 of time or coordinates must be provided.")
  }

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

  plots <- list()

  zoom_region <- NULL
  time_limits <- NULL

  if (!is.null(spatial_zoom)) {
    if (!(is.numeric(spatial_zoom) && length(spatial_zoom) == 4)) {
      stop("spatial_zoom must be a numeric vector of length 4 specifying xmin, xmax, ymin, and ymax.")
    }
    zoom_bbox <- sf::st_bbox(c(xmin = spatial_zoom[1], xmax = spatial_zoom[2],
                               ymin = spatial_zoom[3], ymax = spatial_zoom[4]),
                             crs = sf::st_crs(spatial_region))
    zoom_region <- suppressWarnings(sf::st_crop(spatial_region, zoom_bbox))
  }

  if (!is.null(temporal_zoom)) {
    if (!(is.numeric(temporal_zoom) && length(temporal_zoom) == 2)) {
      stop("temporal_zoom must be a numeric vector of length 2 specifying the time window.")
    }
    time_limits <- sort(temporal_zoom)
    time_limits[1] <- max(time_window[1], time_limits[1])
    time_limits[2] <- min(time_window[2], time_limits[2])
  }

  if (zoom_to_most_recent) {
    if (is.null(time)) {
      stop("time must be provided when zoom_to_most_recent = TRUE.")
    }

    hawkes_prior <- hawkes[hawkes$t <= time, ]
    if (nrow(hawkes_prior) == 0) {
      stop("No events observed prior to the provided time. Unable to determine most recent point.")
    }

    recent_point <- hawkes_prior[which.max(hawkes_prior$t), ]

    if (is.null(recent_spatial_radius)) {
      recent_spatial_radius <- rep(stepsize * 50, 2)
    }
    if (length(recent_spatial_radius) == 1) {
      recent_spatial_radius <- rep(recent_spatial_radius, 2)
    }
    if (!(is.numeric(recent_spatial_radius) && length(recent_spatial_radius) == 2)) {
      stop("recent_spatial_radius must be a numeric vector of length 1 or 2.")
    }

    spatial_zoom <- c(sf::st_coordinates(recent_point)[,1] - recent_spatial_radius[1],
                      sf::st_coordinates(recent_point)[,1] + recent_spatial_radius[1],
                      sf::st_coordinates(recent_point)[,2] - recent_spatial_radius[2],
                      sf::st_coordinates(recent_point)[,2] + recent_spatial_radius[2])

    zoom_bbox <- sf::st_bbox(c(xmin = spatial_zoom[[1]], xmax = spatial_zoom[[2]],
                               ymin = spatial_zoom[[3]], ymax = spatial_zoom[[4]]),
                             crs = sf::st_crs(spatial_region))
    zoom_region <- suppressWarnings(sf::st_crop(spatial_region, zoom_bbox))

    if (is.null(recent_time_window)) {
      recent_time_window <- c(10, 10)
    }
    if (length(recent_time_window) == 1) {
      recent_time_window <- rep(recent_time_window, 2)
    }
    if (!(is.numeric(recent_time_window) && length(recent_time_window) == 2)) {
      stop("recent_time_window must be numeric of length 1 or 2.")
    }

    temporal_zoom <- c(recent_point$t - recent_time_window[1],
                       recent_point$t + recent_time_window[2])
    temporal_zoom[1] <- max(time_window[1], temporal_zoom[1])
    temporal_zoom[2] <- min(time, temporal_zoom[2])

    time_limits <- temporal_zoom
  }

  if (!is.null(zoom_region) && nrow(zoom_region) == 0) {
    stop("The requested spatial zoom does not intersect the observed spatial region.")
  }

  if (!is.null(time_limits) && time_limits[1] >= time_limits[2]) {
    stop("The requested temporal zoom does not overlap with the observed time window.")
  }

  spatial_layer_region <- if (!is.null(zoom_region)) zoom_region else spatial_region
  zoom_bbox <- if (!is.null(zoom_region)) sf::st_bbox(zoom_region) else NULL

  if (!is.null(time)) {
    spatial <- spatial_conditional_intensity(hawkes, est, time, stepsize,
                                             spatial_zoom = spatial_layer_region)

    if (!is.null(zoom_bbox)) {
      spatial <- dplyr::filter(spatial, .data$x >= zoom_bbox[["xmin"]],
                               .data$x <= zoom_bbox[["xmax"]],
                               .data$y >= zoom_bbox[["ymin"]],
                               .data$y <= zoom_bbox[["ymax"]])
    }

    hawkes_points <- dplyr::filter(hawkes, .data$t < time)
    if (!is.null(zoom_region)) {
      suppressWarnings({
        pts_filter <- sf::st_intersects(hawkes_points, zoom_region, sparse = FALSE)
      })
      if (length(pts_filter)) {
        hawkes_points <- hawkes_points[apply(pts_filter, 1, any), ]
      } else {
        hawkes_points <- hawkes_points[FALSE, ]
      }
    }

    plots$spatial <- spatial |>
      ggplot2::ggplot() +
      ggplot2::geom_raster(ggplot2::aes(.data$x, .data$y, fill = {{intensity_type}}), interpolate = interpolate) +
      ggplot2::coord_sf() +
      ggplot2::scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, NA)) +
      ggplot2::geom_sf(data = spatial_layer_region, inherit.aes = FALSE, fill = NA) +
      ggplot2::geom_sf(data = hawkes_points, size = point_size) +
      ggplot2::labs(x = "X", y = "Y", fill = "Intensity",
                    title = paste("Conditional Intensity at t =", time))
  }

  if (!is.null(coordinates)) {
    temporal <- temporal_conditional_intensity(hawkes, est, coordinates, stepsize,
                                               time_window = time_limits)

    if (!is.null(time_limits)) {
      temporal <- dplyr::filter(temporal, .data$t >= time_limits[1], .data$t <= time_limits[2])
    }

    plots$temporal <- temporal |>
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(.data$t, .data$intensity)) +
      ggplot2::labs(x = "Time", y = "Conditional Intensity",
                    title = paste0("Conditional Intensity at (", coordinates[1], ", ", coordinates[2], ")"))
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  } else {
    return(plots)
  }
}
