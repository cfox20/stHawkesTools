#' Compute Conditional Intensity at Each Event
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters. Values usually represent current estimates.
#'
#' @returns A numeric vector.
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
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' conditional_intensity(hawkes, params)
#'
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4,
#'                          event_type = c(a = 1, b = .25, c = .5)),
#'   triggering_rate = matrix(c(.4, .15, .05,
#'                              .2, .05, .02,
#'                              .2, .05, .20),
#'                            nrow = 3,
#'                            dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
#'   spatial = list(mean = 0, sd = 0.1),
#'   temporal = list(rate = 2)
#' )
#' (hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1 + mark(event_type),
#'   spatial_burnin = 1
#' ))
#'
#' conditional_intensity(hawkes, params)
conditional_intensity <- function(hawkes, parameters) {
  if(!inherits(hawkes, "hawkes")) stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  covariate_matrix <- attrs$covariate_matrix
  branching_matrix <- attrs$branching_matrix
  background_formula <- attrs$background_formula
  mark_column <- attrs$mark_column
  time_window <- attrs$time_window
  spatial_region <- attrs$spatial_region
  parameters <- attrs$parameters
  spatial_family <- attrs$spatial_family
  temporal_family <- attrs$temporal_family
  spatial_sampler <- attrs$spatial_sampler
  spatial_pdf <- attrs$spatial_pdf
  spatial_cdf <- attrs$spatial_cdf
  temporal_pdf <- attrs$temporal_pdf
  temporal_cdf <- attrs$temporal_cdf
  temporal_sampler <- attrs$temporal_sampler
  spatial_is_separable <- attrs$spatial_is_separable


  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  mark_effects <- NULL
  if (!is.null(mark_column) && length(mark_column) > 0 &&
      !is.null(background_rate[[mark_column]])) {
    mark_effects <- background_rate[[mark_column]]
    background_rate[[mark_column]] <- NULL

    if (!is.null(mark_effects)) {
      mark_effect_names <- names(mark_effects)
      mark_effects <- as.numeric(mark_effects)
      if (!is.null(mark_effect_names)) {
        names(mark_effects) <- mark_effect_names
      }
    }
  }

  if (!is.null(mark_effects)) {
    background_mark_offset <- mark_effects[hawkes$event_type]
  } else{
    background_mark_offset <- 0
  }

  if (length(background_rate) > 0) {
    background_rate <- as.numeric(background_rate)
  } else {
    background_rate <- numeric(0)
  }

  if(!exists("covariate_columns", inherits = FALSE)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry() |> as.matrix())
  }

  if (is.matrix(triggering_rate)){
    triggering_matrix <- triggering_rate[hawkes[[mark_column]], hawkes[[mark_column]], drop = FALSE]
  } else {
    triggering_matrix <- matrix(triggering_rate, nrow = nrow(hawkes), ncol = nrow(hawkes))
  }
  triggering_matrix[upper.tri(triggering_matrix, diag = TRUE)] <- 0


# Compute Intensity -------------------------------------------------------

  x_diff <- outer(hawkes$x, hawkes$x, `-`)
  x_diff[upper.tri(x_diff)] <- 0
  y_diff <- outer(hawkes$y, hawkes$y, `-`)
  y_diff[upper.tri(y_diff)] <- 0

  time_diff <- outer(hawkes$t, hawkes$t, `-`)
  time_diff[upper.tri(time_diff)] <- 0


  # Compute a matrix of the triggering intensities
  if (!spatial_is_separable) {
    s_diff <- cbind(x_diff, y_diff)
    g_mat <- {triggering_rate *
        do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
        do.call(spatial_pdf, c(list(x = s_diff), parameters$spatial))
    }
  } else {
    g_mat <- {triggering_matrix *
        do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
        do.call(spatial_pdf, c(list(x = x_diff), parameters$spatial)) *
        do.call(spatial_pdf, c(list(x = y_diff), parameters$spatial))
    }
  }
  g_mat[upper.tri(g_mat, diag = TRUE)] <- 0

  # Store the values of the complete likelihood at each point
  as.numeric(exp(as.numeric(X %*% background_rate) + background_mark_offset) + rowSums(g_mat))
}


#' Compute the conditional intensity at a specified time
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters.
#' @param time Time point where the conditional intensity is evaluated.
#' @param stepsize Grid cell size used to evaluate the spatial intensity.
#' @param spatial_zoom Optional `sf` object representing a subset of the spatial region on which to
#'   evaluate the intensity.
#'
#' @returns A numeric vector.
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
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' spatial_conditional_intensity(hawkes, params, 25, 0.5)
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
#' spatial_conditional_intensity(hawkes, params, 25, 0.5)
spatial_conditional_intensity <- function(hawkes, parameters, time, stepsize, spatial_zoom = NULL) {
  if(!inherits(hawkes, "hawkes")) stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

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


  hawkes <- hawkes[hawkes$t < time,]

  if (!is.null(spatial_zoom)) {
    spatial_region <- spatial_zoom

    suppressWarnings({
      spatial_filter <- sf::st_intersects(hawkes, spatial_region, sparse = FALSE)
    })

    if (length(spatial_filter)) {
      hawkes <- hawkes[apply(spatial_filter, 1, any), ]
    } else {
      hawkes <- hawkes[FALSE, ]
    }
  }

  time_window[2] <- time

  point_grid <- spatial_region |>
    sf::st_make_grid(what = "centers", cellsize = stepsize) |>
    sf::st_as_sf() |>
    sf::st_intersection(spatial_region) |>
    suppressWarnings()

  if (nrow(point_grid) == 0) {
    return(data.frame(x = numeric(), y = numeric(), intensity = numeric()))
  }

  x <- sf::st_coordinates(point_grid)[,1]
  y <- sf::st_coordinates(point_grid)[,2]

  X <- point_grid |>
    sf::st_drop_geometry() |>
    dplyr::select(tidyselect::all_of(covariate_columns))

  X <- cbind(1,X) |>
    as.matrix()

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  x_diff <- outer(x, hawkes$x, `-`)
  y_diff <- outer(y, hawkes$y, `-`)

  time_diff <- matrix(rep(outer(time, hawkes$t, `-`), nrow(x_diff)), nrow = nrow(x_diff), byrow = TRUE)

  # Compute a matrix of the triggering intensities
  if (!spatial_is_separable) {
    s_diff <- cbind(x_diff, y_diff)
    g_mat <- {triggering_rate *
             do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
             do.call(spatial_pdf, c(list(x = s_diff), parameters$spatial))
             }
  } else {
    g_mat <- {triggering_rate *
             do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
             do.call(spatial_pdf, c(list(x = x_diff), parameters$spatial)) *
             do.call(spatial_pdf, c(list(x = y_diff), parameters$spatial))
             }
  }
  # g_mat[upper.tri(g_mat, diag = TRUE)] <- 0


  # Store the values of the complete likelihood at each point
  # exp(as.numeric(X %*% background_rate)) + rowSums(g_mat)
  sf::st_geometry(point_grid) <- "geometry"
  point_grid$x <-  x
  point_grid$y <-  y
  point_grid$intensity <- exp(as.numeric(X %*% background_rate)) + rowSums(g_mat)
  point_grid$background <- exp(as.numeric(X %*% background_rate))
  point_grid$triggering <- rowSums(g_mat)


  point_grid[, c(c("x", "y", "intensity"), setdiff(names(point_grid), c("x", "y", "intensity")))]
}


#' Compute the conditional intensity at a location
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters.
#' @param coordinates Numeric vector of length two giving the evaluation location.
#' @param step Grid step size used to build the temporal evaluation grid.
#' @param time_window Optional numeric vector of length two giving the time interval over which to
#'   evaluate the intensity.
#'
#' @returns A numeric vector.
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
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' temporal_conditional_intensity(hawkes, params, c(5, 5))
temporal_conditional_intensity <- function(hawkes, parameters, coordinates, step = .1,
                                           time_window = NULL) {
  if(!inherits(hawkes, "hawkes")) stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  time_window_attr <- attrs$time_window
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


  x <- coordinates[1]
  y <- coordinates[2]

  if (is.null(time_window)) {
    time_window <- time_window_attr
  } else {
    if (length(time_window) != 2 || !is.numeric(time_window)) {
      stop("time_window must be a numeric vector of length 2.")
    }
    time_window <- sort(time_window)
    time_window[1] <- max(time_window_attr[1], time_window[1])
    time_window[2] <- min(time_window_attr[2], time_window[2])
  }

  if (time_window[1] > time_window[2]) {
    stop("time_window must overlap the observed time interval.")
  }

  hawkes <- hawkes[hawkes$t >= time_window[1] & hawkes$t <= time_window[2], ]

  t_points <- hawkes$t
  t_grid <- c(t_points, seq(time_window[1], time_window[2], by = step)) |>
    unique() |>
    sort()

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  if (nrow(hawkes) > 0) {
    time_diff <- outer(t_grid, hawkes$t, `-`)
    x_diff <- matrix(rep(outer(x, hawkes$x, `-`), nrow(time_diff)), nrow = nrow(time_diff), byrow = TRUE)
    y_diff <- matrix(rep(outer(y, hawkes$y, `-`), nrow(time_diff)), nrow = nrow(time_diff), byrow = TRUE)

    # Compute a matrix of the triggering intensities
    if (!spatial_is_separable) {
      s_diff <- cbind(x_diff, y_diff)
      g_mat <- {triggering_rate *
               do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
               do.call(spatial_pdf, c(list(x = s_diff), parameters$spatial))
               }
    } else {
      g_mat <- {triggering_rate *
               do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
               do.call(spatial_pdf, c(list(x = x_diff), parameters$spatial)) *
               do.call(spatial_pdf, c(list(x = y_diff), parameters$spatial))
               }
    }
  } else {
    g_mat <- matrix(0, nrow = length(t_grid), ncol = 0)
  }
  # g_mat[upper.tri(g_mat, diag = TRUE)] <- 0


  if (length(background_rate) == 1) {
    X <- 1
  } else{
    X <- sf::st_as_sf(data.frame(t(coordinates)), coords = c("X1", "X2"), crs = sf::st_crs(spatial_region)) |>
      sf::st_join(spatial_region, join = sf::st_intersects) |>
      dplyr::mutate(.wkt = sf::st_as_text(.data$geometry)) |>
      dplyr::distinct(.data$.wkt, .keep_all = TRUE) |>
      dplyr::select(-.data$.wkt) |>
      sf::st_drop_geometry() |>
      dplyr::select(tidyselect::all_of(covariate_columns)) |>
      as.matrix()
    X <- cbind(1, X)
  }

  # Store the values of the complete likelihood at each point
  # Update this to handle background covariates
  data.frame(x = x,
             y = y,
             t = t_grid,
             intensity = exp(as.numeric(X %*%
                                          background_rate)) + rowSums(g_mat))
}


#' Compute log-likelihood
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters.
#'
#' @returns A numeric vector.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.5),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' log_likelihood(hawkes, params)
#'
#' params <- list(
#'   background_rate = list(intercept = -4,
#'                          event_type = c(a = 1, b = .25, c = .5)),
#'   triggering_rate = matrix(c(.4, .15, .05,
#'                              .2, .05, .02,
#'                              .2, .05, .20),
#'                            nrow = 3,
#'                            dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
#'   spatial = list(mean = 0, sd = 0.1),
#'   temporal = list(rate = 2)
#' )
#' (hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1 + mark(event_type),
#'   spatial_burnin = 1
#' ))
#'
#' log_likelihood(hawkes, params)
log_likelihood <- function(hawkes, parameters) {
  if(!inherits(hawkes, "hawkes")) stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  time_window <- attrs$time_window
  spatial_region <- attrs$spatial_region
  covariate_columns    <- attrs$covariate_columns
  mark_column <- attrs$mark_column
  spatial_family    <- attrs$spatial_family
  temporal_family    <- attrs$temporal_family
  spatial_sampler    <- attrs$spatial_sampler
  temporal_sampler    <- attrs$temporal_sampler
  spatial_pdf  <- attrs$spatial_pdf
  temporal_pdf <- attrs$temporal_pdf
  spatial_cdf  <- attrs$spatial_cdf
  temporal_cdf <- attrs$temporal_cdf
  spatial_is_separable <- isTRUE(attrs$spatial_is_separable)


  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  mark_effects <- NULL
  if (!is.null(mark_column) && length(mark_column) > 0 &&
      !is.null(background_rate[[mark_column]])) {
    mark_effects <- background_rate[[mark_column]]
    background_rate[[mark_column]] <- NULL

    if (!is.null(mark_effects)) {
      mark_effect_names <- names(mark_effects)
      mark_effects <- as.numeric(mark_effects)
      if (!is.null(mark_effect_names)) {
        names(mark_effects) <- mark_effect_names
      }
    }
  } else{
    mark_effects <- 0
  }

  if (length(background_rate) > 0) {
    background_rate <- as.numeric(background_rate)
  } else {
    background_rate <- numeric(0)
  }

  if (is.matrix(triggering_rate)){
    triggering_matrix <- triggering_rate[hawkes$event_type, hawkes$event_type, drop = FALSE]
  } else {
    triggering_matrix <- matrix(triggering_rate, nrow = nrow(hawkes), ncol = nrow(hawkes))
  }
  triggering_matrix[upper.tri(triggering_matrix, diag = TRUE)] <- 0

  time_length <- time_window[2] - time_window[1]

  # Conditional intensity at observed events
  log_lambda <- log(conditional_intensity(hawkes, parameters))

  log_lambda[is.infinite(log_lambda)] <- -1e10  # numerical safeguard
  log_part <- sum(log_lambda)



  # Background integral (spatial + temporal)
  if(!is.null(covariate_columns)){
    covariate_map <- spatial_region |>
      sf::st_drop_geometry() |>
      dplyr::select(tidyselect::all_of(covariate_columns)) |>
      as.matrix()
    covariate_map <- cbind(1, covariate_map)
    area <- spatial_region$area |> as.numeric()

    background_rates <- area * exp(covariate_map %*% background_rate + mark_effects)
    background_integral <- time_length * sum(background_rates)
  } else {
    # No covariates — scalar background
    background_integral <- time_length * as.numeric(sf::st_area(spatial_region)) * sum(exp(background_rate + mark_effects))
  }


  # if (!spatial_is_separable) {
  #   upper <- cbind(region$x[2] - hawkes$x, region$y[2] - hawkes$y)
  #   lower <- cbind(region$x[1] - hawkes$x, region$y[1] - hawkes$y)
  #
  #   spatial_mass <- do.call(spatial_cdf, c(list(q = upper), spatial_params)) -
  #                   do.call(spatial_cdf, c(list(q = lower), spatial_params))
  #
  # } else {
  #   spatial_mass <- (
  #     do.call(spatial_cdf, c(list(q = region$x[2] - hawkes$x), spatial_params)) -
  #     do.call(spatial_cdf, c(list(q = region$x[1] - hawkes$x), spatial_params))
  #   ) * (
  #     do.call(spatial_cdf, c(list(q = region$y[2] - hawkes$y), spatial_params)) -
  #     do.call(spatial_cdf, c(list(q = region$y[1] - hawkes$y), spatial_params))
  #   )
  # }
  # This is an approximation by assuming the the integral over the spatial region integrates to 1.
  # This is done to simplify the computation for non-square spatial regions
  # and should only cause problems if there are many events near the edge of the observed spatial region.
  triggering_integral <- sum(triggering_rate * sum(
    do.call(temporal_cdf, c(list(q = time_window[2] - hawkes$t), temporal_params))
    # spatial_mass
  ))

  return(log_part - background_integral - triggering_integral)
}





