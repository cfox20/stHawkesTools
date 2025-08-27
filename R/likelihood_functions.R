#' Compute Conditional Intensity at Each Event
#'
#' @param hawkes A hawkes object
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, region)
#' conditional_intensity(hawkes, params)
conditional_intensity <- function(hawkes, parameters) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  if (length(background_rate) == 1) {
    X <- 1
  } else{
    if (dim(X)[2] - length(background_rate) == -1) {
      X <- cbind(1,X)
    }
  }


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
    g_mat <- {triggering_rate *
        do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
        do.call(spatial_pdf, c(list(x = x_diff), parameters$spatial)) *
        do.call(spatial_pdf, c(list(x = y_diff), parameters$spatial))
    }
  }
  g_mat[upper.tri(g_mat, diag = TRUE)] <- 0

  # Store the values of the complete likelihood at each point
  exp(as.numeric(X %*% background_rate)) + rowSums(g_mat)
}


#' Compute the Conditional Intensity at a Specified Time t
#'
#' @param hawkes A `hawkes` object.
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#' @param time A time to evaluate the conditional intensity.
#' @param step A numeric value specifying the size of the grid to evaluate the conditional intensity over.
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, region)
#' spatial_conditional_intensity(hawkes, params, 25)
spatial_conditional_intensity <- function(hawkes, parameters, time, step = .1) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  hawkes <- hawkes[hawkes$t < time,]

  region$t[2] <- time

  grid <- expand.grid(seq(region$x[1], region$x[2], step), seq(region$y[1], region$y[2], step))
  x <- grid[,1]
  y <- grid[,2]

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  x_diff <- outer(x, hawkes$x, `-`)
  # x_diff[upper.tri(x_diff)] <- 0
  y_diff <- outer(y, hawkes$y, `-`)
  # y_diff[upper.tri(y_diff)] <- 0

  time_diff <- matrix(rep(outer(time, hawkes$t, `-`), nrow(x_diff)), nrow = nrow(x_diff), byrow = TRUE)
  # time_diff[upper.tri(time_diff)] <- 0

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


  if (length(background_rate) == 1) {
    X <- 1
  } else{
    cov_names <- colnames(cov_map)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]

    X <- sf::st_as_sf(grid, coords = c("Var1", "Var2"), crs = sf::st_crs(cov_map)) |>
      sf::st_join(cov_map, join = sf::st_intersects) |>
      dplyr::mutate(.wkt = sf::st_as_text(geometry)) |>
      dplyr::distinct(.wkt, .keep_all = TRUE) |>
      dplyr::select(-.wkt) |>
      sf::st_drop_geometry() |>
      dplyr::select(all_of(cov_names)) |>
      as.matrix()
    X <- cbind(1, X)
  }

  # Store the values of the complete likelihood at each point
  # exp(as.numeric(X %*% background_rate)) + rowSums(g_mat)
  data.frame(x = x,
             y = y,
             t = time,
             intensity = exp(as.numeric(X %*% background_rate)) + rowSums(g_mat))
}


#' Compute the Conditional Intensity at a Specified Location (x,y)
#'
#' @param hawkes A `hawkes` object.
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#' @param coordinates A numeric vector of length 2 to specify the location to evaluate the conditional intensity.
#' @param step A numeric value specifying the size of the grid to evaluate the conditional intensity over.
#'
#' @returns A numeric vector
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, region)
#' temporal_conditional_intensity(hawkes, params, c(5,5))
temporal_conditional_intensity <- function(hawkes, parameters, coordinates, step = .1) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  x <- coordinates[1]
  y <- coordinates[2]

  t_points <- hawkes$t
  t_grid <- c(t_points, seq(region$t[1], region$t[2], by = step)) |>
    sort()

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

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
  # g_mat[upper.tri(g_mat, diag = TRUE)] <- 0


  if (length(background_rate) == 1) {
    X <- 1
  } else{
    cov_names <- colnames(cov_map)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]

    X <- sf::st_as_sf(data.frame(t(coordinates)), coords = c("X1", "X2"), crs = sf::st_crs(cov_map)) |>
      sf::st_join(cov_map, join = sf::st_intersects) |>
      dplyr::mutate(.wkt = sf::st_as_text(geometry)) |>
      dplyr::distinct(.wkt, .keep_all = TRUE) |>
      dplyr::select(-.wkt) |>
      sf::st_drop_geometry() |>
      dplyr::select(all_of(cov_names)) |>
      as.matrix()
    X <- cbind(1, X)
  }

  # Store the values of the complete likelihood at each point
  # Update this to handle background covariates
  data.frame(x = x,
             y = y,
             t = t_grid,
             intensity = exp(as.numeric(X %*% background_rate)) + rowSums(g_mat))
}


#' Compute Full Log-Likelihood
#'
#' @param hawkes A hawkes object
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#'
#' @returns A numeric vector.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, region)
#' full_log_likelihood(hawkes, params)
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(min = -1, max = 1),temporal = list(shape = 3, rate = 1))
#' hawkes <- rHawkes(params, region, spatial_family = "Uniform", temporal_family = "Gamma")
#' full_log_likelihood(hawkes, params)
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(rate = 5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, region, spatial_family = "Exponential", temporal_family = "Exponential")
#' full_log_likelihood(hawkes, params)
full_log_likelihood <- function(hawkes, parameters) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  # Region bounds
  x_min <- region$x[1]
  x_max <- region$x[2]
  y_min <- region$y[1]
  y_max <- region$y[2]
  t_min <- region$t[1]
  t_max <- region$t[2]

  # Conditional intensity at observed events
  log_lambda <- log(conditional_intensity(hawkes, parameters))

  log_lambda[is.infinite(log_lambda)] <- -1e10  # numerical safeguard
  log_part <- sum(log_lambda)

  # Background integral (spatial + temporal)
  if (exists("cov_map", inherits = FALSE)) {
    cov_names <- colnames(cov_map)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]

    cov_map_X <- cov_map |>
      sf::st_drop_geometry() |>
      dplyr::select(cov_names) |>
      as.matrix()
    cov_map_X <- cbind(1, cov_map_X)
    area <- cov_map$area

    background_rates <- area * exp(cov_map_X %*% background_rate)
    background_integral <- t_max * sum(background_rates)
  } else {
    # No covariates — scalar background
    background_integral <- (t_max - t_min) * (x_max - x_min) * (y_max - y_min) * exp(background_rate)
  }


  if (!spatial_is_separable) {
    upper <- cbind(region$x[2] - hawkes$x, region$y[2] - hawkes$y)
    lower <- cbind(region$x[1] - hawkes$x, region$y[1] - hawkes$y)

    spatial_mass <- do.call(spatial_cdf, c(list(q = upper), spatial_params)) -
                    do.call(spatial_cdf, c(list(q = lower), spatial_params))

  } else {
    spatial_mass <- (
      do.call(spatial_cdf, c(list(q = region$x[2] - hawkes$x), spatial_params)) -
      do.call(spatial_cdf, c(list(q = region$x[1] - hawkes$x), spatial_params))
    ) * (
      do.call(spatial_cdf, c(list(q = region$y[2] - hawkes$y), spatial_params)) -
      do.call(spatial_cdf, c(list(q = region$y[1] - hawkes$y), spatial_params))
    )
  }

  triggering_integral <- triggering_rate * sum(
    do.call(temporal_cdf, c(list(q = region$t[2] - hawkes$t), temporal_params)) *
    spatial_mass
  )

  return(log_part - background_integral - triggering_integral)
}





