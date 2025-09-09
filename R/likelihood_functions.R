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

  if(!exists("covariate_columns", inherits = FALSE)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry() |> as.matrix())
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
  as.numeric(exp(as.numeric(X %*% background_rate)) + rowSums(g_mat))
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
    X <- sf::st_as_sf(grid, coords = c("Var1", "Var2"), crs = sf::st_crs(spatial_region)) |>
      sf::st_join(spatial_region, join = sf::st_intersects) |>
      dplyr::mutate(.wkt = sf::st_as_text(geometry)) |>
      dplyr::distinct(.wkt, .keep_all = TRUE) |>
      dplyr::select(-.wkt) |>
      sf::st_drop_geometry() |>
      dplyr::select(all_of(covariate_columns)) |>
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
    X <- sf::st_as_sf(data.frame(t(coordinates)), coords = c("X1", "X2"), crs = sf::st_crs(spatial_region)) |>
      sf::st_join(spatial_region, join = sf::st_intersects) |>
      dplyr::mutate(.wkt = sf::st_as_text(geometry)) |>
      dplyr::distinct(.wkt, .keep_all = TRUE) |>
      dplyr::select(-.wkt) |>
      sf::st_drop_geometry() |>
      dplyr::select(all_of(covariate_columns)) |>
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


#' Compute Log-Likelihood
#'
#' @param hawkes A hawkes object
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#'
#' @returns A numeric vector.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, c(0,50), spatial_region)
#' log_likelihood(hawkes, params)
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 0)
#' log_likelihood(hawkes, params)
log_likelihood <- function(hawkes, parameters) {
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

  time_length <- time_window[2] - time_window[1]

  # Conditional intensity at observed events
  log_lambda <- log(conditional_intensity(hawkes, parameters))

  log_lambda[is.infinite(log_lambda)] <- -1e10  # numerical safeguard
  log_part <- sum(log_lambda)

  # Background integral (spatial + temporal)
  if (exists("covariate_columns", inherits = FALSE)) {
    covariate_map <- spatial_region |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(covariate_columns)) |>
      as.matrix()
    covariate_map <- cbind(1, covariate_map)
    area <- spatial_region$area

    background_rates <- area * exp(covariate_map %*% background_rate)
    background_integral <- time_length * sum(background_rates)
  } else {
    # No covariates — scalar background
    background_integral <- time_length * sf::st_area(spatial_region) * exp(background_rate)
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
  triggering_integral <- triggering_rate * sum(
    do.call(temporal_cdf, c(list(q = time_window[2] - hawkes$t), temporal_params))
    # spatial_mass
  )

  return(log_part - background_integral - triggering_integral)
}





