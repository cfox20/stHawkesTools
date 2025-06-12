# Functions to generate a spatio-temporal Hawkes process
#


#' Create an empty hawkes object
#'
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list.
##'
#' @returns An empty hawkes object.
.blank_hawkes <- function(params = NULL, region = NULL) {
  new_hawkes(data = data.frame(x = numeric(), y = numeric(), t = numeric(), parent = numeric(), gen = numeric()),
             params = params,
             region = region,  # could be passed in or defaulted inside rhawkes
             spatial_kernel = NULL,
             temporal_kernel = NULL
  )
}


#' Simulate background events
#'
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param background_rate A vector of coefficients for the background covariates.
#' @param cov_map An sf object containing background covariates.
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A hawkes object with only generated background events.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#' background_rate <- -4
#'
#' sim_background_events(background_rate, region)
sim_background_events <- function(background_rate, region, cov_map = NULL) {
  x_length <- region$x[2] - region$x[1]
  y_length <- region$y[2] - region$y[1]
  t_length <- region$t[2] - region$t[1]

  if (length(background_rate) == 1) {
    num_events <- stats::rpois(1, exp(background_rate) * x_length * y_length * t_length)

    if (num_events == 0) {
      df <- tibble::tibble(
                       x = NULL,
                       y = NULL,
                       t = NULL,
                       id = NULL,
                       parent = NULL,
                       gen = NULL,
                       family = NULL)
    }
    df <- tibble::tibble(
                     x = runif(num_events, region$x[1], region$x[2]),
                     y = runif(num_events, region$y[1], region$y[2]),
                     t = runif(num_events, region$t[1], region$t[2]),
                     id = 1:num_events,
                     parent = 0,
                     gen = 0,
                     family = 1:num_events)

  } else{
    X <- cov_map |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::starts_with("X")) |>
      as.matrix()
    X <- cbind(1, X)

    num_events <- stats::rpois(nrow(X), exp(X %*% background_rate) * region$t[2] * cov_map$area)

    df <- sf::st_sample(cov_map, size = num_events, exact = TRUE) |>
      suppressWarnings() |>
      sf::st_coordinates() |>
      tibble::as_tibble() |>
      dplyr::rename(x = X,
             y = Y) |>
      dplyr::mutate(t = stats::runif(sum(num_events), region$t[1], region$t[2]),
             parent = 0,
             gen = 0) |>
      dplyr::ungroup() |>
      dplyr::mutate(family = dplyr::row_number(),
             id = family) |>
      dplyr::relocate(id)
  }

  as_hawkes(df, region = region)
}



#' Generate a Hawkes Process
#'
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param t_burnin Temporal burn-in for simulation. Defaults to 50 if not used.
#' @param s_burnin Spatial burn-in for simulation. Note that this may not work as expected with a background covariate map. Defaults to 0 if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
#' @param spatial_kernel A spatial triggering kernel function to generate data from. Defaults to rnorm if not used.
#' @param temporal_kernel A spatial triggering kernel function to generate data from. Defaults to rexp if not used.
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A hawkes object with a generated Hawkes process.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' rHawkes(params, region)
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(min = -1, max = 1),temporal = list(shape = 3, rate = 1))
#' rHawkes(params, region, spatial_kernel = runif, temporal_kernel = rgamma)
#'
#' params <- list(background_rate = -4,triggering_rate = 0.5,spatial = list(rate = 5),temporal = list(rate = 2))
#' rHawkes(params, region, spatial_kernel = rexp_spatial, temporal_kernel = rexp)
#'
#'
#' #params <- list(background_rate = c(-4.5, 1),triggering_rate = 0.5,spatial = list(rate = 5),temporal = list(rate = 2))
#' #data("example_background_covariates")
#' #example_background_covariates <- example_background_covariates[,c("geoid", "name", "area", "X1")]
#' #rHawkes(params, region, spatial_kernel = rexp_spatial, temporal_kernel = rexp, cov_map = example_background_covariates)
rHawkes <- function(params, region, t_burnin = 50, s_burnin = 0, cov_map = NULL, spatial_kernel = rnorm, temporal_kernel = rexp) {
  # Set burnin regions
  x_max <- region$x[2] + 2*s_burnin
  x_min <- region$x[1]
  y_max <- region$y[2] + 2*s_burnin
  y_min <- region$y[1]
  time_min <- region$t[1]
  time_max <- region$t[2] + t_burnin

  new_region <- list(x = c(x_min, x_max),
                     y = c(y_min, y_max),
                     t = c(time_min, time_max))

  background_rate <- params$background_rate
  triggering_rate <- params$triggering_rate
  temporal_params <- params$temporal
  spatial_params <- params$spatial

  if (!is.numeric(background_rate)) {
    stop("background_rate must be numeric vector stored within named params list.")
  }
  if (!((triggering_rate >= 0) && (triggering_rate < 1) && (is.numeric(triggering_rate)))) {
    stop("triggering_rate must be a numeric value between 0 and 1.")
  }

  # Check to see if covariates are included
  covariates <- !is.null(cov_map)


  # Generate background events
  data <- G <- sim_background_events(background_rate, new_region, cov_map = cov_map)

  # Specify generation l
  l <- 0

  while (TRUE) {
    O <- .blank_hawkes(params = params, region = region)
    l <- l+1


    N <- stats::rpois(nrow(G), triggering_rate)

    if(sum(N) == 0) {
      if (covariates) {
        # Ensure both are sf objects and in same CRS
        data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(cov_map))

        # Spatial join: assign covariates from polygon to each point
        x_values <- sf::st_join(data_sf, cov_map, join = sf::st_within) |>
          sf::st_drop_geometry() |>
          dplyr::select(id, dplyr::starts_with("X"))

        data <- dplyr::left_join(data, x_values, by = "id")
      }

      data <- data[(data$t > t_burnin) & (data$x > x_min + s_burnin) & (data$x < x_max - s_burnin) & (data$y > y_min + s_burnin) & (data$y < y_max - s_burnin),]
      data$t <- data$t - t_burnin
      data$x <- data$x - s_burnin
      data$y <- data$y - s_burnin

      data <- data[order(data$t),] |> as_hawkes(region = region, params = params)

      return(data)
    }

    for (i in 1:nrow(G)) {
      if (N[i] > 0) {
        # Model specified self-exciting kernels. Add parameters here
        spatial_result <- do.call(spatial_kernel, c(list(n = N[i]), params$spatial))

        if (is.data.frame(spatial_result) && all(c("x", "y") %in% names(spatial_result))) {
          # Case A: kernel returns a data.frame with x and y
          x <- spatial_result$x + G$x[i]
          y <- spatial_result$y + G$y[i]
        } else if (is.numeric(spatial_result) && length(spatial_result) == N[i]) {
          # Case B: kernel returns a numeric vector (e.g., rnorm)
          x <- spatial_result + G$x[i]
          y <- do.call(spatial_kernel, c(list(n = N[i]), params$spatial)) + G$y[i]
        } else {
          stop("Invalid return from spatial_kernel: must be either vector or data.frame with x and y")
        }
        t <- do.call(temporal_kernel, c(list(n = N[i]), params$temporal)) + G$t[i]

        parent <- G$id[i]
        family <- G$family[i]

        # Create dataframe of events in the new generation O
        O_i <- data.frame(x = x, y = y, t = t, parent = parent, gen = l, family = family)
        O <- rbind(O, O_i)
      }
    }
    O$id <- (1:nrow(O))+G$id[nrow(G)]

    # Filter events that lie outside of region or time bounds
    O <- O[,c("id", "x", "y", "t", "parent", "gen", "family")]
    O <- subset(O, (x < x_max) & (x > 0) & (y < y_max) & (y > 0) & (t < time_max))

    # Return data ordered by time if no events remain in new generation
    if (nrow(O) == 0) {
      if (covariates) {
        # Ensure both are sf objects and in same CRS
        data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(cov_map))

        # Spatial join: assign covariates from polygon to each point
        x_values <- sf::st_join(data_sf, cov_map, join = sf::st_within) |>
          sf::st_drop_geometry() |>
          dplyr::select(id, dplyr::starts_with("X"))

        data <- dplyr::left_join(data, x_values, by = "id")
      }

      data <- data[(data$t > t_burnin) & (data$x > x_min + s_burnin) & (data$x < x_max - s_burnin) & (data$y > y_min + s_burnin) & (data$y < y_max - s_burnin),]
      data$t <- data$t - t_burnin
      data$x <- data$x - s_burnin
      data$y <- data$y - s_burnin


      data <- data[order(data$t),] |> as_hawkes(region = region, params = params)

      return(data)
    }

    # Add new generation to data and reset current generation to previous new one
    data <- rbind(data, O)
    G <- O
  }
}






