
#' Create Rectangular SF Object
#'
#' @param xmin horizontal minimum of spatial region
#' @param xmax horizontal maximum of spatial region
#' @param ymin vertical minimum of spatial region
#' @param ymax vertical maximum of spatial region
#' @param crs coordinate reference system. Defaults to NA if unused.
#'
#' @returns an sf object with a rectangular region
#' @export
#'
#' @examples
#' create_rectangular_sf(0,10,0,10)
create_rectangular_sf <- function(xmin, xmax, ymin, ymax, crs = NA) {
  sf::st_as_sfc(sf::st_bbox(c(xmin = xmin,
                              ymin = ymin,
                              xmax = xmax,
                              ymax = ymax)), crs = crs) |>
    sf::st_as_sf()
}

#' Simulate background events
#'
#' @param background_rate A vector of coefficients for the background covariates.
#' @param time_window A numeric vector of length 2 specifying the simulated time window.
#' @param spatial_region An sf object defining the spatial region for simulation.
#' @param covariate_columns A character vector of the names of the columns in spatial_region to be used as background covariates.
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A dataframe with generated background events.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#' time_window <- c(0,50)
#' background_rate <- -4
#'
#' sim_background_events(background_rate, time_window, spatial_region)
sim_background_events <- function(background_rate, time_window, spatial_region, covariate_columns = NULL) {

  spatial_area <- spatial_region |> sf::st_area() |> sum()
  t_length <- time_window[2] - time_window[1]

  background_rate <- as.numeric(background_rate)

  if (length(background_rate) == 1) {
    num_events <- stats::rpois(1, exp(background_rate) * spatial_area * t_length)


    if (num_events == 0) {
      background_events <- tibble::tibble(
                       x = NULL,
                       y = NULL,
                       t = NULL,
                       id = NULL,
                       parent = NULL,
                       gen = NULL,
                       family = NULL)
    }

    sample <- sf::st_sf(geometry = sf::st_sample(spatial_region, num_events, type = "random", exact = TRUE)) |>
      sf::st_as_sf()

    background_events <- sample |>
          dplyr::mutate(
                     x = sf::st_coordinates(sample)[,"X"],
                     y = sf::st_coordinates(sample)[,"Y"],
                     t = runif(num_events, time_window[1], time_window[2]),
                     id = 1:num_events,
                     parent = 0,
                     gen = 0,
                     family = 1:num_events) |>
      dplyr::relocate(.data$t, .after = .data$y) |>
      dplyr::arrange(.data$t)

  } else{
    X <- spatial_region |>
      sf::st_drop_geometry()
    X <- X[,covariate_columns] |>
      as.matrix()
    X <- cbind(1, X)

    num_events <- stats::rpois(nrow(X), exp(X %*% background_rate) * t_length * spatial_region$area)

    sample <- sf::st_sfc(crs = sf::st_crs(spatial_region))
    for (i in 1:nrow(spatial_region)) {
      region_sample <- sf::st_sfc(crs = sf::st_crs(spatial_region))
      while (length(region_sample) != num_events[i]) {
        region_sample <- c(region_sample, sf::st_sample(spatial_region[i,], size = 1, type = "random", exact = TRUE))
      }
      sample <- c(sample, region_sample)
    }

    background_events <- sf::st_sf(geometry = sample) |>
      sf::st_as_sf() |>
      dplyr::mutate(
        x = sf::st_coordinates(sample)[,"X"],
        y = sf::st_coordinates(sample)[,"Y"],
        t = runif(sum(num_events), time_window[1], time_window[2]),
        id = 1:sum(num_events),
        parent = 0,
        gen = 0,
        family = 1:sum(num_events)) |>
      dplyr::relocate(t, .after = .data$y) |>
      dplyr::arrange(t)
  }

  background_events
}



#' Generate a Hawkes Process
#'
#' @param params A named list of lists containing the values for the background rate, triggering rate, spatial parameters in a named list, and temporal parameters in a named list. See example for exact specification.
#' @param time_window A numeric vector of length 2 specifying the simulated time window.
#' @param spatial_region An sf object defining the spatial region for simulation.
#' @param covariate_columns A character vector of the names of the columns in spatial_region to be used as background covariates.
#' @param temporal_burnin Temporal burn-in for simulation. Defaults to 1/10 of the length of the time window if not used.
#' @param spatial_burnin Spatial burn-in for simulation. Defaults to the area of spatial_region^.25.
#' @param temporal_family A temporal triggering kernel function to generate data from. Defaults to "Exponential" if not used. Other options include "Power Law", "Uniform", and "Gamma".
#' @param spatial_family A spatial triggering kernel function to generate data from. Defaults to "Gaussian" if not used. Other options include "Uniform" and "Exponential"
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A hawkes object with a generated Hawkes process.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' (hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region, spatial_burnin = 1))
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
rHawkes <- function(params, time_window, spatial_region, covariate_columns = NULL,
                    temporal_burnin = (time_window[2] - time_window[1]) / (10), spatial_burnin = sum(sf::st_area(spatial_region) |> as.numeric())^.25,
                    temporal_family = "Exponential", spatial_family = "Gaussian") {
  # Create empty hawkes object and unpack to assign triggering sampler functions using the hawkes constructor
  hawkes <- hawkes(params = params, time_window = time_window, spatial_region = spatial_region,
         spatial_family = spatial_family, temporal_family = temporal_family)

  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  spatial_sampler    <- attrs$spatial_sampler
  temporal_sampler    <- attrs$temporal_sampler
  spatial_pdf  <- attrs$spatial_pdf
  temporal_pdf <- attrs$temporal_pdf
  spatial_cdf  <- attrs$spatial_cdf
  temporal_cdf <- attrs$temporal_cdf
  spatial_is_separable <- isTRUE(attrs$spatial_is_separable)



  # Check to see if covariates are included
  covariates <- !is.null(covariate_columns)

  # Set burnin regions
  if (spatial_burnin > 0) {

    spatial_region_burnin <- sf::st_buffer(spatial_region |> sf::st_union(), spatial_burnin) |> sf::st_as_sf()
    if (covariates) {
      spatial_region_burnin <- sf::st_difference(spatial_region_burnin, spatial_region |> sf::st_union()) |>
        sf::st_cast("POLYGON") |>
        sf::st_intersection(sf::st_make_grid(spatial_region_burnin, cellsize = spatial_burnin, square = TRUE))

      nearest_regions_ids <- sf::st_nearest_feature(sf::st_centroid(spatial_region_burnin), spatial_region)

      spatial_region_burnin <- spatial_region_burnin |>
        cbind({spatial_region[nearest_regions_ids, covariate_columns, drop = FALSE] |>
            sf::st_drop_geometry()
        }) |>
        sf::st_as_sf()

      spatial_region_burnin <- spatial_region |>
        dplyr::select(tidyselect::all_of(covariate_columns)) |>
        rbind(spatial_region_burnin |> dplyr::rename(geometry = x))

      spatial_region_burnin <- spatial_region_burnin |>
        dplyr::mutate(area = sf::st_area(spatial_region_burnin))
    }
  } else {
    spatial_region_burnin <- spatial_region
  }

  crs <- sf::st_crs(spatial_region_burnin)

  time_window_burnin <- time_window
  time_window_burnin[2] <- time_window[2] + temporal_burnin

  background_rate <- params$background_rate
  triggering_rate <- params$triggering_rate
  temporal_params <- params$temporal
  spatial_params <- params$spatial

  if (!is.list(background_rate)) {
    stop("background_rate must be named list stored within named params list.")
  }
  if (!((triggering_rate >= 0) && (triggering_rate < 1) && (is.numeric(triggering_rate)))) {
    stop("triggering_rate must be a numeric value between 0 and 1.")
  }

  # Generate background events
  data <- G <- sim_background_events(background_rate,
                                     time_window_burnin, spatial_region_burnin,
                                     covariate_columns = covariate_columns)


  # Specify generation l
  l <- 0

  while (TRUE) {
    O <- hawkes(params = params, time_window = time_window_burnin, spatial_region = spatial_region_burnin,
                spatial_family = spatial_family, temporal_family = temporal_family) |>
      dplyr::mutate(parent = numeric(), gen = numeric(), family = numeric(), .after = .data$t)
    sf::st_crs(O) <- crs

    l <- l+1


    N <- stats::rpois(nrow(G), triggering_rate)

    if(sum(N) == 0) {

      # Filter out buffer region and burning period
      data <- data |>
        dplyr::filter(t > temporal_burnin) |>
        dplyr::mutate(t = t - temporal_burnin) |>
        sf::st_filter(spatial_region) |>
        dplyr::arrange(t)

      if (covariates) {
        # Spatial join: assign covariates from polygon to each point
        data <- sf::st_join(data, spatial_region, join = sf::st_within) |>
          dplyr::select(tidyselect::all_of(c(names(data), covariate_columns, "area")))
      }

      data <- as_hawkes(data, time_window = time_window, spatial_region = spatial_region,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        covariate_columns = covariate_columns)

      return(data)
    }

    for (i in 1:nrow(G)) {
      if (N[i] > 0) {
        # Model specified self-exciting kernels. Add parameters here
        spatial_result <- do.call(spatial_sampler, c(list(n = N[i]), params$spatial))

        if (is.data.frame(spatial_result) && all(c("x", "y") %in% names(spatial_result))) {
          # Case A: kernel returns a data.frame with x and y
          x <- spatial_result$x + G$x[i]
          y <- spatial_result$y + G$y[i]
        } else if (is.numeric(spatial_result) && length(spatial_result) == N[i]) {
          # Case B: kernel returns a numeric vector (e.g., rnorm)
          x <- spatial_result + G$x[i]
          y <- do.call(spatial_sampler, c(list(n = N[i]), params$spatial)) + G$y[i]
        } else {
          stop("Invalid return from spatial_kernel: must be either vector or data.frame with x and y")
        }
        t <- do.call(temporal_sampler, c(list(n = N[i]), params$temporal)) + G$t[i]

        parent <- G$id[i]
        family <- G$family[i]

        # Create dataframe of events in the new generation O
        O_i <- data.frame(x = x, y = y, t = t, parent = parent, gen = l, family = family) |>
          sf::st_as_sf(coords = c("x", "y"), crs = crs)
        O <- rbind(O, O_i)
      }
    }
    O$id <- (1:nrow(O))+G$id[nrow(G)]

    # Filter events that lie outside of region or time bounds
    O <- O |>
      sf::st_as_sf(coords = c("x", "y")) |>
      sf::st_filter(spatial_region) |>
      dplyr::filter(t < time_window_burnin[2])

    O <- O |>
      dplyr::relocate(.data$t, .data$id, .before = parent) |>
      dplyr::mutate(
        x = sf::st_coordinates(O)[,1],
        y = sf::st_coordinates(O)[,2],
        .before = t)


    # Return data ordered by time if no events remain in new generation
    if (nrow(O) == 0) {

      # Filter out buffer region and burning period
      data <- data |>
        dplyr::filter(t > temporal_burnin) |>
        dplyr::mutate(t = t - temporal_burnin) |>
        sf::st_filter(spatial_region) |>
        dplyr::arrange(t)


      if (covariates) {
        # Spatial join: assign covariates from polygon to each point
        data <- sf::st_join(data, spatial_region, join = sf::st_within) |>
          dplyr::select(tidyselect::all_of(c(names(data), covariate_columns, "area")))
      }

      data <- as_hawkes(data, time_window = time_window, spatial_region = spatial_region,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        covariate_columns = covariate_columns)

      return(data)
    }

    # Add new generation to data and reset current generation to previous new one
    data <- rbind(data, O)
    G <- O
  }
}






