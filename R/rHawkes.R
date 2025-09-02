
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
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param background_rate A vector of coefficients for the background covariates.
#' @param cov_map An sf object containing background covariates.
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
sim_background_events <- function(background_rate, time_window, spatial_region,
                                  covariate_columns = covariate_columns, cov_map = NULL) {
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
      dplyr::relocate(t, .after = y) |>
      dplyr::arrange(t)

  } else{
    X <- cov_map |>
      sf::st_drop_geometry()
    X <- X[,covariate_columns] |>
      as.matrix()
    X <- cbind(1, X)

    num_events <- stats::rpois(nrow(X), exp(X %*% background_rate) * t_length * cov_map$area)

    sample <- sf::st_sfc(crs = sf::st_crs(cov_map))
    for (i in 1:nrow(cov_map)) {
      region_sample <- sf::st_sfc(crs = sf::st_crs(cov_map))
      while (length(region_sample) != num_events[i]) {
        region_sample <- c(region_sample, sf::st_sample(cov_map[i,], size = 1, type = "random", exact = TRUE))
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
      dplyr::relocate(t, .after = y) |>
      dplyr::arrange(t)
  }

  background_events
}



#' Generate a Hawkes Process
#'
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param t_burnin Temporal burn-in for simulation. Defaults to 10 if not used.
#' @param s_burnin Spatial burn-in for simulation. Note that this may not work as expected with a background covariate map. Defaults to 0 if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
#' @param spatial_family A spatial triggering kernel function to generate data from. Alternatively, a list can be provided to designate a custom kernel. The list must contain the objects named spatial_pdf, spatial_cdf, and spatial_sampler. They should follow the format of the dnorm, pnorm, and rnorm functions and the parameters must match the names. Defaults to NULL if not used.
#' @param temporal_family A spatial triggering kernel function to generate data from. Alternatively, a list can be provided to designate a custom kernel. The list must contain the objects named temporal_pdf, temporal_cdf, and temporal_sampler. They should follow the format of the dnorm, pnorm, and rnorm functions and the parameters must match the names. Defaults to NULL if not used.
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A hawkes object with a generated Hawkes process.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(min = -1, max = 1),temporal = list(shape = 3, rate = 1))
#' rHawkes(params, time_window = c(0,50), spatial_region = spatial_region, spatial_family = "Uniform", temporal_family = "Gamma")
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(rate = 5),temporal = list(rate = 2))
#' rHawkes(params, time_window = c(0,50), spatial_region = spatial_region, spatial_family = "Exponential", temporal_family = "Exponential")
#'
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' rHawkes(params, c(0,50), spatial_region, cov_map = example_background_covariates)
rHawkes <- function(params, time_window, spatial_region,
                    covariate_columns = NULL, cov_map = NULL,
                    temporal_burnin = 10, spatial_burnin = 2,
                    spatial_family = "Gaussian", temporal_family = "Exponential") {
  # Create empty hawkes object and unpack to assign triggering sampler functions using the hawkes constructor
  hawkes(params = params, time_window = time_window, spatial_region = spatial_region,
         spatial_family = spatial_family, temporal_family = temporal_family) |>
    .unpack_hawkes()

  # Set burnin regions
  spatial_region_burnin <- sf::st_buffer(spatial_region |> sf::st_union(), spatial_burnin) |> sf::st_as_sf()

  crs <- sf::st_crs(spatial_region)

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

  # Check to see if covariates are included
  covariates <- !is.null(covariate_columns)

  # Generate background events
  data <- G <- sim_background_events(background_rate,
                                     time_window_burnin, spatial_region_burnin,
                                     covariate_columns = covariate_columns, cov_map = cov_map)


  # Specify generation l
  l <- 0

  while (TRUE) {
    O <- hawkes(params = params, time_window = time_window_burnin, spatial_region = spatial_region_burnin,
                spatial_family = spatial_family, temporal_family = temporal_family) |>
      dplyr::mutate(parent = numeric(), gen = numeric(), family = numeric(), .after = t)
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
        data <- sf::st_join(data, cov_map, join = sf::st_within) |>
          dplyr::select(all_of(c(names(data), covariate_columns, "area")))
      }

      data <- as_hawkes(data, params = params, time_window = time_window, spatial_region = spatial_region,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        covariate_columns = covariate_columns, cov_map = cov_map)

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
      dplyr::relocate(t, id, .before = parent) |>
      dplyr::mutate(
        x = sf::st_coordinates(O)[,"X"],
        y = sf::st_coordinates(O)[,"Y"],
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
        data <- sf::st_join(data, cov_map, join = sf::st_within) |>
          dplyr::select(all_of(c(names(data), covariate_columns, "area")))
      }

      data <- as_hawkes(data, params = params, time_window = time_window, spatial_region = spatial_region,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        covariate_columns = covariate_columns, cov_map = cov_map)

      return(data)
    }

    # Add new generation to data and reset current generation to previous new one
    data <- rbind(data, O)
    G <- O
  }
}






