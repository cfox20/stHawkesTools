
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
sim_background_events <- function(background_rate, time_window, spatial_region, cov_map = NULL) {
  spatial_area <- spatial_region |> sf::st_area()
  t_length <- time_window[2] - time_window[1]

  background_rate <- as.numeric(background_rate)

  if (length(background_rate) == 1) {
    num_events <- stats::rpois(1, exp(background_rate) * spatial_area * t_length)

    sample <- sf::st_sample(spatial_region, num_events, exact = TRUE) |>
      sf::st_coordinates()

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
                     x = sample[,"X"],
                     y = sample[,"Y"],
                     t = runif(num_events, time_window[1], time_window[2]),
                     id = 1:num_events,
                     parent = 0,
                     gen = 0,
                     family = 1:num_events)

  } else{
    cov_names <- colnames(cov_map)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]

    X <- cov_map |>
      sf::st_drop_geometry()
    X <- X[,cov_names] |>
      as.matrix()
    X <- cbind(1, X)

    num_events <- stats::rpois(nrow(X), exp(X %*% background_rate) * region$t[2] * cov_map$area)

    sample <- sf::st_sfc(crs = sf::st_crs(cov_map))
    for (i in 1:nrow(cov_map)) {
      region_sample <- sf::st_sfc(crs = sf::st_crs(cov_map))
      while (length(region_sample) != num_events[i]) {
        region_sample <- c(region_sample, sf::st_sample(cov_map[i,], size = 1, type = "random", exact = TRUE))
      }
      sample <- c(sample, region_sample)
    }

    df <- sample |>
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

  df
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
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2))
#' rHawkes(params, region)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(min = -1, max = 1),temporal = list(shape = 3, rate = 1))
#' rHawkes(params, region, spatial_family = "Uniform", temporal_family = "Gamma")
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(rate = 5),temporal = list(rate = 2))
#' rHawkes(params, region, spatial_family = "Exponential", temporal_family = "Exponential")
#'
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1, X3 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' rHawkes(params, region, cov_map = example_background_covariates)
rHawkes <- function(params, time_window, spatial_region, temporal_burnin = 10, spatial_burnin = 2, cov_map = NULL, spatial_family = "Gaussian", temporal_family = "Exponential") {
  # Create empty hawkes object and unpack to assign triggering sampler functions using the hawkes constructor
  hawkes(params = params, time_window = time_window, spatial_region = spatial_region, spatial_family = spatial_family, temporal_family = temporal_family) |>
    .unpack_hawkes()

  # Set burnin regions
  spatial_region_burnin <- sf::st_buffer(spatial_region, spatial_burnin)

  temporal_window_burnin <- temporal_window
  temporal_window_burnin[2] <- temporal_window[2] + temporal_burnin

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
  covariates <- !is.null(cov_map)
  if (covariates) {
    cov_names <- colnames(X)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]
  }
  # Set X to null. Will be overwritten if there are covariates
  X <- NULL


  # Generate background events
  data <- G <- sim_background_events(background_rate,
                                     temporal_window_burnin, spatial_region_burnin,
                                     cov_map = cov_map)

  # Specify generation l
  l <- 0

  while (TRUE) {
    O <- hawkes(params = params, temporal_window = temporal_window_burnin, spatial_region = spatial_region_burnin,
                spatial_family = spatial_family, temporal_family = temporal_family)
    l <- l+1


    N <- stats::rpois(nrow(G), triggering_rate)

    if(sum(N) == 0) {

      data <- data[(data$t > t_burnin) & (data$x > x_min + s_burnin) & (data$x < x_max - s_burnin) & (data$y > y_min + s_burnin) & (data$y < y_max - s_burnin),]
      data$t <- data$t - t_burnin
      data$x <- data$x - s_burnin
      data$y <- data$y - s_burnin

      data <- data[order(data$t),]

      if (covariates) {
        # Ensure both are sf objects and in same CRS
        data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(cov_map))

        # Spatial join: assign covariates from polygon to each point
        x_values <- sf::st_join(data_sf, cov_map, join = sf::st_within) |>
          sf::st_drop_geometry() |>
          dplyr::select(id, cov_names, geoid, area)

        data <- dplyr::left_join(data, x_values, by = "id")

        X <- as.matrix(data[ ,cov_names])
        # data <- data[ ,!(colnames(data) %in% cov_names)]

      }

      data <- as_hawkes(data, region = region, params = params,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        cov_map = cov_map, X = X)
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

      data <- data[(data$t > t_burnin) & (data$x > x_min + s_burnin) & (data$x < x_max - s_burnin) & (data$y > y_min + s_burnin) & (data$y < y_max - s_burnin),]
      data$t <- data$t - t_burnin
      data$x <- data$x - s_burnin
      data$y <- data$y - s_burnin


      data <- data[order(data$t),]

      if (covariates) {
        # Ensure both are sf objects and in same CRS
        data_sf <- sf::st_as_sf(data, coords = c("x", "y"), crs = sf::st_crs(cov_map))

        # Spatial join: assign covariates from polygon to each point
        x_values <- sf::st_join(data_sf, cov_map, join = sf::st_within) |>
          sf::st_drop_geometry() |>
          dplyr::select(id, cov_names, geoid, area)

        data <- dplyr::left_join(data, x_values, by = "id")

        X <- as.matrix(data[ ,cov_names])
        # data <- data[ ,!(colnames(data) %in% cov_names)]

      }

      data <- as_hawkes(data, region = region, params = params,
                        spatial_family = spatial_family, temporal_family = temporal_family,
                        cov_map = cov_map, X = X)

      return(data)
    }

    # Add new generation to data and reset current generation to previous new one
    data <- rbind(data, O)
    G <- O
  }
}






