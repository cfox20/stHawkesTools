
#' Create Rectangular SF Object
#'
#' @param xmin horizontal minimum of spatial region
#' @param xmax horizontal maximum of spatial region
#' @param ymin vertical minimum of spatial region
#' @param ymax vertical maximum of spatial region
#' @param covariates a matrix, data.frame, or tibble containing columns specifying covariate values for the grid.
#' The number of rows must equal `n_grid[1]*n_grid[2]`. The values fill in the region
#' left to right, then bottom to top.
#' @param n_grid numeric vector of length 2 specifying dimensions for the grid
#' of covariate regions. e.g. `n_grid = c(10,20)`
#' @param crs coordinate reference system. Defaults to NA if unused.
#'
#' @returns an sf object with a rectangular region
#' @export
#'
#' @examples
#' create_rectangular_sf(0,10,0,10)
create_rectangular_sf <- function(xmin, xmax, ymin, ymax, covariates = NULL, n_grid = c(1,1), crs = NA) {
  spatial_region <- sf::st_as_sfc(sf::st_bbox(c(xmin = xmin,
                                                ymin = ymin,
                                                xmax = xmax,
                                                ymax = ymax)), crs = crs) |>
    sf::st_as_sf() |>
    sf::st_make_grid(n = n_grid) |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = .data$x) |>
    dplyr::mutate(geoid = dplyr::row_number())

  spatial_region$area <- sf::st_area(spatial_region) |> as.numeric()

  # no covariates
  if (is.null(covariates)) return(spatial_region)

  spatial_region |>
    cbind(covariates)
}


#' Simulate background events
#'
#' @param background_rate Vector of coefficients for the background covariates.
#' @param time_window Numeric vector of length two giving the simulated time window.
#' @param spatial_region `sf` object defining the simulation region.
#' @param covariate_columns Optional character vector naming background covariates in
#'   `spatial_region`.
#' @param mark_column Optional character string identifying the column used for mark effects.
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A dataframe with generated background events.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#' time_window <- c(0, 50)
#' background_rate <- -4
#'
#' background_rate = list(intercept = -4, event_type = c(a = 2, b = 1))
#'
#' sim_background_events(background_rate, time_window, spatial_region, mark_column = "event_type")
sim_background_events <- function(background_rate, background_formula, time_window, spatial_region, mark_column = NULL) {

  spatial_area <- spatial_region |> sf::st_area() |> sum()
  t_length <- time_window[2] - time_window[1]

  mark_effects <- NULL
  background_rate_list <- background_rate
  if (!is.null(mark_column)) {
    mark_name <- mark_column
    mark_effects <- background_rate_list[[mark_name]]
    if (is.null(mark_effects)) {
      stop("`background_rate` must include mark-specific coefficients named `", mark_name, "`.")
    }
    if (is.null(names(mark_effects)) || anyNA(names(mark_effects))) {
      stop("Mark coefficients must be a named numeric vector specifying event types.")
    }
    background_rate_list[[mark_name]] <- NULL
    mark_effects <- as.numeric(mark_effects)
    names(mark_effects) <- names(background_rate[[mark_name]])
  }

  background_rate <- as.numeric(background_rate_list)

  if (length(background_rate) == 1) {
    if (is.null(mark_effects)) {
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
        return(background_events)
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
    } else {
      mark_levels <- names(mark_effects)
      event_counts <- stats::rpois(length(mark_effects), exp(background_rate + mark_effects) * spatial_area * t_length)
      total_events <- sum(event_counts)

      if (total_events == 0) {
        background_events <- tibble::tibble(
          x = NULL,
          y = NULL,
          t = NULL,
          id = NULL,
          parent = NULL,
          gen = NULL,
          family = NULL
        )
        background_events[[mark_column]] <- factor(levels = mark_levels)
        background_events[[mark_column]] <- background_events[[mark_column]][0]
        background_events <- background_events |>
          dplyr::relocate(tidyselect::all_of(mark_column), .after = .data$t)
        return(background_events)
      }

      sample <- sf::st_sfc(crs = sf::st_crs(spatial_region))
      event_types <- character(total_events)
      idx <- 1L
      for (j in seq_along(mark_effects)) {
        count <- event_counts[j]
        if (count == 0) next
        region_sample <- sf::st_sample(spatial_region, size = count, type = "random", exact = TRUE)
        sample <- c(sample, region_sample)
        event_types[idx:(idx + count - 1)] <- mark_levels[j]
        idx <- idx + count
      }

      sample <- sf::st_sf(geometry = sample) |>
        sf::st_as_sf()

      background_events <- sample |>
        dplyr::mutate(
          x = sf::st_coordinates(sample)[, "X"],
          y = sf::st_coordinates(sample)[, "Y"],
          t = runif(total_events, time_window[1], time_window[2]),
          id = seq_len(total_events),
          parent = 0,
          gen = 0,
          family = seq_len(total_events),
          "{mark_column}" := factor(event_types, levels = mark_levels)
        ) |>
        dplyr::relocate(.data$t, .after = .data$y) |>
        dplyr::relocate(tidyselect::all_of(mark_column), .after = .data$t) |>
        dplyr::arrange(.data$t)
    }

  } else{
    X <- spatial_region |>
      sf::st_drop_geometry()
    X <- X[,covariate_columns] |>
      as.matrix()
    X <- cbind(1, X)

    base_counts <- exp(X %*% background_rate) * t_length * spatial_region$area

    if (is.null(mark_effects)) {
      num_events <- stats::rpois(nrow(X), base_counts)

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
    } else {
      mark_levels <- names(mark_effects)
      num_events <- matrix(0, nrow = nrow(X), ncol = length(mark_effects))
      for (j in seq_along(mark_effects)) {
        num_events[, j] <- stats::rpois(nrow(X), base_counts * exp(mark_effects[j]))
      }

      total_events <- sum(num_events)

      if (total_events == 0) {
        background_events <- tibble::tibble(
          x = NULL,
          y = NULL,
          t = NULL,
          id = NULL,
          parent = NULL,
          gen = NULL,
          family = NULL
        )
        background_events[[mark_column]] <- factor(levels = mark_levels)
        background_events[[mark_column]] <- background_events[[mark_column]][0]
        background_events <- background_events |>
          dplyr::relocate(tidyselect::all_of(mark_column), .after = .data$t)
        return(background_events)
      }

      sample <- sf::st_sfc(crs = sf::st_crs(spatial_region))
      event_types <- character(total_events)
      idx <- 1L
      for (i in seq_len(nrow(spatial_region))) {
        for (j in seq_along(mark_effects)) {
          count <- num_events[i, j]
          if (count == 0) next
          region_sample <- sf::st_sfc(crs = sf::st_crs(spatial_region))
          while (length(region_sample) != count) {
            region_sample <- c(region_sample, sf::st_sample(spatial_region[i,], size = 1, type = "random", exact = TRUE))
          }
          sample <- c(sample, region_sample)
          event_types[idx:(idx + count - 1)] <- mark_levels[j]
          idx <- idx + count
        }
      }

      sample <- sf::st_sf(geometry = sample) |>
        sf::st_as_sf()

      background_events <- sample |>
        dplyr::mutate(
          x = sf::st_coordinates(sample)[, "X"],
          y = sf::st_coordinates(sample)[, "Y"],
          t = runif(total_events, time_window[1], time_window[2]),
          id = seq_len(total_events),
          parent = 0,
          gen = 0,
          family = seq_len(total_events),
          "{mark_column}" := factor(event_types, levels = mark_levels)
        ) |>
        dplyr::relocate(t, .after = .data$y) |>
        dplyr::relocate(tidyselect::all_of(mark_column), .after = .data$t) |>
        dplyr::arrange(t)
    }
  }

  background_events
}



#' Generate a Hawkes process
#'
#' @param hawkes Optional template `hawkes` object supplying process metadata. When
#'   omitted, `time_window`, `spatial_region`, and kernel families must be provided.
#' @param background_process One-sided formula specifying background covariates.
#'   When omitted and `hawkes` is provided, covariate information stored on the
#'   object is reused.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters. See the examples for the expected structure.
#' @param time_window Numeric vector of length two specifying the simulated window.
#'   Defaults to the window stored on `hawkes` when available.
#' @param spatial_region `sf` object defining the spatial region. Defaults to the
#'   region stored on `hawkes` when available.
#' @param temporal_burnin Temporal burn-in duration. Defaults to one tenth of the window
#'   length.
#' @param spatial_burnin Spatial burn-in radius. Defaults to `area(spatial_region)^0.25`.
#' @param temporal_family Temporal triggering kernel. Defaults to the family stored
#'   on `hawkes` when supplied, otherwise "Exponential".
#' @param spatial_family Spatial triggering kernel. Defaults to the family stored on
#'   `hawkes` when supplied, otherwise "Gaussian".
#'
#' @importFrom stats rnorm rpois rexp runif
#'
#' @returns A hawkes object with a generated Hawkes process.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' parameters <- list(
#'   background_rate = list(intercept = -4,
#'                          event_type = c(a = 1, b = .25, c = .5)),
#'   branching_ratio = matrix(c(.4, .15, .05,
#'                              .2, .05, .02,
#'                              .2, .05, .20),
#'                            nrow = 3,
#'                            dimnames = list(c("a", "b", "c"), c("a", "b", "c"))),
#'   spatial = list(mean = 0, sd = 0.1),
#'   temporal = list(rate = 2)
#' )
#' (hawkes <- rHawkes(
#'   parameters = parameters,
#'   time_window = c(0, 50),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1 + mark(event_type),
#'   spatial_burnin = 1
#' ))
#'
#' parameters <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1,
#'                          event_type_a = 1, event_type_b = .25),
#'   branching_ratio = matrix(c(.4, .15,
#'                              .2, .05),
#'                            nrow = 2,
#'                            dimnames = list(c("a", "b"), c("a", "b"))),
#'   spatial = list(mean = 0, sd = 0.25),
#'   temporal = list(rate = 2),
#'
#'   fixed = list(spatial = "mean")
#' )
#'
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   parameters = parameters,
#'   time_window = c(0, 50),
#'   spatial_region = example_background_covariates,
#'   background_formula = ~ X1 + X2 + event_type,
#'   spatial_burnin = 1
#' )
rHawkes <- function(hawkes = NULL, background_formula = ~ 1, mark_column = NULL,
                    parameters, time_window, spatial_region,
                    temporal_burnin = (time_window[2] - time_window[1]) / (10),
                    spatial_burnin = sum(sf::st_area(spatial_region) |> as.numeric())^.25,
                    temporal_family = "Exponential", spatial_family = "Gaussian") {
  if (!is.null(hawkes) && class(hawkes)[1] != "hawkes") {
    stop("hawkes must be a hawkes object or NULL.")
  }

  covariates <- background_formula != (~ 1)

  covariate_columns <- all.vars(background_formula)
  covariate_columns <- setdiff(covariate_columns, mark_column)

  if (missing(time_window) || is.null(time_window)) {
    if (!is.null(hawkes)) {
      time_window <- attr(hawkes, "time_window")
    } else {
      stop("time_window must be provided when hawkes is NULL.")
    }
  }

  if (missing(spatial_region) || is.null(spatial_region)) {
    if (!is.null(hawkes)) {
      spatial_region <- attr(hawkes, "spatial_region")
    } else {
      stop("spatial_region must be provided when hawkes is NULL.")
    }
  }


  # Create empty hawkes object and unpack to assign triggering sampler functions using the hawkes constructor
  hawkes <- hawkes(background_formula = background_formula,
                   mark_column = mark_column,
                   time_window = time_window,
                   spatial_region = spatial_region,
                   spatial_family = spatial_family,
                   temporal_family = temporal_family,
                   parameters = parameters)

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

  background_rate <- parameters$background_rate
  branching_ratio <- parameters$branching_ratio
  temporal_parameters <- parameters$temporal
  spatial_parameters <- parameters$spatial

  if (!is.list(background_rate)) {
    stop("background_rate must be named list stored within named parameters list.")
  }
  mark_levels <- NULL
  if (!is.null(mark_column) && mark_column %in% names(background_rate)) {
    mark_levels <- colnames(branching_ratio)
  }

  # Generate background events
  data <- G <- sim_background_events(background_rate, background_formula = background_formula,
                                     time_window_burnin, spatial_region_burnin,
                                     mark_column = mark_column)


  # Specify generation l
  l <- 0

  while (TRUE) {
    O <- hawkes(parameters = parameters, time_window = time_window_burnin, spatial_region = spatial_region_burnin,
                spatial_family = spatial_family, temporal_family = temporal_family, mark_column = mark_column) |>
      dplyr::mutate(parent = numeric(), gen = numeric(), family = numeric(), .after = .data$t)
    if (!is.null(mark_column)) {
      O <- O |>
        dplyr::mutate("{mark_column}" := factor(character(), levels = mark_levels), .after = .data$t)
    }
    sf::st_crs(O) <- crs

    l <- l+1

    if (is.null(mark_column)) {
      N <- stats::rpois(nrow(G), branching_ratio)
      total_children <- sum(N)
    } else {
      parent_marks <- as.character(G[[mark_column]])
      N <- vector("list", length = nrow(G))
      total_children <- 0L
      for (i in seq_len(nrow(G))) {
        child_rates <- branching_ratio[parent_marks[i], mark_levels]
        child_counts <- stats::rpois(length(child_rates), child_rates)
        N[[i]] <- child_counts
        total_children <- total_children + sum(child_counts)
      }
    }

    if(total_children == 0) {

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
                        covariate_columns = covariate_columns, mark_column = mark_column)

      return(data)
    }

    for (i in 1:nrow(G)) {
      if (is.null(mark_column)) {
        if (N[i] > 0) {
          spatial_result <- do.call(spatial_sampler, c(list(n = N[i]), parameters$spatial))

          if (is.data.frame(spatial_result) && all(c("x", "y") %in% names(spatial_result))) {
            x <- spatial_result$x + G$x[i]
            y <- spatial_result$y + G$y[i]
          } else if (is.numeric(spatial_result) && length(spatial_result) == N[i]) {
            x <- spatial_result + G$x[i]
            y <- do.call(spatial_sampler, c(list(n = N[i]), parameters$spatial)) + G$y[i]
          } else {
            stop("Invalid return from spatial_kernel: must be either vector or data.frame with x and y")
          }
          t <- do.call(temporal_sampler, c(list(n = N[i]), parameters$temporal)) + G$t[i]

          parent <- G$id[i]
          family <- G$family[i]

          O_i <- data.frame(x = x, y = y, t = t, parent = parent, gen = l, family = family) |>
            sf::st_as_sf(coords = c("x", "y"), crs = crs)
          O <- rbind(O, O_i)
        }
      } else {
        child_counts <- N[[i]]
        if (sum(child_counts) == 0) next

        for (j in seq_along(child_counts)) {
          count <- child_counts[j]
          if (count == 0) next

          spatial_result <- do.call(spatial_sampler, c(list(n = count), parameters$spatial))

          if (is.data.frame(spatial_result) && all(c("x", "y") %in% names(spatial_result))) {
            x <- spatial_result$x + G$x[i]
            y <- spatial_result$y + G$y[i]
          } else if (is.numeric(spatial_result) && length(spatial_result) == count) {
            x <- spatial_result + G$x[i]
            y <- do.call(spatial_sampler, c(list(n = count), parameters$spatial)) + G$y[i]
          } else {
            stop("Invalid return from spatial_kernel: must be either vector or data.frame with x and y")
          }
          t <- do.call(temporal_sampler, c(list(n = count), parameters$temporal)) + G$t[i]

          parent <- G$id[i]
          family <- G$family[i]
          marks <- rep(mark_levels[j], count)

          O_i <- data.frame(x = x, y = y, t = t, parent = parent, gen = l, family = family,
                            mark = marks)
          names(O_i)[names(O_i) == "mark"] <- mark_column
          O_i[[mark_column]] <- factor(O_i[[mark_column]], levels = mark_levels)
          O_i <- O_i |>
            sf::st_as_sf(coords = c("x", "y"), crs = crs)
          O <- rbind(O, O_i)
        }
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
                        covariate_columns = covariate_columns, mark_column = mark_column)

      return(data)
    }

    # Add new generation to data and reset current generation to previous new one
    data <- rbind(data, O)
    G <- O
  }
}






