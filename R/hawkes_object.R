# Functions to cosntruct a hawkes object
#

#' Constructor for a hawkes object
#'
#' @param data A dataframe containing the event locations in the columns x, y, and t. Defaults to NULL if not used.
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Defaults to NULL if not used.
#' @param time_window A numeric vector of length 2 specifying the time window.
#' @param spatial_region An sf object defining the spatial region and covariate regions.
#' @param covariate_columns A character vector containing the names of covariates columns to be used for the model.
#' @param spatial_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param temporal_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#'
#' @returns A hawkes object containing a tibble with the events.
#' @export
#'
hawkes <- function(data = NULL, params = NULL,
                   time_window = NULL, spatial_region = NULL,
                   spatial_family = NULL, temporal_family = NULL,
                   covariate_columns = NULL) {
  if (is.null(data)) {
    data <- data.frame(x = numeric(), y = numeric(), t = numeric()) |>
      sf::st_as_sf(coords = c("x", "y"), crs = NA) |>
      suppressWarnings() |>
      dplyr::mutate(x = numeric(), y = numeric(), .before = t)
  }


# Argument Checks ---------------------------------------------------------

  if (class(data)[1] != "sf" | !all(c("x", "y", "t") %in% names(data))) {
    stop("data must be a sf object and contain columns x, y, and t.")
  }

  if (is.null(spatial_family)) {
    stop("Provide spatial triggering family.")
  }

  if (is.null(temporal_family)) {
    stop("Provide temporal triggering family.")
  }

  if (!is.null(spatial_region) && !(class(spatial_region)[1] == "sf")) {
    print(spatial_region)
    stop("'spatial_region' must be a sf object defining the spatial region for the observed Hawkes process.")
  }

  if (!is.null(time_window) && !is.numeric(time_window) && !(length(time_window) == 2)) {
    stop("'time_window' must be a numeric vector defining the observed time window (e.g. c(0, 100))")
  }

  if (!is.null(params) && (!is.list(params) | is.null(params$background_rate) | is.null(params$triggering_rate) | is.null(params$spatial) | is.null(params$temporal))) {
    stop("Missing parameters. Make sure all components are provided in the params named list.")
  }

  # Assign Kernel Density and Sampler Functions -----------------------------------------
  if (!is.null(spatial_family) && class(spatial_family) == "character") {
    spatial_pdf <- switch (spatial_family,
                           "Gaussian" = stats::dnorm,
                           "Uniform" = stats::dunif,
                           "Exponential" = dexp_spatial,
                           stop("Spatial family is not supported.\nUse one of the provided spatial kernels (Guassian, Uniform, Exponential) or pass a list of kernel functions to the family arguments.")
    )
    spatial_cdf <- switch (spatial_family,
                           "Gaussian" = stats::pnorm,
                           "Uniform" = stats::punif,
                           "Exponential" = pexp_spatial,
                           stop("Spatial family is not supported.\nUse one of the provided spatial kernels (Guassian, Uniform, Exponential) or provide a density function to the spatial_family argument.")
    )
    spatial_sampler <- switch (spatial_family,
                               "Gaussian" = stats::rnorm,
                               "Uniform" = stats::runif,
                               "Exponential" = rexp_spatial,
                               stop("Spatial family is not supported.\nUse one of the provided spatial kernels (Guassian, Uniform, Exponential) or provide a density function to the spatial_family argument.")
    )
    spatial_is_separable <- if (spatial_family %in% c("Exponential")) {
      FALSE
    } else{
      TRUE
    }
  } else {
    message("Custom spatial kernel is being used\nEnsure the provided list includes a pdf, cdf, sampling function with the correct structure, and a variable specifying if the triggering intensity is separable. (Add help)")
    if (class(spatial_family) == "list") {
      spatial_pdf <- spatial_family$pdf
      spatial_cdf <- spatial_family$cdf
      spatial_sampler <- spatial_family$sampler
      spatial_is_separable <- spatial_family$is_separable
    }
  }

  # Assign the appropriate sampling method for the specified spatial kernel function.
  if (!is.null(temporal_family) && class(temporal_family) == "character") {
    temporal_pdf <- switch (temporal_family,
                            "Exponential" = stats::dexp,
                            "Gamma" = stats::dgamma,
                            "Uniform" = stats::dunif,
                            "Power Law" = dpower_law,
                            stop("Temporal family is not supported.\nUse one of the provided temporal kernels (Exponential, Gamma, Uniform, Power Law) or provide a density function to the temporal_family argument.")
    )
    temporal_cdf <- switch (temporal_family,
                            "Exponential" = stats::pexp,
                            "Gamma" = stats::pgamma,
                            "Unifrom" = stats::punif,
                            "Power Law" = ppower_law,
                            stop("Temporal family is not supported.\nUse one of the provided temporal kernels (Exponential, Gamma, Uniform, Power Law) or provide a density function to the spatial_family argument.")
    )
    temporal_sampler <- switch (temporal_family,
                                "Exponential" = stats::rexp,
                                "Gamma" = stats::rgamma,
                                "Uniform" = stats::runif,
                                "Power Law" = rpower_law,
                                stop("Temporal family is not supported.\nUse one of the provided temporal kernels (Exponential, Gamma, Uniform, Power Law) or provide a density function to the temporal_family argument.")
    )
  } else {
    # Add in support to provide all the custom functions in 1 list
    message("Custom temporal kernel is being used\nEnsure the provided list includes a pdf, cdf, and sampling fucntion with the correct structure. (Add help)")
    if (class(spatial_family) == "list") {
      temporal_pdf <- temporal_family$pdf
      temporal_cdf <- temporal_family$cdf
      temporal_sampler <- temporal_family$sampler
    }
  }

  if (!is.null(params) && !all(names(params$spatial) %in% formalArgs(spatial_pdf))) {
    stop(paste("Spatial parameter names are missing in spatial density function arguments."))
  }
  if (!is.null(params) && !all(names(params$temporal) %in% formalArgs(temporal_pdf))) {
    stop(paste("Spatial parameter names are missing in temporal density function arguments."))
  }

  if (!is.null(params) && !all(names(params$spatial) %in% formalArgs(spatial_sampler))) {
    stop(paste("Spatial parameter names are missing in spatial sampler function arguments."))
  }
  if (!is.null(params) && !all(names(params$temporal) %in% formalArgs(temporal_sampler))) {
    stop(paste("Spatial parameter names are missing in temporal sampler function arguments."))
  }


# Output object -----------------------------------------------------------

  structure(
    # data[,c(covariate_columns)],
    data |> dplyr::arrange(t),
    time_window = time_window,
    spatial_region = spatial_region,
    # params = params,
    covariate_columns = covariate_columns,
    spatial_family = spatial_family,
    temporal_family = temporal_family,
    spatial_sampler = spatial_sampler,
    spatial_pdf = spatial_pdf,
    spatial_cdf = spatial_cdf,
    temporal_pdf = temporal_pdf,
    temporal_cdf = temporal_cdf,
    temporal_sampler = temporal_sampler,
    spatial_is_separable = spatial_is_separable,
    class = c("hawkes", class(data))
  )
}

# Functions to generate a spatio-temporal Hawkes process
#


#' Convert an Object to Type Hawkes
#'
#' @param data A dataframe containing the event locations in the columns x, y, and t or an sf object containing the event locations in geometry and the event times in a column named t.
#' @param time_window A numeric vector of length 2 specifying the time window.
#' @param spatial_region An sf object defining the spatial region and covariate regions.
#' @param spatial_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param temporal_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param covariate_columns A character vector containing the names of covariates columns to be used for the model. Defaults to NULL if not used
#'
#' @returns A hawkes object.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   x = runif(100, 0, 10),
#'   y = runif(100, 0, 10),
#'   t = runif(100, 0, 50)
#' )
#'
#' # Convert to hawkes object
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#' hawkes_df <- as_hawkes(df, c(0,50), spatial_region, spatial_family = "Gaussian", temporal_family = "Exponential")
#' print(hawkes_df)
#'
as_hawkes <- function(data, time_window, spatial_region, spatial_family, temporal_family, covariate_columns = NULL) {
  if (class(data)[1] == "sf") {
    data <- data |>
      dplyr::mutate(
        x = sf::st_coordinates(data)[,1],
        y = sf::st_coordinates(data)[,2],
        .before = t
      )
  } else{
    data <- sf::st_as_sf(data, coords = c("x", "y"))
    data <- data |>
      dplyr::mutate(
        x = sf::st_coordinates(data)[,1],
        y = sf::st_coordinates(data)[,2],
        .before = t
      )
  }

  hawkes(data = data,
         time_window = time_window, spatial_region = spatial_region,
         spatial_family = spatial_family, temporal_family = temporal_family,
         covariate_columns = covariate_columns)
}


#' Print hawkes object
#'
#' @param x a hawkes object to be printed
#' @param n the number of events to print
#' @param ... arguments passed to or from other methods
#'
#' @returns The input hawkes object, invisibly.
#' @export
#'
print.hawkes <- function(x, n = 10, ...) {
  cat("<hawkes object>\n\n")
  cat("Number of events:", nrow(x), "\n\n")

  spatial_region <- attr(x, "spatial_region")
  time_window <- attr(x, "time_window")
  cat("Spatial Region:\n")
  print(spatial_region)

  cat("\nTime Window:\n")
  print(time_window)


    params <- attr(x, "params")
  if (!is.null(params)) {
    cat("Triggering Parameters:\n")

    cat("Background Rate (\u03B2):   \n")
    for (nm in names(params$background_rate)) {
      cat(sprintf("    %s: %s\n", nm, toString(round(params$background_rate[[nm]], 3))))
    }

    cat(sprintf(" Triggering Rate (\u03b8):   %s\n", params$triggering_rate))

  spatial_family <- attr(x, "spatial_family")
  if (is.null(spatial_family)) cat("\nSpatial kernel: not specified\n") else{
    cat("\nSpatial kernel: ",spatial_family)
  }

  temporal_family <- attr(x, "temporal_family")
  if (is.null(temporal_family)) cat("Temporal kernel: not specified\n") else{
    cat("\nTemporal kernel: ",temporal_family,"\n")
  }

    cat("\nSpatial Triggering Parameters:\n")
    for (nm in names(params$spatial)) {
      cat(sprintf("    %s: %s\n", nm, toString(params$spatial[[nm]])))
    }

    cat("Temporal Triggering Parameters:\n")
    for (nm in names(params$temporal)) {
      cat(sprintf("    %s: %s\n", nm, toString(params$temporal[[nm]])))
    }
  }


  total_rows <- nrow(x)
  shown_rows <- min(n, total_rows)

  cat("\nEvent data (first", shown_rows, "rows):\n")
  print.data.frame(x[seq_len(shown_rows), , drop = FALSE], ...)

  if (shown_rows < total_rows) {
    cat(sprintf("\033[90m# %d more rows\n", total_rows - shown_rows))
    cat("\033[90m# Use `print(n = ...)` to see more rows\033[39m\n")
  }

  invisible(x)
}



