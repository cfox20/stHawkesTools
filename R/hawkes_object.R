

# Functions to cosntruct a hawkes object
#

#' Constructor for a hawkes object
#'
#' @param data Data frame or `sf` object with columns `x`, `y`, and `t`. Defaults to `NULL`.
#' @param params Optional named list holding background, triggering, spatial, and temporal
#'   parameters.
#' @param time_window Numeric vector of length two defining the observation window.
#' @param spatial_region `sf` object describing the spatial domain and covariate regions.
#' @param covariate_columns Optional character vector naming background covariates.
#' @param mark_column Optional character string naming the mark column for multivariate processes.
#' @param spatial_family Spatial triggering kernel or list of custom kernel functions.
#' @param temporal_family Temporal triggering kernel or list of custom kernel functions.
#'
#' @returns A hawkes object containing a tibble of events.
#' @export
#'
hawkes <- function(data = NULL,
                   location_time_columns = c("x", "y", "t"),
                   background_formula = ~ 1,
                   mark_column = NULL,
                   time_window = NULL,
                   spatial_region = NULL,
                   spatial_family = NULL,
                   temporal_family = NULL,
                   parameters = NULL) {

  if (is.null(data)) {
    data <- data.frame(x = numeric(), y = numeric(), t = numeric()) |>
      sf::st_as_sf(coords = c("x", "y"), crs = NA) |>
      suppressWarnings() |>
      dplyr::mutate(x = numeric(), y = numeric(), .before = t)
  }

  if (!all(location_time_columns %in% names(data))) {
    stop("data must contain columns matching the location_time_columns argument.")
  }

  if (inherits(data, "sf")) {
    data <- data |>
      dplyr::mutate(
        x = sf::st_coordinates(data)[,1],
        y = sf::st_coordinates(data)[,2],
        .before = all_of(location_time_columns[3])
      ) |>
      dplyr::arrange(dplyr::all_of(location_time_columns[3]))
  } else{
    data <- sf::st_as_sf(data, coords = c(location_time_columns[1], location_time_columns[2]))
    data <- data |>
      dplyr::mutate(
        x = sf::st_coordinates(data)[,1],
        y = sf::st_coordinates(data)[,2],
        .before = all_of(location_time_columns[3])
      ) |>
      dplyr::arrange(dplyr::all_of(location_time_columns[3]))
  }


# Argument Checks ---------------------------------------------------------

  if (is.null(spatial_family)) {
    stop("Provide spatial self-excitation family.")
  }

  if (is.null(temporal_family)) {
    stop("Provide temporal self-excitation family.")
  }

  if (!is.null(spatial_region) && !(class(spatial_region)[1] == "sf")) {
    print(spatial_region)
    stop("'spatial_region' must be a sf object defining the spatial region for the observed Hawkes process.")
  }

  if (!is.null(time_window) && !is.numeric(time_window) && !(length(time_window) == 2)) {
    stop("'time_window' must be a numeric vector defining the observed time window (e.g. c(0, 100))")
  }

  if (!is.null(parameters) && (!is.list(parameters) | is.null(parameters$background_rate) | is.null(parameters$branching_ratio) | is.null(parameters$spatial) | is.null(parameters$temporal))) {
    stop("Missing parameters. Make sure all components are provided in the params named list.")
  }

  # Assign Kernel Density and Sampler Functions -----------------------------------------
  if (!is.null(spatial_family) && class(spatial_family)[1] == "character") {
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
    if (class(spatial_family)[1] == "list") {
      spatial_pdf <- spatial_family$pdf
      spatial_cdf <- spatial_family$cdf
      spatial_sampler <- spatial_family$sampler
      spatial_is_separable <- spatial_family$is_separable
    }
  }

  # Assign the appropriate sampling method for the specified spatial kernel function.
  if (!is.null(temporal_family) && class(temporal_family)[1] == "character") {
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
                            "Uniform" = stats::punif,
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
    if (class(spatial_family)[1] == "list") {
      temporal_pdf <- temporal_family$pdf
      temporal_cdf <- temporal_family$cdf
      temporal_sampler <- temporal_family$sampler
    }
  }

  if (!is.null(parameters) && !all(names(parameters$spatial) %in% methods::formalArgs(spatial_pdf))) {
    stop(paste("Spatial parameter names are missing in spatial density function arguments."))
  }
  if (!is.null(parameters) && !all(names(parameters$temporal) %in% methods::formalArgs(temporal_pdf))) {
    stop(paste("Spatial parameter names are missing in temporal density function arguments."))
  }

  if (!is.null(parameters) && !all(names(parameters$spatial) %in% methods::formalArgs(spatial_sampler))) {
    stop(paste("Spatial parameter names are missing in spatial sampler function arguments."))
  }
  if (!is.null(parameters) && !all(names(parameters$temporal) %in% methods::formalArgs(temporal_sampler))) {
    stop(paste("Spatial parameter names are missing in temporal sampler function arguments."))
  }

  # --- Validate parameters structure -----------------------------------------
  if (!is.null(parameters)) {

    # 1. Must be a list
    if (!is.list(parameters)) {
      stop("`parameters` must be a list.", call. = FALSE)
    }

    # 2. Required top-level names
    required_names <- c("background_rate", "branching_ratio", "spatial", "temporal")
    missing_names <- setdiff(required_names, names(parameters))
    if (length(missing_names) > 0) {
      stop("`parameters` is missing required components: ",
           paste(missing_names, collapse = ", "), call. = FALSE)
    }

    # 4. Match kernel argument names ------------------------------------------
    spatial_args  <- methods::formalArgs(spatial_pdf)
    temporal_args <- methods::formalArgs(temporal_pdf)

    bad_spatial  <- setdiff(names(parameters$spatial),  spatial_args)
    bad_temporal <- setdiff(names(parameters$temporal), temporal_args)

    if (length(bad_spatial) > 0) {
      stop("Unknown spatial parameter name(s): ",
           paste(bad_spatial, collapse = ", "), "\nValid names are: ",
           paste(spatial_args, collapse = ", "), call. = FALSE)
    }
    if (length(bad_temporal) > 0) {
      stop("Unknown temporal parameter name(s): ",
           paste(bad_temporal, collapse = ", "), "\nValid names are: ",
           paste(temporal_args, collapse = ", "), call. = FALSE)
    }

    # 5. Optional: check values are numeric and positive
    if (!all(sapply(parameters$spatial, is.numeric))) {
      stop("All spatial parameters must be numeric.", call. = FALSE)
    }
    if (!all(sapply(parameters$temporal, is.numeric))) {
      stop("All temporal parameters must be numeric.", call. = FALSE)
    }

    # (optional) branching ratio sanity check
    if (any(parameters$branching_ratio >= 1)) {
      warning("Triggering rate (branching ratio) ≥ 1 may lead to an unstable process.", call. = FALSE)
    }
  }


  if (nrow(data) > 1) {
  # Make matrix of covariate values for the observed data
    covariate_matrix <- .construct_background_covariate_matrix(background_formula = background_formula, data)
  } else{
    covariate_matrix <- NULL
  }

  if (nrow(data) > 1) {
  # Make matrix of covariate values for the observed data
    if (is.matrix(triggering_rate)){
      branching_matrix <- triggering_rate[data[[mark_column]], data[[mark_column]], drop = FALSE]
    } else {
      branching_matrix <- matrix(branching_ratio, nrow = nrow(branching_ratio), ncol = nrow(branching_ratio))
    }
    branching_matrix[upper.tri(branching_matrix, diag = TRUE)] <- 0
  } else{
    branching_matrix <- NULL
  }


# Output object -----------------------------------------------------------

  structure(
    data,
    covariate_matrix = covariate_matrix,
    branching_matrix = branching_matrix,
    background_formula = background_formula,
    mark_column = mark_column,
    time_window = time_window,
    spatial_region = spatial_region,
    parameters = parameters,
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

  background_formula <- attr(x, "background_formula")
  cat("Background Formula:")
  print(background_formula)
  # cat(
  #   "\nConditional Intensity λ(s, t):\n",
  #   sprintf("   λ(s, t) = exp{%s} + Σ g_t(t - tᵢ) · g_s(||s - sᵢ||)\n",
  #           deparse(attr(x, "background_formula"))),
  #   sep = ""
  # )

  spatial_region <- attr(x, "spatial_region")
  time_window <- attr(x, "time_window")
  # cat("Spatial Region:\n")
  # print(spatial_region)

  cat("\nTime Window:\n")
  print(time_window)


    params <- attr(x, "params")
  if (!is.null(params)) {
    cat("Triggering Parameters:\n")

    cat("Background Rate (\u03B2):   \n")
    for (nm in names(params$background_rate)) {
      cat(sprintf("    %s: %s\n", nm, toString(round(params$background_rate[[nm]], 3))))
    }

    cat(sprintf(" Branching Ratio (\u03b8):   %s\n", params$branching_ratio))

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



