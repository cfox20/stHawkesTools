# Functions to cosntruct a hawkes object
#

#' Constructor for a hawkes object
#'
#' @param data A dataframe containing the event locations in the columns x, y, and t. Defaults to NULL if not used.
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Defaults to NULL if not used.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)). Defaults to NULL if not used.
#' @param spatial_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param temporal_family A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
#' @param X A dataframe with nrow(data) rows containing the covariates observed at each event. Defaults to NULL if not used.
#'
#' @returns A hawkes object containing a tibble with the events.
#' @export
#'
hawkes <- function(data = NULL, params = NULL, region = NULL, spatial_family = NULL, temporal_family = NULL, cov_map = NULL, X = NULL) {
  if (is.null(data)) {
    data <- data.frame(x = numeric(), y = numeric(), t = numeric())
  }


# Argument Checks ---------------------------------------------------------

  if (!all(c("x", "y", "t") %in% names(data))) {
    stop("Data must contain columns 'x', 'y', and 't'.")
  }

  if (is.null(spatial_family)) {
    stop("Provide spatial triggering family.")
  }

  if (is.null(temporal_family)) {
    stop("Provide temporal triggering family.")
  }

  if (!is.null(region) && !is.list(region) && !is.numeric(region)) {
    stop("'region' should be a spatial window defined as a named list (e.g., list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax))).")
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
                            "Gamma" = stats::dgamma,
                            "Uniform" = stats::dunif,
                            "Exponential" = stats::dexp,
                            stop("Temporal family is not supported.\nUse one of the provided spatial kernels (Exponential, Gamma, Uniform) or provide a density function to the temporal_family argument.")
    )
    temporal_cdf <- switch (temporal_family,
                            "Gamma" = stats::pgamma,
                            "Unifrom" = stats::punif,
                            "Exponential" = stats::pexp,
                            stop("Spatial family is not supported.\nUse one of the provided spatial kernels (Guassian, Uniform, Exponential) or provide a density function to the spatial_family argument.")
    )
    temporal_sampler <- switch (temporal_family,
                                "Gamma" = stats::rgamma,
                                "Uniform" = stats::runif,
                                "Exponential" = stats::rexp,
                                stop("Temporal family is not supported.\nUse one of the provided spatial kernels (Exponential, Gamma, Uniform) or provide a density function to the temporal_family argument.")
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

  if (!is.null(X) & !is.null(cov_map)) {
    missing_covariates <- setdiff(colnames(X), colnames(cov_map))
    if (length(missing_covariates) > 0) {
      stop("The following covariates in X are missing from cov_map: ",
           paste(missing, collapse = ", "))
    }
  }

# Output object -----------------------------------------------------------

  structure(
    data[c("x", "y", "t", setdiff(names(data), c("x", "y", "t")))],
    params = params,
    region = region,
    spatial_family = spatial_family,
    temporal_family = temporal_family,
    cov_map = cov_map,
    X = X,
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
#' @param data A dataframe containing the event locations in the columns x, y, and t. Defaults to NULL if not used.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)). Defaults to NULL if not used.
#' @param spatial_family A spatial triggering kernel function to generate data from. Alternatively, a list can be provided to designate a custom kernel. The list must contain the objects named spatial_pdf, spatial_cdf, and spatial_sampler. They should follow the format of the dnorm, pnorm, and rnorm functions and the parameters must match the names. Defaults to NULL if not used.
#' @param temporal_family A spatial triggering kernel function to generate data from. Alternatively, a list can be provided to designate a custom kernel. The list must contain the objects named temporal_pdf, temporal_cdf, and temporal_sampler. They should follow the format of the dnorm, pnorm, and rnorm functions and the parameters must match the names. Defaults to NULL if not used.
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Defaults to NULL if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
#' @param X A dataframe with nrow(data) rows containing the covariates observed at each event. Defaults to NULL if not used.
#'
#' @returns A hawkes object.
#' @export
#'
#' @examples
#' set.seed(123)
#'
#' df <- data.frame(
#'   x = runif(100, 0, 10),
#'   y = runif(100, 0, 10),
#'   t = runif(100, 0, 50)
#' )
#'
#' # Convert to hawkes object
#' region <- list(x = c(0, 10), y = c(0, 10), t = c(0, 50))
#' hawkes_df <- as_hawkes(df, region, spatial_family = "Gaussian", temporal_family = "Exponential")
#' print(hawkes_df)
#'
as_hawkes <- function(data, region, spatial_family, temporal_family, params = NULL, cov_map = NULL, X = NULL) {
  stopifnot(is.data.frame(data))

  hawkes(data = data, params = params, region = region, spatial_family = spatial_family, temporal_family = temporal_family,
         cov_map = cov_map, X = X)
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

  region <- attr(x, "region")
  cat("Region:\n")
  cat(sprintf("  $x: [%g, %g]\n", region$x[1], region$x[2]))
  cat(sprintf("  $y: [%g, %g]\n", region$y[1], region$y[2]))
  cat(sprintf("  $t: [%g, %g]\n\n", region$t[1], region$t[2]))

  # cov_map <- attr(x, "cov_map")
  # if (!is.null(cov_map)) {
  #   print(cov_map)
  # }

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



