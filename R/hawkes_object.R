# Functions to cosntruct a hawkes object
#

#' Create hawkes object
#'
#' @param data A dataframe containing the event locations in the columns x, y, and t.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Defaults to NULL if not used.
#' @param spatial_kernel A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param temporal_kernel A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
#'
#' @returns A hawkes object containing a tibble with the events.
#' @export
#'
#' @examples
#' new_hawkes()
new_hawkes <- function(data = NULL, params = NULL, region = NULL, spatial_kernel = NULL, temporal_kernel = NULL, cov_map = NULL) {
  if (is.null(data)) {
    data <- data.frame(x = numeric(), y = numeric(), t = numeric())
  }

  structure(
    data[c("x", "y", "t", setdiff(names(data), c("x", "y", "t")))],
    params = params,
    region = region,
    kernel_functions = list(spatial_kernel = spatial_kernel, temporal_kernel = temporal_kernel),
    cov_map = cov_map,
    class = c("hawkes", class(data))
  )
}

#' Convert an Object to Type Hawkes
#'
#' @param data A dataframe containing the event locations in the columns x, y, and t.
#' @param region A list containing the spatial and temporal windows of the form list(x = c(xmin, xmax), y = c(ymin, ymax), t = c(tmin, tmax)).
#' @param params A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Defaults to NULL if not used.
#' @param spatial_kernel A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param temporal_kernel A spatial triggering kernel function to generate data from. Defaults to NULL if not used.
#' @param cov_map An sf object containing spatial polygons with associated covariate values named X1, X2, .... Defaults to NULL if not used.
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
#' hawkes_df <- as_hawkes(df, region)
#' print(hawkes_df)
as_hawkes <- function(data, region, params = NULL, spatial_kernel = NULL, temporal_kernel = NULL, cov_map = NULL) {
  stopifnot(is.data.frame(data))

  if (!all(c("x", "y", "t") %in% names(data))) {
    stop("Data must contain columns 'x', 'y', and 't'.")
  }

  if (!is.list(region) && !is.numeric(region)) {
    stop("'region' should be a bounding box or spatial window (e.g., named list).")
  }

  new_hawkes(data, params, region, spatial_kernel, temporal_kernel, cov_map)
}


#' Title
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

  cov_map <- attr(x, "cov_map")
  if (!is.null(cov_map)) {
    print(cov_map)
  }

    params <- attr(x, "params")
  if (!is.null(params)) {
    cat("Triggering Parameters:\n")
    cat(sprintf(" Background Rate (\u03B2):   %s\n", params$background_rate))
    cat(sprintf(" Triggering Rate (\u03b8):   %s\n", params$triggering_rate))

    cat(" Spatial Triggering Parameters:\n")
    for (nm in names(params$spatial)) {
      cat(sprintf("    %s: %s\n", nm, toString(params$spatial[[nm]])))
    }

    cat(" Temporal Triggering Parameters:\n")
    for (nm in names(params$temporal)) {
      cat(sprintf("    %s: %s\n", nm, toString(params$temporal[[nm]])))
    }
  }

  kernels <- attr(x, "kernel_functions")
  if (!is.null(kernels$spatial_kernel)) cat("\nSpatial kernel: defined\n")
  if (!is.null(kernels$temporal_kernel)) cat("Temporal kernel: defined\n")

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
