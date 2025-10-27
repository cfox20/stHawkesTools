#' Example Background Covariate sf Object
#'
#' A simple `sf` object with a rectangular region and a covariate column
#' `z` with higher values in the center of the region.
#'
#' @importFrom graphics plot
#'
#' @format An `sf` object with geometries and covariates named X1, X2, and X3
#' @source Bounding box and covariates derived from U.S. Census Bureau data accessed via `tidycensus`
#'
#' @examples
#' library("sf")
#'
#' data(example_background_covariates)
#' plot(example_background_covariates)
"example_background_covariates"




#' Example Simulated Hawkes Data
#'
#' A simple `hawkes` object simulated from the `example_background_covariates`
#' spatial region.
#'
#' @importFrom graphics plot
#'
#' @format A `hawkes` object with a covariate named `z`.
#'
#' @examples
#'
#' data(example_data)
#' plot_hawkes(example_data)
#'
#' inits <- list(
#'   background_rate = list(intercept = -5, z = .5),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' hawkes_mle(example_data, inits = inits, boundary = .5)
"example_data"

