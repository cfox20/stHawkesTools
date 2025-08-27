#' Example Background Covariate sf Object
#'
#' A simple `sf` object with a rectangular spatial region surrounding Waco, TX.
#' Three scaled covariate columns are included from U.S. Census Bureau data,
#' useful for simulating background intensity in Hawkes processes.
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
