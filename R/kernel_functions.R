#' Sample events with exponential(rate) distance in random direction from the origin
#'
#' @param n number of observations
#' @param rate vector of rates
#'
#' @returns a dataframe of points located exp(rate) away from the origin
#' @export
#'
#' @examples
#' rexp_spatial(5, 2)
rexp_spatial <- function(n, rate) {
  if (is.null(rate)) {
    stop("params must include 'rate' for the exponential kernel")
  }

  r <- rexp(n, rate = rate)
  theta <- runif(n, 0, 2 * pi)

  data.frame(
    x = r * cos(theta),
    y = r * sin(theta)
  )
}
