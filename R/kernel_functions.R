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


#' Density of events with exponential(rate) distance in random direction from the origin
#'
#' @param x numeric vector of x coordinates
#' @param y numeric vector of y coordinates
#' @param rate rate parameter of the exponential distribution
#'
#' @returns numeric vector of densities at each (x, y) point
#' @export
#'
#' @examples
#' dexp_spatial(c(1, 0), c(0, 1), rate = 2)
dexp_spatial <- function(x, rate) {
  if (is.null(rate)) {
    stop("rate must be specified")
  }

  r <- sqrt(x[,1]^2 + x[,2]^2)
  density <- (rate / (2 * pi)) * exp(-rate * r)
  return(density)
}



#' CDF of events with exponential(rate) distance in random direction from the origin
#'
#' @param x numeric vector of x coordinates
#' @param y numeric vector of y coordinates
#' @param rate rate parameter of the exponential distribution
#'
#' @returns numeric vector of cumulative probabilities at each (x, y) point
#' @export
#'
#' @examples
#' pexp_spatial(c(1, 0), c(0, 1), rate = 2)
pexp_spatial <- function(q, rate) {
  if (is.null(rate)) {
    stop("rate must be specified")
  }

  r <- sqrt(q[,1]^2 + q[,2]^2)
  cdf <- 1 - exp(-rate * r)
  return(cdf)
}
