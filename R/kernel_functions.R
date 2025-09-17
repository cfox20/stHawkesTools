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
#' @param x a nx2 matrix of coordinates
#' @param rate rate parameter of the exponential distribution
#'
#' @returns numeric vector of densities at each (x, y) point
#' @export
#'
#' @examples
#' dexp_spatial(matrix(c(1, 0, 0, 1), ncol = 2), rate = 2)
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
#' @param q a nx2 matrix of coordinates
#' @param rate rate parameter of the exponential distribution
#'
#' @returns numeric vector of cumulative probabilities at each (x, y) point
#' @export
#'
#' @examples
#' pexp_spatial(matrix(c(1, 0, 0, 1), ncol = 2), rate = 2)
pexp_spatial <- function(q, rate) {
  if (is.null(rate)) {
    stop("rate must be specified")
  }

  r <- sqrt(q[,1]^2 + q[,2]^2)
  cdf <- 1 - exp(-rate * r)
  return(cdf)
}



#' Density function for power law (Lomax)
#'
#' @param x numeric vector
#' @param shape shape parameter for the Lomax distribution
#' @param scale scale parameter for the Lomax distribution
#'
#' @returns a numeric vector of densities
#' @export
#'
#' @examples
#' dpower_law(1, shape = 2, scale = 1)
#' dpower_law(matrix(c(1:48, NA), nrow = 7), shape = 2, scale = 1)
dpower_law <- function(x, shape = 2, scale = 1) {
  if (any(x[!is.na(x)] < 0)) stop("x must be nonnegative")
  if (any(shape <= 0)) stop("shape must be > 1")
  if (any(scale <= 0)) stop("scale must be > 0")

  dens <- (shape - 1) / scale * (1 + x / scale)^(-shape)
  return(dens)
}

#' CDF for power law (Lomax)
#'
#' @param q numeric vector
#' @param shape shape parameter for the Lomax distribution
#' @param scale scale parameter for the Lomax distribution
#' @param lower.tail logical; if TRUE (default), probabilities are $P[X< x]$ otherwise, $P[X>x]$.
#'
#' @returns a numeric vector of cumulative probabilities at each value of q
#' @export
#'
#' @examples
#' ppower_law(1, shape = 2, scale = 1)
ppower_law <- function(q, shape = 2, scale = 1, lower.tail = TRUE) {
  if (any(q[!is.na(q)] < 0)) stop("q must be nonnegative")
  if (any(shape <= 0)) stop("shape must be > 1")
  if (any(scale <= 0)) stop("scale must be > 0")

  cdf <- 1 - (1 + q / scale)^(-(shape - 1))
  if (!lower.tail) cdf <- 1 - cdf
  return(cdf)
}

#' Sample events with power law distance in time
#'
#' @param n number of observations
#' @param shape shape parameter for the Lomax distribution
#' @param scale scale parameter for the Lomax distribution
#'
#' @returns a numeric vector of values of the Lomax distribution
#' @export
#'
#' @examples
#' rpower_law(5, shape = 2, scale = 1)
rpower_law <- function(n, shape = 2, scale = 1) {
  if (is.null(shape) | is.null(scale)) {
    stop("params must include 'shape' and 'scale' for the power law kernel")
  }

  u <- runif(n)

  scale * ((1-u)^(-1/(shape-1))-1)
}
#






























