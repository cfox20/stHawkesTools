#' Compute Time Scaled Residuals
#'
#' @param hawkes a `hawkes` object
#' @param est a `hawkes_mle` object produced from `hawkes_mle()`.
#'
#' @returns A numeric vector of transformed waiting times.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,100))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 5))
#' residuals <- time_scaled_residuals(hawkes, est)
#'
#' plot(density(rexp(length(residuals)), from = 0), lwd = 2)
#' lines(density(residuals, from = 0), col = "red")
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(min = -1, max = 1),temporal = list(shape = 3, rate = 1))
#' hawkes <- rHawkes(params, region, spatial_family = "Uniform", temporal_family = "Gamma")
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 5))
#' residuals <- time_scaled_residuals(hawkes, est)
#'
#' plot(density(rexp(length(residuals)), from = 0), lwd = 2)
#' lines(density(residuals, from = 0), col = "red")
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1, X3 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, region, spatial_family = "Gaussian", temporal_family = "Exponential", cov_map = example_background_covariates)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' residuals <- time_scaled_residuals(hawkes, est)
#'
#' curve(dexp(x, rate = 1), from = 0, to = 10, lwd = 2, xlab = "Residual", ylab = "Density", main = "Residual vs Exponential(1)")
#' lines(density(residuals, from = 0), col = "red")
time_scaled_residuals <- function(hawkes, est) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  background_rate <- est$background_rate
  triggering_rate <- est$triggering_rate
  temporal_params <- est$temporal
  spatial_params <- est$spatial

  background_rate <- as.numeric(background_rate)

  if(!exists("X", inherits = FALSE)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, X)
  }

  x_min <- region$x[1]
  x_max <- region$x[2]
  y_min <- region$y[1]
  y_max <- region$y[2]
  t_min <- region$t[1]
  t_max <- region$t[2]

  x_diff <- outer(hawkes$x, hawkes$x, `-`)
  x_diff[upper.tri(x_diff, diag = FALSE)] <- 0
  y_diff <- outer(hawkes$y, hawkes$y, `-`)
  y_diff[upper.tri(y_diff, diag = FALSE)] <- 0

  t_i   <- hawkes$t
  t_im1 <- dplyr::lag(t_i, default = t_min)

  time_diff     <- outer(t_i, t_i, `-`)
  time_diff_lag <- outer(t_im1, t_i, `-`)

  # Then mask both:
  time_diff[!lower.tri(time_diff)] <- 0
  time_diff_lag[!lower.tri(time_diff_lag)] <- 0

  w <- hawkes$t - dplyr::lag(hawkes$t, default = t_min)

  # Store integral of background rate
  if (exists("cov_map")) {
    cov_cols <- colnames(X)[-1]

    cov_df <- cov_map |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(c(cov_cols, "area")))

    X_mat <- cbind(1, as.matrix(cov_df[cov_cols]))
    eta   <- drop(X_mat %*% background_rate)

    bg_space_integral <- sum(exp(eta) * cov_df$area)

    background_term <- w * bg_space_integral
  } else{
    background_term <- (x_max - x_min) * (y_max - y_min) * as.numeric(exp(X %*% background_rate) * w)
  }


  (x_max - x_min) * (y_max - y_min) * as.numeric(exp(X %*% background_rate) * w)

  # Integrate over entire spatial region and over waiting time
  triggering_term <-
    triggering_rate * {
      (do.call(spatial_cdf, c(list(q = x_max - hawkes$x), spatial_params)) -
        do.call(spatial_cdf, c(list(q = x_min - hawkes$x), spatial_params))) *
      (do.call(spatial_cdf, c(list(q = y_max - hawkes$y), spatial_params)) -
        do.call(spatial_cdf, c(list(q = y_min - hawkes$y), spatial_params))) *
      rowSums(do.call(temporal_cdf, c(list(q = time_diff), temporal_params)) -
        do.call(temporal_cdf, c(list(q = time_diff_lag), temporal_params)))
    }

  background_term + triggering_term
}

