#' Compute Time Scaled Residuals
#'
#' @param hawkes a `hawkes` object
#' @param est a `hawkes_mle` object produced from `hawkes_mle()`.
#'
#' @returns A numeric vector of transformed waiting times.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,100), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' residuals <- time_scaled_residuals(hawkes, est)
#'
#' curve(dexp(x, rate = 1), from = 0, to = 10, lwd = 2, xlab = "Residual", ylab = "Density", main = "Residual vs Exponential(1)")
#' lines(density(residuals, from = 0), col = "red")
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' curve(dexp(x, rate = 1), from = 0, to = 10, lwd = 2, xlab = "Residual", ylab = "Density", main = "Residual vs Exponential(1)")
#' lines(density(residuals, from = 0), col = "red")
time_scaled_residuals <- function(hawkes, est) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if (class(est)[1] == "hawkes_fit") {
    est <- est$est
  }

  background_rate <- est$background_rate
  triggering_rate <- est$triggering_rate
  temporal_params <- est$temporal
  spatial_params <- est$spatial

  background_rate <- as.numeric(background_rate)

  if(!exists("covariate_columns", inherits = FALSE)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry()) |>
      as.matrix()
  }

  spatial_area <- sf::st_area(spatial_region) |> sum()

  x_diff <- outer(hawkes$x, hawkes$x, `-`)
  x_diff[upper.tri(x_diff, diag = FALSE)] <- 0
  y_diff <- outer(hawkes$y, hawkes$y, `-`)
  y_diff[upper.tri(y_diff, diag = FALSE)] <- 0

  t_i   <- hawkes$t
  t_im1 <- dplyr::lag(t_i, default = time_window[1])

  time_diff     <- outer(t_i, t_i, `-`)
  time_diff_lag <- outer(t_im1, t_i, `-`)

  # Then mask both:
  time_diff[!lower.tri(time_diff)] <- 0
  time_diff_lag[!lower.tri(time_diff_lag)] <- 0

  w <- hawkes$t - dplyr::lag(hawkes$t, default = time_window[1])

  # Store integral of background rate
  if (exists("covariate_columns", inherits = FALSE)) {
    cov_df <- spatial_region |>
      sf::st_drop_geometry() |>
      dplyr::select(dplyr::all_of(c(covariate_columns, "area")))

    X_mat <- cbind(1, as.matrix(cov_df[covariate_columns]))
    eta   <- drop(X_mat %*% background_rate)

    bg_space_integral <- sum(exp(eta) * cov_df$area)

    background_term <- w * bg_space_integral
  } else{
    background_term <- spatial_area * as.numeric(exp(X %*% background_rate) * w)
  }


  spatial_area * as.numeric(exp(X %*% background_rate) * w)

  # Integrate over entire spatial region and over waiting time
  triggering_term <-
    triggering_rate * {
      # (do.call(spatial_cdf, c(list(q = x_max - hawkes$x), spatial_params)) -
      #   do.call(spatial_cdf, c(list(q = x_min - hawkes$x), spatial_params))) *
      # (do.call(spatial_cdf, c(list(q = y_max - hawkes$y), spatial_params)) -
      #   do.call(spatial_cdf, c(list(q = y_min - hawkes$y), spatial_params))) *
      rowSums(do.call(temporal_cdf, c(list(q = time_diff), temporal_params)) -
        do.call(temporal_cdf, c(list(q = time_diff_lag), temporal_params)))
    }

  as.numeric(background_term) + triggering_term
}

