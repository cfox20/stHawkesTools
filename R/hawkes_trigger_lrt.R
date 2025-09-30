
#' Likelihood Ratio Test to Detect Triggering
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Parameter values stored in a named list or a `hawkes_fit` object from `hawkes_mle()`.
#' @param alpha Type-1 error rate for the hypothesis test. Defaults to 0.05 if unused.
#'
#' @returns
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.01,
#'   spatial = list(mean = 0, sd = 0.5),
#'   temporal = list(rate = 1)
#' )
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 100),
#'   spatial_region = spatial_region,
#'   background_process = ~ 1
#' )
#' est <- hawkes_mle(hawkes, ~ 1, inits = params)
#' hawkes_trigger_lrt(hawkes, est)
#'
#'
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.1,
#'   spatial = list(mean = 0, sd = 0.25),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params = params,
#'   time_window = c(0, 50),
#'   spatial_region = example_background_covariates,
#'   background_process = ~ X1 + X2,
#'   spatial_burnin = 1
#' )
#' est <- hawkes_mle(hawkes, ~ X1 + X2, inits = params, boundary = 1)
#' hawkes_trigger_lrt(hawkes, est)
#'
hawkes_trigger_lrt <- function(hawkes, parameters, alpha = .05) {
  # Extract all hawkes object attributes
  attrs <- attributes(hawkes)

  # Assign all attributes to variables in the function environment
  time_window <- attrs$time_window
  spatial_region <- attrs$spatial_region
  covariate_columns    <- attrs$covariate_columns
  spatial_family    <- attrs$spatial_family
  temporal_family    <- attrs$temporal_family
  spatial_sampler    <- attrs$spatial_sampler
  temporal_sampler    <- attrs$temporal_sampler
  spatial_pdf  <- attrs$spatial_pdf
  temporal_pdf <- attrs$temporal_pdf
  spatial_cdf  <- attrs$spatial_cdf
  temporal_cdf <- attrs$temporal_cdf
  spatial_is_separable <- isTRUE(attrs$spatial_is_separable)


  if(is.null(covariate_columns)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry()) |>
      as.matrix()
  }

  if (inherits(parameters, "hawkes_fit")) {
    parameters <- parameters$est
  }

  background_rate <- parameters$background_rate |> as.numeric()

  null_background_rate <- nleqslv::nleqslv(background_rate, background_covariates_function, jac = NULL,
                   hawkes, parent_est_mat = diag(nrow = nrow(hawkes)), time_window,  spatial_region, covariate_columns, X)$x

  null_parameters <- parameters

  null_parameters$background_rate <- null_background_rate
  null_parameters$triggering_rate <- 0

  like_null <- log_likelihood(hawkes, null_parameters)
  like_alt <- log_likelihood(hawkes, parameters)

  test_statistic <- 2*(like_alt - like_null)

  p_val <- .5 * pchisq(test_statistic, df = 1, lower.tail = FALSE)

  tibble::tibble(
    test_statistic = test_statistic,
    alpha = alpha,
    p_value = p_val,
    reject_H0 = p_value < alpha
  )
}
