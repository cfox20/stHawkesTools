
#' Parametric bootstrap for Hawkes MLE confidence intervals
#'
#' @param hawkes A `hawkes` object.
#' @param est A `hawkes_mle` object providing initial parameter estimates.
#' @param B Number of bootstrap iterations. At least 1000 is recommended.
#' @param alpha Significance level for confidence intervals. Defaults to 0.05.
#' @param parallel Logical flag for \pkg{furrr} based parallel computation. Configure a
#'   multisession plan before calling when `TRUE`.
#' @param max_iters Maximum EM iterations for each refit. Defaults to 500.
#' @param boundary Optional boundary width used for edge correction.
#' @param temporal_burnin Temporal burn-in used when simulating bootstrap samples.
#' @param spatial_burnin Spatial burn-in radius used when simulating bootstrap samples.
#'
#' @returns A tibble of bootstrap summary statistics.
#' @export
#'
#' @examples
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.25),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params,
#'   c(0, 50),
#'   example_background_covariates,
#'   covariate_columns = c("X1", "X2"),
#'   spatial_burnin = 1
#' )
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' # Demonstration with few bootstrap iterations (use at least 1000 in practice)
#' parametric_bootstrap(hawkes, est, B = 5, boundary = c(0.5, 3))
#'
#'
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.25),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params,
#'   c(0, 50),
#'   example_background_covariates,
#'   covariate_columns = c("X1", "X2"),
#'   spatial_burnin = 1
#' )
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' parametric_bootstrap(hawkes, est, B = 2, parallel = TRUE, boundary = c(0.5, 3))
#'
#' future::plan(future::sequential)
parametric_bootstrap <- function(hawkes, est, B, alpha = 0.05, parallel = FALSE, max_iters = 500,
                                 boundary = NULL, temporal_burnin = NULL, spatial_burnin = NULL) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

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


  if (is.null(spatial_burnin)) {
    spatial_burnin <- sum(sf::st_area(spatial_region))^.25
  }
  if (is.null(temporal_burnin)) {
    temporal_burnin <- (time_window[2] - time_window[1]) / (10)
  }

  if (parallel) {
  n_failed <- 0

  boot_ests <- furrr::future_map_dfr(1:B, ~ tryCatch({
    sample <- rHawkes(est$est, time_window, spatial_region,
                      covariate_columns = covariate_columns,
                      temporal_burnin = temporal_burnin,
                      spatial_burnin = spatial_burnin,
                      spatial_family = spatial_family,
                      temporal_family = temporal_family)

    # Estimate the parameters on each bootstrapped sample and transform to long tibble
    boot_est <- hawkes_mle(sample, inits = est$est, boundary = boundary, max_iters = max_iters)

    .hawkes_mle_to_dataframe(boot_est) |>
      dplyr::mutate(B = .x)

  }, error = function(e){
    warning(paste("Parameter estimation failed on bootstrap iteration:", .x,
                  "Due to error:", e))
    # If estimation failed return NA for the estimates and warn the user
    n_failed <- n_failed + 1
    .hawkes_mle_to_dataframe(est) |>
      dplyr::mutate(value = NA,
                    B = .x)

  }), .options = furrr::furrr_options(seed = TRUE), .progress = TRUE, packages = "stHawkesTools")

  } else{
    boot_samples <- purrr::map(1:B, ~ {
      sample <-  rHawkes(est$est, time_window, spatial_region,
                         covariate_columns = covariate_columns,
                         temporal_burnin = temporal_burnin, spatial_burnin = spatial_burnin,
                         spatial_family = spatial_family, temporal_family = temporal_family)
    })
    n_failed <- 0
    boot_ests <- purrr::imap_dfr(boot_samples, ~ tryCatch({
      # Estimate the parameters on each bootstrapped sample and transform to long tibble
      boot_est <- hawkes_mle(.x, inits = est$est, boundary = boundary, max_iters = max_iters)
      .hawkes_mle_to_dataframe(boot_est) |>
        dplyr::mutate(B = .y)
    }, error = function(e){
      warning(paste("Parameter estimation failed on bootstrap iteration:", .y,
                    "Due to error:", e))
      # If estimation failed return NA for the estimates and warn the user
      n_failed <<- n_failed + 1
      .hawkes_mle_to_dataframe(est) |>
        dplyr::mutate(value = NA,
                      B = .y)

    }))
  }

  if (nrow(boot_ests |> tidyr::drop_na()) == 0) {
    return({
      .hawkes_mle_to_dataframe(est) |>
        dplyr::mutate(lower = NA,
                      upper = NA,
                      width = NA)
    })
  }

  boot_ests |>
    tidyr::drop_na() |>
    dplyr::group_by(.data$parameter_type, .data$parameter) |>
    dplyr::summarise(
              boot_est_median = stats::median(.data$value),
              lower = stats::quantile(.data$value, .env$alpha/2),
              upper = stats::quantile(.data$value, 1-(.env$alpha/2)),
              width = .data$upper - .data$lower,
              sd = stats::sd(.data$value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
