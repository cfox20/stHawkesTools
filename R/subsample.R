#' Resample data with moving blocls
#'
#' @param hawkes A `hawkes` object.
#' @param length Length of the temporal subsample.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.75,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(params, time_window = c(0, 50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' sample_subregion(hawkes, 25)
#'
sample_subregion <- function(hawkes, length) {
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


  sub_samp_t <- runif(1, time_window[1], time_window[2] - length)

  time_window <- c(sub_samp_t, sub_samp_t + length)

  sample_points <- ((hawkes$t >= time_window[1]) & (hawkes$t <= time_window[2]))
  new_hawkes <- hawkes[sample_points,]
  new_hawkes$t <- new_hawkes$t - sub_samp_t

  time_window <- c(0, length)

  new_hawkes |>
    as_hawkes(time_window = time_window,
              spatial_region = spatial_region,
              spatial_family = spatial_family,
              temporal_family = temporal_family)
}

#' Block bootstrap for Hawkes MLEs via subsampling
#'
#' @param hawkes A `hawkes` object.
#' @param est A `hawkes_est` object containing parameter estimates.
#' @param B Number of bootstrap iterations.
#' @param length Length of each temporal subsample.
#' @param alpha Type I error rate for confidence intervals. Defaults to 0.05.
#' @param parallel Logical flag for `furrr` based parallel computation.
#' @param max_iters Maximum EM iterations per fit. Defaults to 500.
#' @param boundary Optional boundary width used for edge correction.
#'
#' @returns A tibble of bootstrap estimates.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.75,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(params, time_window = c(0, 50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' subsample(hawkes, est, 5, 25, alpha = 0.05)
subsample <- function(hawkes, est, B, length, alpha, parallel = FALSE, max_iters = 500, boundary = NULL) {
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


  if(!is.null(covariate_columns)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry()) |>
      as.matrix()
  }

  if (parallel) {
    n_failed <- 0

    boot_ests <- furrr::future_map_dfr(1:B, ~ tryCatch({
      sample <- sample_subregion(hawkes, length)

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
      sample <- sample_subregion(hawkes, length)
    })
    n_failed <<- 0
    boot_ests <- purrr::imap_dfr(boot_samples, ~ tryCatch({
      # Estimate the parameters on each bootstrapped sample and transform to long tibble
      boot_est <- hawkes_mle(.x, inits = est$est, boundary = boundary, max_iters = max_iters)
      .hawkes_mle_to_dataframe(boot_est) |>
        dplyr::mutate(B = .y)
    }, error = function(e){
      warning(paste("Parameter estimation failed on bootstrap iteration:", .y,
                    "Due to error:", e))
      # If estimation failed return NA for the estimates and warn the user
      n_failed <- n_failed + 1
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
    dplyr::summarise(subsamp_est_median = stats::median(.data$value),
                     lower = stats::quantile(.data$value, alpha/2),
                     upper = stats::quantile(.data$value, 1-(alpha/2)),
                     width = .data$upper - .data$lower,
                     sd = stats::sd(.data$value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
