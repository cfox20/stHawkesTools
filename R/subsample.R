#' Resample data with moving blocls
#'
#' @param hawkes A `hawkes` object
#' @param length A numeric value to set the length of the subsample.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' sample_subregion(hawkes, 25)
#'
sample_subregion <- function(hawkes, length) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

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

#' Block Bootstrap for COnfidence Intervals of Hawkes MLEs
#'
#' @param hawkes A `hawkes` object
#' @param est A `hawkes_est` object containing the parameter estimates of `hawkes`
#' @param B number of bootstrap iterations
#' @param length A numeric value to set the length of the subsample
#' @param alpha type-1 error rate for constructing confidence intervals. Defaults to .05 if unused
#' @param parallel a logical specifying if parallel computation should be used. Parallel computation is implemented with the `furrr` package.
#' @param max_iters maximum number of iterations to use in maximum likelihood estimation. Defaults to 500 if unused.
#' @param boundary size of boundary to use for border correction. Defaults to NULL if unused
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' subsample(hawkes, est, 5, 25, alpha = 0.05)
subsample <- function(hawkes, est, B, length, alpha, parallel = FALSE, max_iters = 500, boundary = NULL) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if(!exists("covariate_columns", inherits = FALSE)){
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

    }), .options = furrr::furrr_options(seed = !is.null(seed)), .progress = TRUE, packages = "stHawkesTools")

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
    dplyr::group_by(parameter_type, parameter) |>
    dplyr::summarise(subsamp_est_median = median(value),
                     lower = quantile(value, alpha/2),
                     upper = quantile(value, 1-(alpha/2)),
                     width = upper - lower,
                     sd = sd(value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
