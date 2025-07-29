#' Resample data with moving blocls
#'
#' @param hawkes A `hawkes` object
#' @param length A numeric value to set the length of the subsample.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' sample_subregion(hawkes, 25)
#'
sample_subregion <- function(hawkes, length) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  sub_samp_t <- runif(1, region$t[1], region$t[2] - length)

  region$t <- c(sub_samp_t, sub_samp_t + length)

  sample_points <- ((hawkes$t >= region$t[1]) & (hawkes$t <= region$t[2]))
  new_hawkes <- hawkes[sample_points,]
  new_hawkes$t <- new_hawkes$t - sub_samp_t

  region$t <- c(0, length)

  attr(new_hawkes, "region") <- region

  if (exists("X", inherits = FALSE)) {
    attr(new_hawkes, "X") <- X[sample_points,]
  }

  new_hawkes
}

#' Block Bootstrap for COnfidence Intervals of Hawkes MLEs
#'
#' @param hawkes A `hawkes` object
#' @param length A numeric value to set the length of the subsample
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,150))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' subsample(hawkes, est, 5, 25, alpha = 0.05)
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1, X3 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, region, spatial_family = "Gaussian", temporal_family = "Exponential", cov_map = example_background_covariates)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#'
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' subsample(hawkes, est, B = 5, length = 20, alpha = .05, parallel = TRUE, seed = 123, boundary = c(.5,3))
#'
#' future::plan(future::sequential)
subsample <- function(hawkes, est, B, length, alpha, parallel = FALSE, seed = NULL, max_iters = 500, boundary = NULL, t_burnin = 10, s_burnin = 0) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if (!exists("cov_map", inherits = FALSE)) {
    cov_map <- NULL
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (parallel) {
    n_failed <- 0

    boot_ests <- furrr::future_map_dfr(1:B, ~ tryCatch({
      sample <- sample_subregion(hawkes, length)

      # Estimate the parameters on each bootstrapped sample and transform to long tibble
      boot_est <- hawkes_mle(sample, est, boundary = boundary, max_iters = max_iters)

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
      boot_est <- hawkes_mle(.x, est, boundary = boundary, max_iters = max_iters)
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
    dplyr::summarise(est = mean(value),
                     lower = quantile(value, alpha/2),
                     upper = quantile(value, 1-(alpha/2)),
                     width = upper - lower,
                     sd = sd(value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
