
#' Wrap data for block sampling
#'
#' @param hawkes A `hawkes` object
#' @param block_length_t A numeric value to set the length of blocks to wrap the process.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' extend_data_t_only(hawkes, 5)
#'
extend_data_t_only <- function(hawkes, block_length_t) {
  time_window <- attr(hawkes, "time_window")
  t_length <- time_window[2] - time_window[1]

  covariate_columns <- attr(hawkes, "covariate_columns")

  # create dataframes to extend points to hand block wrapping
  first_block <- transform(hawkes[hawkes$t <= block_length_t,], t = t + t_length)

  new_hawkes <- rbind(hawkes, first_block)

  time_window[2] <- time_window[2] + block_length_t

  attr(new_hawkes, "time_window") <- time_window
  new_hawkes
}


#' Resample data with moving blocls
#'
#' @param hawkes A `hawkes` object
#' @param block_length_t A numeric value to set the length of blocks to wrap the process.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#'
#' sample_blocks(hawkes, 25)
#'
sample_blocks <- function(hawkes, num_blocks) {

  .unpack_hawkes(hawkes)

  block_length <- (time_window[2] - time_window[1]) / num_blocks
  block_start <- runif(num_blocks, time_window[1], time_window[2])

  new_hawkes <- hawkes[0,]

  n_block <- 1

  # Sample a subinterval and rearrange to form a new block sample
  for (t_start in block_start) {
    block <- hawkes[(hawkes$t >= t_start & hawkes$t <= (t_start + block_length)),]
    block$t <- (block$t - t_start) + (n_block-1) * block_length

    n_block <- n_block + 1
    new_hawkes <- rbind(new_hawkes, block)
  }

  as_hawkes(new_hawkes, time_window, spatial_region, spatial_family, temporal_family, covariate_columns = covariate_columns)
}

#' Block Bootstrap for COnfidence Intervals of Hawkes MLEs
#'
#' @param hawkes A `hawkes` object
#' @param block_length_t A numeric value to set the length of blocks to wrap the process.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#'
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' block_bootstrap(hawkes, est, B = 2, num_blocks = 25, alpha = .05, parallel = FALSE, boundary = c(.5,3))
#'
#' future::plan(future::sequential)
block_bootstrap <- function(hawkes, est, B, num_blocks, alpha = .05, parallel = FALSE, max_iters = 500, boundary = NULL, t_burnin = 10, s_burnin = 0) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  # Wrap data to sample blocks beyond end
  t_length <- time_window[2] - time_window[1]
  block_length <- t_length / num_blocks
  extended_hawkes <- extend_data_t_only(hawkes, block_length)

  if (parallel) {
    n_failed <- 0

    boot_ests <- furrr::future_map_dfr(1:B, ~ tryCatch({
      sample <- sample_blocks(extended_hawkes, num_blocks)

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
      sample <- sample_blocks(extended_hawkes, num_blocks)
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
    dplyr::group_by(parameter_type, parameter) |>
    dplyr::summarise(boot_est_median = median(value),
                     lower = quantile(value, alpha/2),
                     upper = quantile(value, 1-(alpha/2)),
                     width = upper - lower,
                     sd = sd(value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
