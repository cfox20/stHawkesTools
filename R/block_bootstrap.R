
#' Wrap data for block sampling
#'
#' @param hawkes A `hawkes` object
#' @param block_length_t A numeric value to set the length of blocks to wrap the process.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' extend_data_t_only(hawkes, 5)
#'
extend_data_t_only <- function(hawkes, block_length_t) {
  region <- attr(hawkes, "region")

  t_length <- region$t[2] - region$t[1]
  X <- attr(hawkes, "X")

  # create dataframes to extend points to hand block wrapping
  first_block <- transform(hawkes[hawkes$t <= block_length_t,], t = t + t_length)
  first_X <- X[hawkes$t <= block_length_t,]

  new_X <- rbind(X, first_X)
  new_hawkes <- rbind(hawkes, first_block)

  attr(new_hawkes, "X") <- new_X
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
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' sample_blocks(hawkes, 25)
#'
sample_blocks <- function(hawkes, num_blocks) {
  region <- attr(hawkes, "region")

  block_length <- (region$t[2] - region$t[1]) / num_blocks
  block_start <- runif(num_blocks, region$t[1], region$t[2])

  X <- attr(hawkes, "X")
  new_X <- X[0,]

  new_hawkes <- hawkes[0,]

  n_block <- 1

  # Sample a subinterval and rearrange to form a new block sample
  for (t_start in block_start) {
    block <- hawkes[(hawkes$t >= t_start & hawkes$t <= (t_start + block_length)),]
    X_block <- X[(hawkes$t >= t_start & hawkes$t <= (t_start + block_length)),]
    block$t <- (block$t - t_start) + (n_block-1) * block_length

    n_block <- n_block + 1
    new_hawkes <- rbind(new_hawkes, block)
    new_X <- rbind(new_X, X_block)
  }

  attr(new_hawkes, "X") <- new_X
  new_hawkes
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
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,150))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#' block_bootstrap(hawkes, est, 5, 25, alpha = 0.05)
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1, X3 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, region, spatial_family = "Gaussian", temporal_family = "Exponential", cov_map = example_background_covariates)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#'
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' block_bootstrap(hawkes, est, B = 1000, num_blocks = 25, alpha = .05, parallel = TRUE, boundary = c(.5,3))
#'
#' future::plan(future::sequential)
block_bootstrap <- function(hawkes, est, B, num_blocks, alpha, parallel = FALSE, max_iters = 500, boundary = NULL, t_burnin = 10, s_burnin = 0) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if (!exists("cov_map", inherits = FALSE)) {
    cov_map <- NULL
  }

  # Wrap data to sample blocks beyond end
  t_length <- region$t[2] - region$t[1]
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
