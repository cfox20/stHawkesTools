
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


  t_length <- time_window[2] - time_window[1]

  # create dataframes to extend points to hand block wrapping
  first_block <- transform(hawkes[hawkes$t <= block_length_t,], t = t + t_length)

  new_hawkes <- rbind(hawkes, first_block)

  time_window[2] <- time_window[2] + block_length_t

  as_hawkes(data = new_hawkes,
            time_window = time_window,
            spatial_region = spatial_region,
            spatial_family = spatial_family,
            temporal_family = temporal_family,
            covariate_columns = covariate_columns)
}


#' Resample data with moving blocls
#'
#' @param hawkes A `hawkes` object
#' @param num_blocks A numeric value to set the number of blocks to be used for resampling.
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

  as_hawkes(data = new_hawkes,
            time_window = time_window,
            spatial_region = spatial_region,
            spatial_family = spatial_family,
            temporal_family = temporal_family,
            covariate_columns = covariate_columns)
}

#' Block Bootstrap for Confidence Intervals of Hawkes MLEs
#'
#' @param hawkes A `hawkes` object
#' @param est A `hawkes_est` object containing the parameter estimates of `hawkes`
#' @param B number of bootstrap iterations
#' @param num_blocks number of blocs to sample for each bootstap iteration.
#' @param alpha type-1 error rate for constructing confidence intervals. Defaults to .05 if unused
#' @param parallel a logical specifying if parallel computation should be used. Parallel computation is implemented with the `furrr` package.
#' @param max_iters maximum number of iterations to use in maximum likelihood estimation. Defaults to 500 if unused.
#' @param boundary size of boundary to use for border correction. Defaults to NULL if unused
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
#' block_bootstrap(hawkes, est, B = 2, num_blocks = 25, alpha = .05, parallel = FALSE, boundary = 1)
#'
block_bootstrap <- function(hawkes, est, B, num_blocks, alpha = .05, parallel = FALSE, max_iters = 500, boundary = NULL) {
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
    dplyr::group_by(.data$parameter_type, .data$parameter) |>
    dplyr::summarise(boot_est_median = stats::median(.data$value),
                     lower = stats::quantile(.data$value, alpha/2),
                     upper = stats::quantile(.data$value, 1-(alpha/2)),
                     width = .data$upper - .data$lower,
                     sd = stats::sd(.data$value)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      percent_failed = n_failed / B
    )
}
