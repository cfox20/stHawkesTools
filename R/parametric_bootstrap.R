
#' Parametric bootstrap for confidence intervals of Hawkes MLEs
#'
#' @param hawkes A `hawkes` object
#' @param est A `hawkes_mle` object to estimate the parameters for the `hawkes` object passed to the function.
#' @param B The number of bootstrap iterations. At least 1000 iterations are recommended, but computation time can be quite long.
#' @param alpha Specified alpha value. Defaults to 0.05 if not used.
#' @param parallel Logical to specify use of parallel computation via the \pkg{furrr} package. Setup multisession before calling function. See package for details.
#' @param seed A numeric to set seed for simulation. Defaults to NULL if not used. Note that a seed is recommended when using `parallel = TRUE`.
#' @param max_iters A numeric value for the maximum number of iteration in the EM-algorithm. Defaults to 500 if not used.
#' @param boundary A boundary region to correct for the bpundary bias. Defaults to NULL if not used.
#' @param t_burnin A numeric value to set the temporal burn-in region for simulation
#' @param s_burnin A numeric value to set the spatial burn-in region for simulation
#'
#' @returns A `hawkes_fit` object that is a list of lists contatiing the MLEs.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 5))
#' parametric_bootstrap(hawkes, est, B = 5, boundary = c(.5,3))
#'
#'
#' data("example_background_covariates")
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1, X3 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .1),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' hawkes <- rHawkes(params, region, spatial_family = "Gaussian", temporal_family = "Exponential", cov_map = example_background_covariates)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 3))
#'
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' parametric_bootstrap(hawkes, est, B = 5, parallel = TRUE, seed = 123, boundary = c(.5,3))
#'
#' future::plan(future::sequential)
parametric_bootstrap <- function(hawkes, est, B, alpha = 0.05, parallel = FALSE, seed = NULL, max_iters = 500, boundary = NULL, t_burnin = 10, s_burnin = 0) {
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
    sample <- rHawkes(est, region, t_burnin = t_burnin, s_burnin = s_burnin,
                      cov_map = cov_map, spatial_family = spatial_family, temporal_family = temporal_family)

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
      sample <- rHawkes(est, region, t_burnin = t_burnin, s_burnin = s_burnin,
                        cov_map = cov_map, spatial_family = spatial_family, temporal_family = temporal_family)
    })
    n_failed <- 0
    boot_ests <- purrr::imap_dfr(boot_samples, ~ tryCatch({
      # Estimate the parameters on each bootstrapped sample and transform to long tibble
      boot_est <- hawkes_mle(.x, est, boundary = boundary, max_iters = max_iters)
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
