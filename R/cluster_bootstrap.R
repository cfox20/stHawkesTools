#' Resample Data by EM-ALgorithm Clustering
#'
#' @param hawkes A `hawkes` object.
#' @param parent_mat An estimate of the parent matrix produced by `parent_est`.
#' @param boundary size of boundary to use for border correction. Defaults to NULL if unused.
#'
#' @returns A `hawkes` object.
#' @export
#'
#' @examples
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .25),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 1)
#' est <- hawkes_mle(hawkes, inits = params, boundary = 1)
#' parent_mat <- parent_est(hawkes, est)
#' sample_clusters(hawkes, parent_mat, boundary = c(.5,3))
#'
sample_clusters <- function(hawkes, parent_mat, boundary = NULL) {
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

  # Find which event is the parent of every event
  parent_id <- apply(parent_mat, 1, function(probs) {
    if (sum(probs) == 0) return(NA_integer_)  # no valid parent
    # This samples the parent of each event at each resample
    sample(seq_along(probs), size = 1, prob = probs)
    # To fix the parent struct use which.max
  })

  hawkes$parent_est <- parent_id

  # Recursively assign parents to estimate family origin
  repeat {
    grandparent_id <- hawkes$parent_est[hawkes$parent_est]
    if (all(hawkes$parent_est == grandparent_id)) break
    hawkes$parent_est <- grandparent_id
  }

  # Rename the parent_est column to cluster and adjust each event
  # relative to the origin of the cluster
  new_hawkes <- hawkes |>
    dplyr::rename(cluster = parent_est) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::mutate(x = .data$x - dplyr::first(.data$x),
                  y = .data$y - dplyr::first(.data$y),
                  t = .data$t - dplyr::first(.data$t)) |>
    dplyr::ungroup()

  # if (!is.null(cov_map)) {
  sampled_clusters <- hawkes |>
      dplyr::rename(cluster = parent_est) |>
      dplyr::group_by(.data$cluster) |>
      dplyr::slice(1) |>
      dplyr::select(.data$x, .data$y, .data$cluster) |>
      sf::st_as_sf(coords = c("x", "y"), crs = NA) |>
      sf::st_join(spatial_region) |>
      dplyr::select(.data$cluster, .data$geoid) |>
      sf::st_drop_geometry() |>
      dplyr::sample_n(size = rpois(1, dplyr::n()), replace = TRUE) |>
      dplyr::left_join(spatial_region |> dplyr::select(.data$geoid, .data$geometry), by = "geoid") |>
      sf::st_as_sf() |>
      dplyr::rowwise() |>
      dplyr::mutate(
        geometry = sf::st_sample(.data$geometry, size = 1, exact = TRUE) |> sf::st_cast("POINT")
      ) |>
      dplyr::ungroup()

  sampled_clusters <- sampled_clusters |>
      sf::st_drop_geometry() |>
      dplyr::select(.data$cluster) |>
      cbind(sf::st_coordinates(sampled_clusters)) |>
      dplyr::rename(x_shift = .data$X, y_shift = .data$Y) |>
      dplyr::mutate(t_shift = runif(length(.data$cluster), time_window[1], time_window[2]))

  new_hawkes <- new_hawkes |>
    dplyr::inner_join(sampled_clusters, by = "cluster", relationship = "many-to-many") |>
    dplyr::mutate(x = .data$x + .data$x_shift,
           y = .data$y + .data$y_shift,
           t = .data$t + .data$t_shift) |>
    dplyr::arrange(.data$t) |>
    dplyr::filter(t <= .env$time_window[2]) |>
    dplyr::select(-.data$x_shift, -.data$y_shift, -.data$t_shift) |>
    sf::st_drop_geometry() |>
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE) |>
    sf::st_set_crs(sf::st_crs(spatial_region)) |>
    dplyr::select(!tidyselect::all_of(c(.env$covariate_columns))) |>
    sf::st_intersection(spatial_region) |>
    suppressWarnings()


  as_hawkes(new_hawkes,
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
#' future::plan(future::multisession, workers = future::availableCores())
#'
#' cluster_bootstrap(hawkes, est, B = 2, alpha = .05, parallel = TRUE, boundary = c(.5,3))
#'
#' future::plan(future::sequential)
cluster_bootstrap <- function(hawkes, est, B, alpha = .05, parallel = FALSE, max_iters = 500, boundary = NULL) {
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


  # Estimate parent matrix using the original estimate
  parent_mat <- parent_est(hawkes, est)

  if (parallel) {
    n_failed <- 0

    boot_ests <- furrr::future_map_dfr(1:B, ~ tryCatch({
      # Sample estimated clusters to make a bootstrap sample
      sample <- sample_clusters(hawkes, parent_mat, boundary)

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
    n_failed <- 0

    boot_samples <- purrr::map(1:B, ~ {
      sample <- sample_clusters(hawkes, parent_mat, boundary)
      # sample <- sample_clusters2(hawkes, parent_mat, boundary)
    })

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
