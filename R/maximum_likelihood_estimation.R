
#' Parent Matrix Estimation
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters. Values usually represent current estimates.
#'
#' @returns A numeric matrix
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
#'
#' (parent_est_mat <- parent_est(hawkes, params))
#'
#' params <- list(
#'   background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.75),
#'   temporal = list(rate = 2),
#'   fixed = list(spatial = "mean")
#' )
#' data("example_background_covariates")
#' hawkes <- rHawkes(
#'   params,
#'   c(0, 50),
#'   example_background_covariates,
#'   covariate_columns = c("X1", "X2"),
#'   spatial_burnin = 0
#' )
#'
#' (parent_est_mat <- parent_est(hawkes, params))
parent_est <- function(hawkes, parameters) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

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


  if(is.null(covariate_columns)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry()) |>
      as.matrix()
  }

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate <- as.numeric(background_rate)

  # Store the time and space differences between each point in diagonal matrices
  x_diff <- outer(hawkes$x, hawkes$x, `-`)
  x_diff[upper.tri(x_diff, diag = TRUE)] <- NA
  y_diff <- outer(hawkes$y, hawkes$y, `-`)
  y_diff[upper.tri(y_diff, diag = TRUE)] <- NA
  time_diff <- outer(hawkes$t, hawkes$t, `-`)
  time_diff[upper.tri(time_diff, diag = TRUE)] <- NA

  # Compute a matrix of the triggering intensities
  if (!spatial_is_separable) {
    s_diff <- cbind(x_diff, y_diff)
    g_mat <- {triggering_rate *
        do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
        do.call(spatial_pdf, c(list(x = s_diff), parameters$spatial))
    }
  } else {
    g_mat <- {triggering_rate *
        do.call(temporal_pdf, c(list(x = time_diff), parameters$temporal)) *
        do.call(spatial_pdf, c(list(x = x_diff), parameters$spatial)) *
        do.call(spatial_pdf, c(list(x = y_diff), parameters$spatial))
    }
  }
  g_mat[upper.tri(g_mat, diag = TRUE)] <- 0

  # Store the values of the complete likelihood at each point
  lambda_i <- exp(as.numeric(X %*% background_rate)) + rowSums(g_mat)
  lambda_mat <- matrix(lambda_i, nrow = length(lambda_i), ncol = length(lambda_i), byrow = FALSE)

  # Estimate the parent matrix
  parent_mat <- g_mat / lambda_mat
  diag(parent_mat) <-exp(as.numeric(X %*% background_rate)) / diag(lambda_mat)
  parent_mat[upper.tri(parent_mat)] <- 0

  colnames(parent_mat) <- 1:nrow(hawkes)
  parent_mat
}

#' Compute parameter estimate updates
#'
#' @param hawkes A `hawkes` object.
#' @param parameters Named list containing background, triggering, spatial, and temporal
#'   parameters. To fix kernel parameters, include a `fixed` list such as
#'   `fixed$spatial = c("mean")`.
#' @param parent_est_mat Matrix returned by `parent_est()`.
#' @param boundary Optional boundary width used to correct edge bias.
#' @param fixed_spatial Optional character vector naming spatial parameters to hold fixed.
#' @param fixed_temporal Optional character vector naming temporal parameters to hold
#'   fixed.
#'
#' @returns A named list mirroring `parameters` with updated estimates.
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
#' parent_est_mat <- parent_est(hawkes, params)
#' est_params(hawkes, params, parent_est_mat)
#'
est_params <- function(hawkes, parameters, parent_est_mat, boundary = NULL, fixed_spatial = NULL, fixed_temporal = NULL) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")
  if(class(parent_est_mat)[[1]] != "matrix") stop("parent_est_mat must be a matrix")

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


  if(is.null(covariate_columns)){
    X <- matrix(rep(1,nrow(hawkes)), ncol = 1)
  } else {
    X <- cbind(1, hawkes[,covariate_columns, .drop = FALSE] |> sf::st_drop_geometry()) |>
      as.matrix()
  }

  spatial_area <- spatial_region |>
    sf::st_area() |>
    sum()

  background_rate <- parameters$background_rate
  triggering_rate <- parameters$triggering_rate
  temporal_params <- parameters$temporal
  spatial_params <- parameters$spatial

  background_rate_names <- names(background_rate)
  background_rate <- as.numeric(background_rate)

  # Create matrices with the differences in location and times between each point
  x_diff <- outer(hawkes$x, hawkes$x, `-`)
  x_diff[upper.tri(x_diff)] <- 0
  y_diff <- outer(hawkes$y, hawkes$y, `-`)
  y_diff[upper.tri(y_diff)] <- 0

  time_diff <- outer(hawkes$t, hawkes$t, `-`)
  time_diff[upper.tri(time_diff)] <- 0


  # Estimation if a boundary is being used
  if (!is.null(boundary)) {
    spatial_region_boundary <- spatial_region |>
      sf::st_union() |>
      sf::st_buffer(-boundary)

    spatial_region <- spatial_region |>
      sf::st_intersection(spatial_region_boundary) |>
      suppressWarnings()


    # Find the index of values that fall within the boundary region S0
    outside_S <- lengths(sf::st_within(hawkes, spatial_region)) <= 0
    # Set the rows of the events in S0 to 0 in the parent matrix
    parent_est_mat[outside_S,] <- 0

    hawkes <- hawkes[!outside_S,]
  }

  if (length(background_rate) > 1) {
    background_rate_update <- nleqslv::nleqslv(background_rate, background_covariates_function, jac = NULL, hawkes, parent_est_mat, time_window,  spatial_region, covariate_columns, X)$x
  } else {
    background_rate_update <- log(sum(diag(parent_est_mat)) / ((time_window[2] - time_window[1]) * spatial_area))
  }
  background_rate_update <- as.list(background_rate_update)
  names(background_rate_update) <- background_rate_names

  diag(parent_est_mat) <- 0

  triggering_rate_update <- sum(parent_est_mat) / sum(do.call(temporal_cdf, c(list(q = time_window[2] - hawkes$t), temporal_params)))

  # Hardcode analytic solutions if available for the MLE

  # switch(spatial_family,
  #   "Gaussian" = .spatial_update_gaussian(),
  #   "Uniform" = "no",
  #   stop("Spatial parameter update not found")
  # )
  #
  # switch(temporal_family,
  #   "Exponential" = .temporal_update_gaussian(hawkes),
  #   "Uniform" = "no",
  #   stop("Spatial parameter update not found")
  # )


  if (spatial_family == "Gaussian") {
    # triggering_rate_update <- sum(parent_est_mat) / (sum(1 - exp(-temporal_params$rate * (t_max - hawkes$t))))

    spatial_param_updates <- list(mean = 0, sd = sqrt(sum(parent_est_mat * (x_diff^2 + y_diff^2)) / (2 * sum(parent_est_mat))))

  } else {
    # Numerically solve for MLE for general kernels

    spatial_param_updates <- .optim_named(.spatial_parameter_likelihood, parameters, "spatial", fixed = fixed_spatial, hawkes = hawkes,
                                          parent_est_mat = parent_est_mat, x_diff = x_diff, y_diff = y_diff)
  }

  if (temporal_family == "Exponential") {
    temporal_param_updates <- list(rate = sum(parent_est_mat) / (sum(parent_est_mat * time_diff) + triggering_rate_update * sum((time_window[2] - hawkes$t)*exp(-temporal_params$rate * (time_window[2] - hawkes$t)))))
  } else{
    # stop("The shape and scale are going negative for power law")
    temporal_param_updates <- .optim_named(.temporal_parameter_likelihood, parameters, "temporal",
                                           fixed = fixed_temporal, hawkes = hawkes, time_diff = time_diff, parent_est_mat = parent_est_mat, triggering_rate = triggering_rate_update)
  }

  list(background_rate = background_rate_update,
       triggering_rate = triggering_rate_update,
       spatial = spatial_param_updates,
       temporal = temporal_param_updates)
}




#' Function to estimate MLEs
#'
#' @param hawkes A `hawkes` object.
#' @param inits Named list of initial background, triggering, spatial, and temporal
#'   parameters. To fix parameters, include a `fixed` list such as
#'   `fixed$spatial = c("mean")`.
#' @param boundary Optional boundary width used for edge correction.
#' @param max_iters Maximum number of EM iterations. Defaults to 500.
#' @param verbose Logical flag to print parameter updates during fitting.
#'
#' @returns A `hawkes_fit` object containing the MLEs.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0, 10, 0, 10)
#'
#' params <- list(
#'   background_rate = list(intercept = -4),
#'   triggering_rate = 0.5,
#'   spatial = list(mean = 0, sd = 0.1),
#'   temporal = list(rate = 2)
#' )
#' hawkes <- rHawkes(params, time_window = c(0, 100), spatial_region = spatial_region)
#' hawkes_mle(hawkes, inits = params)
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
#' hawkes_mle(hawkes, inits = params, boundary = 1)
hawkes_mle <- function(hawkes, inits, boundary = NULL, max_iters = 500, verbose = FALSE) {
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


  background_rate <- inits$background_rate
  triggering_rate <- inits$triggering_rate
  temporal_params <- inits$temporal
  spatial_params <- inits$spatial

  fixed_spatial <- if (!is.null(inits$fixed[["spatial"]])) inits$fixed[["spatial"]] else NULL
  fixed_temporal <- if (!is.null(inits$fixed[["temporal"]])) inits$fixed[["temporal"]] else NULL

  # Get initial estimate of the parent matrix
  # Function returns both the parent matrix in parent_mat
  # and the conditional intensity in cond_intensity
  parent_est_mat <- parent_est(hawkes, inits)

  # Initial estimate of the parameters
  est1 <- est_params(hawkes, inits, parent_est_mat, boundary = boundary,
                     fixed_spatial = fixed_spatial, fixed_temporal = fixed_temporal)

  likelihood_1 <- log_likelihood(hawkes, est1)

  for (i in 1:max_iters) {
    parent_est_mat <- parent_est(hawkes, est1)

    est2 <- est_params(hawkes, est1, parent_est_mat, boundary = boundary,
                       fixed_spatial = fixed_spatial, fixed_temporal = fixed_temporal)

    if (verbose) {
      print(est2)
    }

    likelihood_2 <- log_likelihood(hawkes, est2)

    # Check if the full likelihood has converged
    if (abs(likelihood_2 - likelihood_1) < 10e-10) {
      est2$fixed <- inits$fixed
      return(new_hawkes_fit(hawkes, est2))
    }
    est1 <- est2
    likelihood_1 <- likelihood_2
  }

  est2$fixed <- inits$fixed
  new_hawkes_fit(hawkes, est2)
}




#' Estimate the covariance matrix of the MLE using observed information
#'
#' Computes the observed information matrix (negative Hessian of the log-likelihood) and
#' returns its inverse as an estimated covariance matrix for Hawkes model parameters.
#'
#' @param hawkes A `hawkes` object.
#' @param est A `hawkes_fit` object from `hawkes_mle()`.
#'
#' @return A named covariance matrix corresponding to the estimated parameters.
#'
#' @details Requires the \pkg{numDeriv} package. The Hessian of the full log-likelihood is
#'   approximated with central differences. The negative Hessian is inverted and its row and
#'   column names are matched to the flattened parameter vector.
#'
#' @importFrom numDeriv hessian
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
#' hessian_est(hawkes, est$est)
#'
hessian_est <- function(hawkes, est) {
  est_vec <- .flatten_free_params(est)

  hessian <- numDeriv::hessian(.vector_input_log_likelihood, est_vec, hawkes = hawkes, param_template = est)
  cov_est <- solve(-1 * hessian)

  colnames(cov_est) <- names(est_vec)
  rownames(cov_est) <- names(est_vec)

  cov_est
}


