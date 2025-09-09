
#' Parent Matrix Estimation
#'
#' @param hawkes A `hawkes` object
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function.
#'
#' @returns A numeric matrix
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#'
#' (parent_est_mat <- parent_est(hawkes, params))
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 0)
parent_est <- function(hawkes, parameters) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  if (class(parameters)[1] == "hawkes_fit") {
    parameters <- parameters$est
  }

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if(!exists("covariate_columns", inherits = FALSE)){
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

#' Compute Parameter Estimate Updates
#'
#' @param hawkes A `hawkes` object
#' @param parameters A named list of lists containing the values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function. If some parameters in the kernel function should be fixed in estimation, pass their names as another named list within params as fixed$spatial = c("mean").
#' @param parent_est_mat A matrix produced from the parent_est_mat() function.
#' @param boundary A boundary region to correct for the bpundary bias. Defaults to NULL if not used.
#'
#' @returns A named list of the form of parameters with udpated parameter estimates.
#' @export
#'
#' @examples
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.5),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' parent_est_mat <- parent_est(hawkes, params)
#' est_params(hawkes, params, parent_est_mat)
#'
est_params <- function(hawkes, parameters, parent_est_mat, boundary = NULL, fixed_spatial = NULL, fixed_temporal = NULL) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")
  if(class(parent_est_mat)[[1]] != "matrix") stop("parent_est_mat must be a matrix")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

  if (is.null(parameters)) {
    message("True parameters stored in Hawkes object are being used for parent estimation. If estimated parameters are desired, pass the estimates to the parameters argument.")
    parameters <- params
  }

  if(!exists("covariate_columns", inherits = FALSE)){
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
    spatial_region <- spatial_region |>
      sf::st_buffer(-boundary)


    # Find the index of values that fall within the boundary region S0
    outside_S <- lengths(sf::st_within(hawkes, spatial_region_boundary)) <= 0
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


  if (spatial_family == "Gaussian" && temporal_family == "Exponential") {
    # triggering_rate_update <- sum(parent_est_mat) / (sum(1 - exp(-temporal_params$rate * (t_max - hawkes$t))))

    temporal_param_updates <- list(rate = sum(parent_est_mat) / (sum(parent_est_mat * time_diff) + triggering_rate_update * sum((time_window[2] - hawkes$t)*exp(-temporal_params$rate * (time_window[2] - hawkes$t)))))

    spatial_param_updates <- list(mean = 0, sd = sqrt(sum(parent_est_mat * (x_diff^2 + y_diff^2)) / (2 * sum(parent_est_mat))))

  } else {
  # Numerically solve for MLE for general kernels
    temporal_param_updates <- .optim_named(.temporal_parameter_likelihood, parameters, "temporal",
                                           fixed = fixed_temporal, hawkes = hawkes, time_diff = time_diff, parent_est_mat = parent_est_mat)

    spatial_param_updates <- .optim_named(.spatial_parameter_likelihood, parameters, "spatial", fixed = fixed_spatial, hawkes = hawkes,
                                          parent_est_mat = parent_est_mat, x_diff = x_diff, y_diff = y_diff)
  }


  list(background_rate = background_rate_update,
       triggering_rate = triggering_rate_update,
       spatial = spatial_param_updates,
       temporal = temporal_param_updates)
}




#' Function to Estimate MLEs
#'
#' @param hawkes A `hawkes` object
#' @param inits A named list of lists containing the initial values for the background rate, triggering ratio, spatial parameters in a named list, and temporal parameters in a named list. Note that these values are of the same form as the true values in the Hawkes object but are often estimates passed to the function. If some parameters in the kernel function should be fixed in estimation, pass their names as another named list within params as fixed$spatial = c("mean").
#' @param boundary A boundary region to correct for the bpundary bias. Defaults to NULL if not used.
#' @param max_iters A numeric value for the maximum number of iteration in the EM-algorithm. Defaults to 500 if not used.
#'
#' @returns A `hawkes_fit` object that is a list of lists containing the MLEs.
#' @export
#'
#' @examples
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#' hawkes_mle(hawkes, inits = params)
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 0)
#' hawkes_mle(hawkes, inits = params)
hawkes_mle <- function(hawkes, inits, boundary = NULL, max_iters = 500, verbose = FALSE) {
  if(class(hawkes)[1] != "hawkes") stop("hawkes must be a hawkes object")

  .sanity_check(hawkes)

  .unpack_hawkes(hawkes)

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




#' Estimate the covariance matrix of the MLE using the observed information
#'
#' Computes the observed information matrix (negative Hessian of the log-likelihood)
#' and returns its inverse as an estimate of the covariance matrix for the MLE parameters
#' in a spatio-temporal Hawkes process model.
#'
#' @param hawkes An object containing the data and model specification for a Hawkes process,
#'   passed to the log-likelihood function.
#' @param est A list of estimated parameter values (typically the output of an EM or MLE routine),
#'   used as the starting point for computing the Hessian. Must match the structure expected by
#'   `.flatten_free_params()` and `.vector_input_log_likelihood()`.
#'
#' @return A named covariance matrix corresponding to the estimated parameters.
#'
#' @details This function requires the \pkg{numDeriv} package. It numerically approximates the
#' Hessian of the full log-likelihood using central differences and inverts the negative Hessian
#' to obtain the covariance estimate. The names of the matrix rows and columns are set to match
#' the flattened parameter vector.
#'
#' @importFrom numDeriv hessian
#' @export
#'
#' @examples
#' set.seed(123)
#' region <- list(x = c(0,10), y = c(0,10), t = c(0,50))
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.5,spatial = list(mean = 0, sd = 0.1),temporal = list(rate = 2), fixed = list(spatial = "mean", temporal = NULL))
#' hawkes <- rHawkes(params, region)
#' est <- hawkes_mle(hawkes, inits = params, boundary = c(.5, 5))
#' hessian_est(hawkes, est)
#'
hessian_est <- function(hawkes, est) {
  est_vec <- .flatten_free_params(est)

  hessian <- numDeriv::hessian(.vector_input_log_likelihood, est_vec, hawkes = hawkes, param_template = est)
  cov_est <- solve(-1 * hessian)

  colnames(cov_est) <- names(est_vec)
  rownames(cov_est) <- names(est_vec)

  cov_est
}


#' Wald Confidence Intervals for Hawkes Model Parameters
#'
#' Computes Wald-style confidence intervals for the parameters of a spatio-temporal Hawkes process
#' using the estimated covariance matrix from the observed information (negative Hessian).
#'
#' @param hawkes_fit A `hawkes` object containing the observed point pattern and model setup.
#' @param est A nested list of estimated parameters, typically the output from `hawkes_mle()`.
#' @param conf_level Confidence level for the interval, e.g., 0.95 for 95% intervals. Default is 0.95.
#'
#' @return A data frame with the point estimates and corresponding lower and upper confidence bounds.
#' The column names match the default format of `confint()`, e.g., `"2.5 %"`, `"97.5 %"`.
#'
#' @details The function uses the estimated Hessian of the log-likelihood to compute the covariance matrix.
#' Fixed parameters (if any) are excluded from the interval computation.
#'
#' @seealso [hessian_est()], [confint()]
#' @export
#'
#' @examples
#' set.seed(123)
#' spatial_region <- create_rectangular_sf(0,10,0,10)
#'
#' params <- list(background_rate = list(intercept = -4),triggering_rate = 0.75,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2))
#' hawkes <- rHawkes(params, time_window = c(0,50), spatial_region = spatial_region)
#' est <- hawkes_mle(hawkes, inits = params)
#' confint(est)
#'
#'
#' params <- list(background_rate = list(intercept = -4.5, X1 = 1, X2 = 1),triggering_rate = 0.5,spatial = list(mean = 0, sd = .75),temporal = list(rate = 2), fixed = list(spatial = "mean"))
#' data("example_background_covariates")
#' hawkes <- rHawkes(params, c(0,50), example_background_covariates, covariate_columns = c("X1", "X2"), spatial_burnin = 0)
#' est <- hawkes_mle(hawkes, inits = params)
#' confint(est)
confint.hawkes_fit <- function(hawkes_fit, conf_level = 0.95) {
  hawkes <- hawkes_fit$hawkes

  est <- hawkes_fit$est

  alpha <- 1 - conf_level

  # if (method == "rathbun") {
  #   cov_mat <- rathbun_est(hawkes, est, X)
  # }else{
  #   cov_mat <- hessian_est(hawkes, est, region)
  # }
  cov_est <- hessian_est(hawkes, est)

  est_vec <- .flatten_free_params(est)

  z <- -qnorm((alpha)/2)

  lower_name <- paste0(formatC(100 * alpha / 2, format = "f", digits = 1), " %")
  upper_name <- paste0(formatC(100 * (1 - alpha / 2), format = "f", digits = 1), " %")

  data.frame(
    Variable = names(est_vec),
    Estimate = est_vec,
    setNames(list(est_vec - z * sqrt(diag(cov_est))), lower_name),
    setNames(list(est_vec + z * sqrt(diag(cov_est))), upper_name),
    check.names = FALSE,
    row.names = NULL
  )
}
