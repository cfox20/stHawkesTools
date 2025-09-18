
#' Safe logarithm to avoid log(0)
#'
#' Computes the natural logarithm of a numeric vector, replacing values too close to zero with a small positive constant to avoid `-Inf` or `NaN` results.
#'
#' @param x A numeric vector.
#'
#' @returns A numeric vector of the same length as `x` with `log(pmax(x, .Machine$double.eps))`.
#' @keywords internal
.safe_log <- function(x) log(pmax(x, .Machine$double.eps))

#' Spatial triggering kernel parameter likelihood
#'
#' Computes the (negative) expected complete-data log-likelihood for spatial triggering kernel parameters, assuming separable spatial kernels.
#'
#' @param p A named list of parameters for the spatial kernel.
#' @param hawkes A `hawkes` object.
#' @param parent_est_mat Matrix of estimated parent probabilities.
#' @param x_diff A matrix of pairwise x-coordinate differences.
#' @param y_diff A matrix of pairwise y-coordinate differences.
#'
#' @returns The negative log-likelihood contribution of the spatial triggering parameters.
#' @keywords internal
.spatial_parameter_likelihood <- function(p, hawkes, parent_est_mat, x_diff, y_diff) {
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


  -sum({parent_est_mat *
      (.safe_log(do.call(.env$spatial_pdf, c(list(x = .env$x_diff), p))) +
         .safe_log(do.call(.env$spatial_pdf, c(list(x = .env$y_diff), p))))
  })
}


#' Temporal triggering kernel parameter likelihood
#'
#' Computes the (negative) expected complete-data log-likelihood for temporal triggering kernel parameters.
#'
#' @param p A named list of parameters for the temporal kernel.
#' @param hawkes A `hawkes` object.
#' @param parent_est_mat Matrix of estimated parent probabilities.
#' @param time_diff A matrix of pairwise time differences.
#'
#' @returns The negative log-likelihood contribution of the temporal triggering parameters.
#' @keywords internal
.temporal_parameter_likelihood <- function(p, hawkes, parent_est_mat, time_diff, triggering_rate) {
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


  tryCatch(-{sum(parent_est_mat *
          .safe_log(do.call(.env$temporal_pdf, c(list(x = .env$time_diff), p)))) -
      triggering_rate *
      sum(do.call(.env$temporal_cdf, c(list(q = .env$time_window[2] - .env$hawkes$t), .env$p)))},
    error = function(e){
      stop(paste("Error in temporal parameter optimization:", e$message, "\n Last parameter values: ", p, "\n"))
    })
}


#' Optimize kernel parameters with named list support
#'
#' Optimizes a likelihood function using `optim()` with parameters passed as a named list. Only non-fixed parameters are estimated.
#'
#' @param fn The objective function to minimize, taking a list of named parameters.
#' @param param_list A named list of parameter lists (e.g., with components `spatial`, `temporal`).
#' @param component A character string indicating which component of `param_list` to optimize.
#' @param fixed A character vector of parameter names to keep fixed during optimization.
#' @param verbose Logical; if `TRUE`, prints convergence messages and allows warnings to print.
#' @param ... Additional arguments passed to `fn`.
#'
#' @returns A named list with updated parameter estimates for the specified component.
#' @keywords internal
.optim_named <- function(fn, param_list, component, fixed = character(), verbose = FALSE, ...) {
  sublist <- param_list[[component]]
  to_opt <- setdiff(names(sublist), fixed)
  par_vec <- unlist(sublist[to_opt])

  wrapped_fn <- function(par_vec_input) {
    updated <- sublist
    updated[to_opt] <- as.list(par_vec_input)

    withCallingHandlers(
      fn(updated, ...),
      warning = function(w) {
        if (!verbose) invokeRestart("muffleWarning")
      }
    )
  }

  opt <- stats::optim(par = par_vec, fn = wrapped_fn, method = "BFGS")

  if (opt$convergence != 0) {
    warning(sprintf("Optimization did not converge for %s: %s", component, opt$message))
  } else if (verbose) {
    message(sprintf("Optimization converged for %s", component))
  }

  param_list[[component]][to_opt] <- as.list(opt$par)
  return(param_list[[component]])
}



background_covariates_function <- function(background_rate, hawkes, parent_est_mat, time_window, spatial_region, covariate_columns, X) {
  # Make matrix of background covariate values for each region
  cov_map_X <- spatial_region |>
    sf::st_drop_geometry() |>
    dplyr::select(tidyselect::all_of(covariate_columns)) |>
    as.matrix()

  cov_map_X <- cbind(1, cov_map_X)
  area <- spatial_region$area
  t_length <- time_window[2] - time_window[1]

  background_term <- as.numeric(t(diag(parent_est_mat)) %*% X)

  intensity_term <- t_length * colSums(area * cov_map_X * as.numeric(exp(cov_map_X %*% background_rate)))

  background_term - intensity_term
}
