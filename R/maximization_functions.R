
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
  .unpack_hawkes(hawkes)

  -sum({parent_est_mat *
      (.safe_log(do.call(spatial_pdf, c(list(x = x_diff), p))) +
         .safe_log(do.call(spatial_pdf, c(list(x = y_diff), p))))
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
.temporal_parameter_likelihood <- function(p, hawkes, parent_est_mat, time_diff) {
  .unpack_hawkes(hawkes)
  triggering_rate <- params$triggering_rate

  -{sum(parent_est_mat *
          .safe_log(do.call(temporal_pdf, c(list(x = time_diff), p)))) -
      triggering_rate *
      sum(do.call(temporal_cdf, c(list(q = region$t[2] - hawkes$t), p)))
  }
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

  opt <- optim(par = par_vec, fn = wrapped_fn, method = "BFGS")

  if (opt$convergence != 0) {
    warning(sprintf("Optimization did not converge for %s: %s", component, opt$message))
  } else if (verbose) {
    message(sprintf("Optimization converged for %s", component))
  }

  param_list[[component]][to_opt] <- as.list(opt$par)
  return(param_list[[component]])
}



background_covariates_function <- function(background_rate, hawkes, parent_est_mat, X, cov_map, region) {
  cov_names <- colnames(cov_map)[!(colnames(cov_map) %in% c("geoid", "name", "area", "geometry"))]

  cov_map_X <- cov_map |>
    sf::st_drop_geometry() |>
    dplyr::select(all_of(cov_names)) |>
    as.matrix()
  cov_map_X <- cbind(1, cov_map_X)
  area <- cov_map$area
  t_length <- region$t[[2]] - region$t[[1]]

  background_term <- as.numeric(t(diag(parent_est_mat)) %*% X)

  intensity_term <- t_length * colSums(area * cov_map_X * as.numeric(exp(cov_map_X %*% background_rate)))

  background_term - intensity_term
}
