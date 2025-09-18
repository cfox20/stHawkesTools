#' Construct a hawkes_fit object from estimated parameters
#'
#' @param hawkes The original hawkes object used to fit the model
#' @param est A named nested list of estimated parameters, e.g. from `hawkes_mle()`
#'
#' @return An object of class `hawkes_fit`
#' @export
new_hawkes_fit <- function(hawkes, est) {
  stopifnot(is.list(est))

  structure(
    list(
      estimate = est,
      hawkes_object = hawkes,
      residuals = time_scaled_residuals(hawkes, est)
    ),
    # est,
    class = c("hawkes_fit", class(est))
  )
}

#' Wald Confidence Intervals for Hawkes Model Parameters
#'
#' Computes Wald marginal confidence intervals for the parameters of a spatio-temporal Hawkes process
#' using the estimated covariance matrix from the observed information (negative Hessian).
#'
#' @param object A `hawkes_fit` object
#' @param parm Parameters to create confidence interval for. (currently not used)
#' @param level Confidence level for the interval. Default is 0.95.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data frame with the point estimates and corresponding lower and upper confidence bounds.
#' The column names match the default format of `confint()`, e.g., `"2.5 %"`, `"97.5 %"`.
#'
#' @details The function uses the estimated Hessian of the log-likelihood to compute the covariance matrix.
#' Fixed parameters (if any) are excluded from the interval computation.
#'
#' @importFrom stats confint
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
confint.hawkes_fit <- function(object, parm = NULL, level = 0.95, ...) {
  hawkes <- object$hawkes

  est <- object$est

  alpha <- 1 - level

  # if (method == "rathbun") {
  #   cov_mat <- rathbun_est(hawkes, est, X)
  # }else{
  #   cov_mat <- hessian_est(hawkes, est, region)
  # }
  cov_est <- hessian_est(hawkes, est)

  est_vec <- .flatten_free_params(est)

  z <- -stats::qnorm((alpha)/2)

  lower_name <- paste0(formatC(100 * alpha / 2, format = "f", digits = 1), " %")
  upper_name <- paste0(formatC(100 * (1 - alpha / 2), format = "f", digits = 1), " %")

  tibble::tibble(
    Variable = names(est_vec),
    Estimate = est_vec,
    std_error = sqrt(diag(cov_est)),
    Lower = est_vec - z * sqrt(diag(cov_est)),
    Upper = est_vec + z * sqrt(diag(cov_est))
  )
}



#' Print hawkes fit object
#'
#' @param x a hawkes fit object to be printed
#' @param ... Further arguments passed to or from other methods.
#'
#' @returns The input hawkes fit object, invisibly.
#' @export
#'
print.hawkes_fit <- function(x, ...) {
  est <- x$est
  cat("Triggering Parameter Estimates:\n")

  cat("Background Rate (\u03B2):   \n")
  for (nm in names(est$background_rate)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(est$background_rate[[nm]], 3))))
  }
  cat(sprintf(" Triggering Rate (\u03b8):   %s\n", round(est$triggering_rate, 3)))

  cat("\nSpatial Triggering Parameter Estimates:\n")
  for (nm in names(est$spatial)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(est$spatial[[nm]], 3))))
  }

  cat("Temporal Triggering Parameter Estimates:\n")
  for (nm in names(est$temporal)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(est$temporal[[nm]], 3))))
  }

  invisible(est)
}


#' Print hawkes fit object
#'
#' @param object a hawkes fit object to be printed
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#'
summary.hawkes_fit <- function(object, ...){
  conf_int <- confint(object)

  cat("Background Coefficients:\n")
  conf_int |>
    dplyr::filter(stringr::str_detect(.data$Variable, "^background_rate")) |>
    dplyr::mutate(Variable = stringr::str_remove(.data$Variable, "^background_rate.")) |>
    dplyr::rename(Coefficient = .data$Variable,
                  `Std. Error` = .data$std_error,
                  `95% Lower Bound` = .data$Lower,
                  `95% Upper Bound` = .data$Upper) |>
    print()
  cat("\n")
  cat("Trigering Coefficients:\n")
  conf_int |>
    dplyr::filter(!stringr::str_detect(.data$Variable, "^background_rate")) |>
    dplyr::mutate(Variable = stringr::str_remove(.data$Variable, "^background_rate.")) |>
    dplyr::rename(Coefficient = .data$Variable,
                  `Std. Error` = .data$std_error,
                  `95% Lower Bound` = .data$Lower,
                  `95% Upper Bound` = .data$Upper) |>
    print()
}
#





































