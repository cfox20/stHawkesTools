#' Construct a hawkes_fit object from estimated parameters
#'
#' @param hawkes The original `hawkes` object used to fit the model.
#' @param est Named nested list of estimated parameters, e.g., from `hawkes_mle()`.
#'
#' @return An object of class `hawkes_fit`.
#' @export
new_hawkes_fit <- function(hawkes, est) {
  stopifnot(is.list(est))

  structure(
    list(
      est = est,
      estimate = est,
      hawkes = hawkes,
      hawkes_object = hawkes,
      residuals = time_scaled_residuals(hawkes, est)
    ),
    class = c("hawkes_fit", class(est))
  )
}

#' Wald confidence intervals for Hawkes model parameters
#'
#' Computes Wald confidence intervals for Hawkes model parameters using the observed
#' information (negative Hessian) as a covariance estimate.
#'
#' @param object A `hawkes_fit` object.
#' @param parm Optional parameter subset (currently unused).
#' @param level Confidence level for the interval. Default is 0.95.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A tibble with point estimates and lower/upper bounds using `confint()` style
#'   column names.
#'
#' @details Fixed parameters, if any, are excluded when computing intervals.
#'
#' @importFrom stats confint
#'
#' @seealso [hessian_est()], [confint()]
#' @export
#'
#' @examples
#' set.seed(123)
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
#' confint(est)
#'
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
#' est <- hawkes_mle(hawkes, inits = params)
#' confint(est)
confint.hawkes_fit <- function(object, parm = NULL, level = 0.95, ...) {
  message("Asymptotic confidence inervals are known to undercover. Consider bootstrap methods.")

  hawkes <- if (!is.null(object$hawkes_object)) object$hawkes_object else object$hawkes

  est <- if (!is.null(object$estimate)) object$estimate else object$est

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
  est <- if (!is.null(x$estimate)) x$estimate else x$est
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
#' @param object a hawkes fit object to be summarized.
#' @param level Confidence level used for interval estimates. Default is 0.95.
#' @param digits Minimum number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `summary.hawkes_fit` containing the coefficient table,
#'   residual summary, and the confidence level used.
#'
#' @export
#'
summary.hawkes_fit <- function(object, level = 0.95, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(object, "hawkes_fit")) {
    stop("`object` must be a `hawkes_fit`.")
  }

  est <- if (!is.null(object$estimate)) object$estimate else object$est

  if (is.null(est)) {
    stop("Hawkes fit object does not contain parameter estimates.")
  }

  conf_tbl <- confint(object, level = level, ...)

  if (is.null(object$residuals) && !is.null(object$hawkes_object)) {
    residuals <- time_scaled_residuals(object$hawkes_object, est)
  } else {
    residuals <- object$residuals
  }

  if (!is.null(residuals)) {
    residual_summary <- stats::quantile(residuals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    names(residual_summary) <- c("Min", "1Q", "Median", "3Q", "Max")
  } else {
    residual_summary <- NULL
  }

  conf_df <- as.data.frame(conf_tbl)

  parts <- stringr::str_split_fixed(conf_df$Variable, "\\.", 2)

  component_raw <- parts[, 1]
  term_raw <- ifelse(parts[, 2] == "", parts[, 1], parts[, 2])

  label_fun <- function(x) {
    tools::toTitleCase(gsub("_", " ", x, fixed = TRUE))
  }

  row_labels <- ifelse(
    component_raw == term_raw,
    label_fun(component_raw),
    paste(label_fun(component_raw), label_fun(term_raw), sep = ": ")
  )

  component_levels <- c("background_rate", "triggering_rate", "spatial", "temporal")
  component_rank <- match(component_raw, component_levels)
  if (any(is.na(component_rank))) {
    next_rank <- max(component_rank, na.rm = TRUE)
    if (!is.finite(next_rank)) {
      next_rank <- 0
    }
    component_rank[is.na(component_rank)] <- next_rank + seq_len(sum(is.na(component_rank)))
  }

  order_idx <- order(component_rank, row_labels)
  conf_df <- conf_df[order_idx, , drop = FALSE]
  row_labels <- row_labels[order_idx]

  alpha <- 1 - level
  lower_col <- paste0(formatC(100 * alpha / 2, format = "f", digits = 1), "%")
  upper_col <- paste0(formatC(100 * (1 - alpha / 2), format = "f", digits = 1), "%")

  coef_mat <- as.matrix(conf_df[, c("Estimate", "std_error", "Lower", "Upper")])
  colnames(coef_mat) <- c("Estimate", "Std. Error", lower_col, upper_col)
  rownames(coef_mat) <- row_labels

  n_obs <- NA_integer_
  if (!is.null(object$hawkes_object)) {
    n_obs <- nrow(object$hawkes_object)
  } else if (!is.null(object$hawkes)) {
    n_obs <- nrow(object$hawkes)
  }

  out <- list(
    coefficients = coef_mat,
    conf.level = level,
    residual_summary = residual_summary,
    digits = digits,
    n = n_obs
  )

  class(out) <- "summary.hawkes_fit"
  out
}

#' @rdname summary.hawkes_fit
#' @param x A `summary.hawkes_fit` object produced by `summary.hawkes_fit()`.
#' @param ... Further arguments passed to or from other methods.
#' @return For `print.summary.hawkes_fit()`, returns `x` (invisibly).
#' @export
print.summary.hawkes_fit <- function(x, digits = x$digits, ...) {
  if (!is.null(x$n) && is.finite(x$n)) {
    cat(sprintf("Number of events: %s\n\n", x$n))
  }

  if (!is.null(x$residual_summary)) {
    cat("Residuals:\n")
    print.default(x$residual_summary, digits = digits)
    cat("\n")
  }

  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, has.Pvalue = FALSE)
  cat(sprintf("\nConfidence level used: %.1f%%\n", 100 * x$conf.level))


  invisible(x)
}
#





































