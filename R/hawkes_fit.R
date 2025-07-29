#' Construct a hawkes_fit object from estimated parameters
#'
#' @param est A named nested list of estimated parameters, e.g. from `hawkes_mle()`
#' @param hawkes The original hawkes object used to fit the model
#' @param residual_fn Function to compute residuals (default = `compute_residuals`)
#'
#' @return An object of class `hawkes_fit`
#' @export
new_hawkes_fit <- function(est) {
  stopifnot(is.list(est))

  structure(
    # list(
    #   estimate = est,
    #   hawkes = hawkes,
    #   residuals = residual_fn(hawkes, est)
    # ),
    est,
    class = c("hawkes_fit", class(est))
  )
}


#' Print hawkes fit object
#'
#' @param x a hawkes fit object to be printed
#'
#' @returns The input hawkes fit object, invisibly.
#' @export
#'
print.hawkes_fit <- function(x) {

  cat("Triggering Parameter Estimates:\n")

  cat("Background Rate (\u03B2):   \n")
  for (nm in names(x$background_rate)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(x$background_rate[[nm]], 3))))
  }
  cat(sprintf(" Triggering Rate (\u03b8):   %s\n", round(x$triggering_rate, 3)))

  cat("\nSpatial Triggering Parameter Estimates:\n")
  for (nm in names(x$spatial)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(x$spatial[[nm]], 3))))
  }

  cat("Temporal Triggering Parameter Estimates:\n")
  for (nm in names(x$temporal)) {
    cat(sprintf("    %s: %s\n", nm, toString(round(x$temporal[[nm]], 3))))
  }

  invisible(x)
}
