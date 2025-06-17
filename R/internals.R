#' (Internal) Unpack all Attributes in Hawkes Object
#'
#' @param hawkes A hawkes object
#'
#' @returns Loads attributes to function environment
#'
#'
.unpack_hawkes <- function(hawkes) {
  list2env(as.list(hawkes), envir = parent.frame()) |> invisible()
  list2env(attributes(hawkes), envir = parent.frame()) |> invisible()
}


#' (Internal) Check to make sure necessary attributes are in the Hawkes object
#'
#' @param hawkes A hawkes object
#'
#' @returns TRUE
#'
.sanity_check <- function(hawkes) {
  if (class(hawkes)[[1]] != "hawkes") {
    stop("Object must be a hawkes object.")
  }

  required_attrs <- c("spatial_pdf", "temporal_pdf", "spatial_cdf", "temporal_cdf", "spatial_sampler", "temporal_sampler", "params", "region")

  missing_attrs <- required_attrs[!required_attrs %in% names(attributes(hawkes))]

  if (length(missing_attrs) > 0) {
    stop(glue::glue("The following required attributes are missing from the hawkes object: {paste(missing_attrs, collapse = ', ')}"))
  }

  invisible(TRUE)
}


