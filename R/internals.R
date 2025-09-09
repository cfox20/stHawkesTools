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
#' @keywords internal
.sanity_check <- function(hawkes) {
  if (class(hawkes)[[1]] != "hawkes") {
    stop("Object must be a hawkes object.")
  }

  required_attrs <- c("spatial_pdf", "temporal_pdf", "spatial_cdf", "temporal_cdf", "spatial_region")

  missing_attrs <- required_attrs[!required_attrs %in% names(attributes(hawkes))]

  if (length(missing_attrs) > 0) {
    stop(glue::glue("The following required attributes are missing from the hawkes object: {paste(missing_attrs, collapse = ', ')}"))
  }

  invisible(TRUE)
}


#' Flatten a nested parameter list into a numeric vector of free parameters with names
#'
#' @param params A nested list of parameter values, including a `$fixed` sublist.
#'
#' @returns A list with two elements:
#'   \describe{
#'     \item{`values`}{Numeric vector of non-fixed parameter values.}
#'     \item{`names`}{Character vector of names, using dot notation for nested entries.}
#'   }
#' @keywords internal
.flatten_free_params <- function(params) {
  fixed <- params$fixed
  params$fixed <- NULL

  values <- c()
  names_out <- c()

  for (name1 in names(params)) {
    if (is.list(params[[name1]])) {
      fixed_sub <- fixed[[name1]] %||% character(0)
      for (name2 in names(params[[name1]])) {
        if (!(name2 %in% fixed_sub)) {
          values <- c(values, params[[name1]][[name2]])
          names_out <- c(names_out, paste0(name1,".", name2))
        }
      }
    } else {
      if (!(name1 %in% names(fixed))) {
        values <- c(values, params[[name1]])
        names_out <- c(names_out, rep(name1, length(params[[name1]])))
      }
    }
  }

  names(values) <- names_out
  values
}


#' Reconstruct parameter list from flat vector and template with fixed values
#'
#' @param flat_vec Named numeric vector of free parameter values
#' @param template Parameter template list, must contain `$fixed`
#'
#' @return Full parameter list including fixed values, matching `template` structure
#' @keywords internal
.vector_to_params <- function(flat_vec, template) {
  fixed <- template$fixed
  fixed_keys <- paste(names(fixed), fixed, sep = ".")
  template$fixed <- NULL

  # Collapse duplicated names
  collapsed_vec <- split(flat_vec, names(flat_vec))

  # Helper: remove names from numeric vectors
  unname_vec <- function(x) {
    x <- as.numeric(x)
    unname(x)
  }

  new_params <- list()

  for (param_type in names(template)) {
    if (is.numeric(template[[param_type]])) {
      new_params[[param_type]] <- unname_vec(collapsed_vec[[param_type]])
    }

    if (is.list(template[[param_type]])) {
      new_param_group <- list()

      for (param in names(template[[param_type]])) {
        full_name <- paste0(param_type, ".", param)

        if (full_name %in% fixed_keys) {
          new_param_group[[param]] <- unname_vec(template[[param_type]][[param]])
        } else {
          new_param_group[[param]] <- unname_vec(collapsed_vec[[full_name]])
        }
      }

      new_params[[param_type]] <- new_param_group
    }
  }

  new_params
}



#' Convert Hawkes Parameter List to Data Frame
#'
#' Converts a fitted Hawkes parameter list into a tidy data frame with one row per parameter.
#'
#' @param est A list of fitted parameter values from a Hawkes model, typically returned by `hawkes_mle()`. Must include components like `background_rate`, `triggering_rate`, `spatial`, and `temporal`.
#'
#' @returns A data frame
#'
#' @keywords internal
.hawkes_mle_to_dataframe <- function(est) {
  if (class(est)[1] == "hawkes_fit") {
    est <- est$est
  }
  est <- est[names(est) != "fixed"]

  purrr::imap_dfr(est, function(value, type) {
    if (is.list(value)) {
      flat <- unlist(value)
    } else {
      flat <- setNames(value, type)
    }

    tibble::tibble(
      parameter_type = type,
      parameter = names(flat),
      value = as.numeric(flat)
    )
  })
}


#' Wrapper for likelihood for numerical optimization
#'
#' @param par_vec a numeric vector of the parameter initial values
#' @param hawkes a hawkes object
#' @param param_template a template to reconstruct the parameter object
#'
#' @returns
#' @export
#'
#' @keywords internal
.vector_input_log_likelihood <- function(par_vec, hawkes, param_template) {
  parameters <- .vector_to_params(par_vec, param_template)
  log_likelihood(hawkes, parameters)
}















