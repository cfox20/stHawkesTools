#' Gaussian spatial parameter update
#'
#'  This is an internal helper and not
#' meant to be called directly by users.
#'
#' @returns A list of updated sd estimate
#'
#' @keywords internal
.spatial_update_gaussian <- function() {
  list(mean = 0, sd = sqrt(sum(parent_est_mat * (x_diff^2 + y_diff^2)) / (2 * sum(parent_est_mat))))
}


#' Exponential temporal parameter update
#'
#'  This is an internal helper and not
#' meant to be called directly by users.
#'
#' @returns A list of the updated rate estimate
#'
#' @keywords internal
.temporal_update_exponential <- function(hawkes) {
  list(rate = sum(parent_est_mat) / (sum(parent_est_mat * time_diff) + triggering_rate_update * sum((t_max - hawkes$t)*exp(-temporal_params$rate * (t_max - hawkes$t)))))
}


# #' Power law parameter update
# #'
# #' This is an internal helper and not
# #' meant to be called directly by users.
# #'
# #' @returns A list of the updated rate estimate
# #'
# #' @keywords internal
# .temporal_update_exponential <- function(hawkes) {
#   list(rate = sum(parent_est_mat) / (sum(parent_est_mat * time_diff) + triggering_rate_update * sum((t_max - hawkes$t)*exp(-temporal_params$rate * (t_max - hawkes$t)))))
# }
