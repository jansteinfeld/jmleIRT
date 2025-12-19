#' Coefficients for jmleIRT objects
#'
#' @param object An object of class \code{"jmleIRT"}.
#' @param ... Currently unused.
#' @return A list with components \code{theta} and \code{beta}.
#' @export
coef.jmleIRT <- function(object, ...) {
  stopifnot(inherits(object, "jmleIRT"))
  list(theta = object$theta, beta = object$beta)
}
