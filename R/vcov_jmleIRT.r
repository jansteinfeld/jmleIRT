#' Approximate covariance matrix for jmleIRT objects
#'
#' Uses diagonal JMLE-based standard errors and returns a block-diagonal
#' covariance matrix (person block and item block).
#'
#' @param object An object of class \code{"jmleIRT"}.
#' @param ... Currently unused.
#' @return A covariance matrix with row/column names.
#' @export
vcov.jmleIRT <- function(object, ...) {
  stopifnot(inherits(object, "jmleIRT"))
  se_theta <- object$se_theta
  se_beta <- object$se_beta

  if (is.null(se_theta) || is.null(se_beta)) {
    stop("Standard errors not available in object.")
  }

  V_theta <- diag(se_theta^2)
  V_beta <- diag(se_beta^2)
  V <- Matrix::bdiag(V_theta, V_beta)
  rn <- c(
    paste0("theta_", seq_along(se_theta)),
    paste0("beta_", seq_along(se_beta))
  )
  colnames(V) <- rownames(V) <- rn
  as.matrix(V)
}
