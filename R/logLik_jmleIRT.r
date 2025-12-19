#' Log-likelihood for jmleIRT objects
#'
#' Computes the joint log-likelihood of the Rasch model at the JMLE estimates.
#'
#' @param object An object of class \code{"jmleIRT"}.
#' @param ... Currently unused.
#' @return An object of class \code{"logLik"}.
#' @export
logLik.jmleIRT <- function(object, ...) {
  stopifnot(inherits(object, "jmleIRT"))
  X <- object$data
  theta <- object$theta
  beta <- object$beta

  ll <- 0
  for (p in seq_len(nrow(X))) {
    for (i in seq_len(ncol(X))) {
      x <- X[p, i]
      if (is.na(x)) next
      eta <- theta[p] - beta[i]
      p_ni <- 1 / (1 + exp(-eta))
      if (x == 1) {
        ll <- ll + log(p_ni)
      } else {
        ll <- ll + log(1 - p_ni)
      }
    }
  }
  attr(ll, "df") <- length(theta) + length(beta)
  attr(ll, "nobs") <- sum(!is.na(X))
  class(ll) <- "logLik"
  ll
}
