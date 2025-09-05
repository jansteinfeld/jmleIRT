#' Summary of JML estimation
#'
#' Print a brief summary for a JML fit, showing means of person and item parameters.
#' Expects an object with numeric components `theta` and `beta`; designed for objects of class "jmleIRT".
#'
#' @param object An object of class "jmleIRT" containing numeric vectors `theta` and `beta`. If not a strict fit object, it must still provide `theta` and `beta` accessible as `object$theta` and `object$beta`.
#' @param digits Number of digits to display (default is 4).
#' @param ... Further arguments passed to or used by methods (currently unused).
#'
#' @return The input object, returned invisibly.
#'
#' # Minimal example with an adapter object
#' X <- matrix(c(
#'  1, 0, 1, NA,
#'  0, 1, 1, 1,
#'  1, 1, 0, 0,
#'  NA, NA, 0, 1,
#'  0, 0, 0, NA
#'), nrow=5, byrow=TRUE)
#'
#' res <- jmle_estimation(X, max_iter=200, estimatewle=TRUE, verbose=FALSE)
#' summary(res)
#' @importFrom stats sd
#' @seealso [jmle_estimation]
#' @export
#' @method summary jmleIRT
summary.jmleIRT <- function(object, digits = 4, ...) {
  stopifnot(inherits(object, "jmleIRT"))
  
  # Extract components
  theta <- object$theta
  beta <- object$beta
  wle <- object$wle_estimate
  conv <- object$converged
  iter <- object$iterations
  bias_corr <- object$bias_correction
  centered <- object$centered
  
  n_persons <- length(theta)
  n_items <- length(beta)
  
  # Helper function: print stats
  print_stats <- function(x) {
    c(
      Mean = mean(x, na.rm = TRUE),
      SD = sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  }
  
  cat("Joint Maximum Likelihood Estimation of Rasch Model\n")
  cat("-----------------------------------------------\n")
  cat("Number of persons:", n_persons, "\n")
  cat("Number of items  :", n_items, "\n")
  cat("Converged       :", conv, "\n")
  cat("Iterations      :", iter, "\n")
  cat("Bias correction :", bias_corr, "\n")
  cat("Centered (mean beta=0):", centered, "\n\n")
  
  cat("Person parameters (theta):\n")
  print(round(print_stats(theta), digits=digits))
  cat("\nItem difficulties (beta):\n")
  print(round(print_stats(beta), digits=digits))
  
  if (!all(is.na(wle))) {
    cat("\nPerson WLE estimates:\n")
    print(round(print_stats(wle), digits=digits))
  }
  invisible(object)
}
