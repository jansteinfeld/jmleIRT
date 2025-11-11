#' Summary of JMLE Rasch Model Estimation
#'
#' Print a summary of the Joint Maximum Likelihood Estimation (JMLE) Rasch model results,
#' including number of iterations, item difficulties, and person ability parameters.
#'
#' @param object An object of class \code{"jmleIRT"} containing numeric vectors 
#'   \code{theta} and \code{beta}.
#' If not a strict fit object, it must still provide \code{theta} and 
#'   \code{beta} accessible as \code{object$theta} and \code{object$beta}.
#' @param digits Number of digits for numeric rounding (default 4).
#' @param ... Additional arguments passed to or from other methods (currently unused).
#'
#' @return Returns the input \code{object} invisibly.
#'
#' @examples
#' \dontrun{
#' \examples{
#' X <- matrix(
#'   c(1, 0, 1, NA, 0,
#'     1, 1, 1, 1, 0,
#'     0, NA, NA, 0, 1,
#'     0, 0, 0, NA),
#'   nrow = 5,
#'   byrow = TRUE
#' )
#' }
#' res <- jmle_estimation(X, max_iter = 200, estimatewle = TRUE, verbose = FALSE)
#' summary(res)
#' }
#'
#' @export
#' @method summary jmleIRT
summary.jmleIRT <- function(object, digits = 4, ...) {
  stopifnot(inherits(object, "jmleIRT"))

  theta <- object$theta
  beta <- object$beta
  wle <- object$wle_estimate
  conv <- object$converged
  iter <- object$iterations
  bias_corr <- object$bias_correction
  centered <- object$centered
  n_persons <- length(theta)
  n_items <- length(beta)

  # Helper for statistics printing
  print_stats <- function(x) {
    c(
      Mean = mean(x, na.rm = TRUE),
      SD = stats::sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  }

  cat("Joint Maximum Likelihood Estimation of Rasch Model\n")
  cat("-----------------------------------------------\n")
  cat("Number of persons:", n_persons, "\n")
  cat("Number of items :", n_items, "\n")
  cat("Converged :", conv, "\n")
  cat("Iterations :", iter, "\n")
  cat("Bias correction :", bias_corr, "\n")
  cat("Centered (mean beta=0):", centered, "\n\n")

  cat("Person parameters (theta):\n")
  print(round(print_stats(theta), digits = digits))

  cat("\nItem difficulties (beta):\n")
  print(round(print_stats(beta), digits = digits))

  if (!all(is.na(wle))) {
    cat("\nPerson WLE estimates:\n")
    print(round(print_stats(wle), digits = digits))
  }

  invisible(object)
}

#' Summary for Bias Corrected JMLE Estimates
#'
#' Print a summary of bias-corrected item difficulty and person ability estimates.
#'
#' @param object Object of class \code{"biasCorrection"} containing \code{corrected_b} and \code{corrected_theta}.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input \code{object}, invisibly.
#'
#' @export
summary.biasCorrection <- function(object, ...) {
  cat("Summary of Bias Corrected Estimates\n")
  cat("-----------------------------------\n")

  cat("Corrected item difficulties:\n")
  print(summary(object$corrected_b))

  cat("Corrected person abilities:\n")
  print(summary(object$corrected_theta))

  invisible(object)
}
