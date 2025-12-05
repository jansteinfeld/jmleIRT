#' @export summary.jmleIRT
#' @export
summary.jmleIRT <- function(object, ...) {
  # Print summary of estimation results
  cat("Summary of JMLE Rasch Model Estimation\n")
  cat("-------------------------------------\n")
  cat(sprintf("Number of iterations: %d\n", object$iter))
  cat("Item difficulties (b):\n")
  print(summary(object$b))
  cat("Person abilities (theta):\n")
  print(summary(object$theta))

  invisible(object)
}

#' @export summary.biasCorrection
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
