#' Print method for JMLE Rasch Model Objects
#'
#' This function prints a brief summary of a Joint Maximum Likelihood Estimation (JMLE)
#' Rasch model object of class \code{jmleIRT}.
#'
#' @param x An object of class \code{jmleIRT} returned from a JMLE Rasch model estimation function.
#' @param ... Additional arguments passed to or from other methods (currently ignored).
#'
#' @details
#' The print method displays basic information including:
#' - The type of the model object
#' - Number of iterations used in the estimation procedure
#' - Prompts the user to use \code{summary()} for more detailed output.
#'
#' @return Invisibly returns the input object \code{x}.
#' @export
#' @method print jmleIRT
print.jmleIRT <- function(x, ...) {
  cat("JMLE Rasch Model Object\n")
  cat(sprintf("Iterations: %d\n", x$iterations))
  cat("Use summary() for detailed output.\n")

  invisible(x)
}
