#' @export print.jmleIRT
#' @export
print.jmleIRT <- function(x, ...) {
  cat("JMLE Rasch Model Object\n")
  cat(sprintf("Iterations: %d\n", x$iter))
  cat("Use summary() for detailed output.\n")
  invisible(x)
}