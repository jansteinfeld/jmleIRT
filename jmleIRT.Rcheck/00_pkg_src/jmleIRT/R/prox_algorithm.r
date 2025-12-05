#' PROX Estimation for the Rasch Model
#'
#' This function implements the PROX routine for Rasch models,
#' ported from Alexander Robitzsch's R implementation to C++ for speed.
#'
#' @param dat Numeric matrix or data frame of binary item responses.
#' @param dat.resp Numeric matrix indicating missing data (1 = observed, 0 = missing).
#' If missing, it is generated internally.
#' @param freq Numeric vector of frequencies for each response pattern.
#' Defaults to 1 for all respondents.
#' @param conv Numeric value specifying the convergence criterion for iteration.
#' Defaults to 0.001.
#' @param maxiter Integer specifying the max number of iterations.
#' Defaults to 30.
#' @return A list with elements:
#' \describe{
#'   \item{b}{Estimated item difficulties (logits).}
#'   \item{theta}{Estimated person abilities (logits).}
#'   \item{iter}{Number of iterations performed.}
#'   \item{sigma.i}{Estimated standard deviations of abilities for persons responding to items.}
#'   \item{sigma.n}{Estimated standard deviations of difficulties for items responded to by persons.}
#' }
#' @details
#' The function starts with initial estimates using standardized raw scores and applies
#' a bias adjustment to item difficulties based on Robitzsch's formula:
#' \deqn{
#' d_i = \mu_i - \sqrt{1 + \frac{\sigma_i^2}{2.9}} \times \text{logit}(\text{item proportion})
#' }
#' to correct bias in finite samples. Then, person and item parameters are updated iteratively
#' using logistic transformations until convergence or the maximum number of iterations is reached.
#'
#' This implementation is primarily for quick and robust approximate Rasch model estimation,
#' useful as a starting point for joint maximum likelihood estimation.
#'
#' @references
#' Robitzsch, A. (2020). sirt: Supplementary Item Response Theory Models. R package version 3.10-9.
#' \url{https://CRAN.R-project.org/package=sirt}
#'
#' Linacre, J.M. (1994). Many-Facet Rasch Measurement. Chicago: MESA Press.
#' @examples
#' \dontrun{
#' data <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10, ncol = 10)
#' result <- prox_algorithm(data)
#' print(result$b)
#' }
#' @export
prox_algorithm <- function(dat, dat.resp = NULL, freq = NULL,
                           conv = 0.001, maxiter = 30) {
  if (is.data.frame(dat)) dat <- as.matrix(dat)
  if (is.null(dat.resp)) {
    dat.resp <- ifelse(!is.na(dat), 1, 0)
  }
  if (is.null(freq)) {
    freq <- rep(1, nrow(dat))
  }
  prox_rasch(dat, dat.resp, freq, conv, maxiter)
}
