#' Analytical First-Order Bias Correction for Joint Maximum Likelihood Estimators in the Rasch Model
#'
#' This function performs an analytical bias correction for the Joint Maximum Likelihood (JML)
#' estimators of person (theta) and item (beta) parameters in the Rasch model.
#'
#' The JML estimators are well-known to suffer from incidental parameter bias when the number of items (\eqn{I})
#' is finite, resulting in biased parameter estimates. This bias is generally of order \eqn{1/I}.
#'
#' The bias correction implemented here follows the first-order analytic correction approach,
#' which approximates the bias term \eqn{B} for each parameter as:
#' \deqn{
#' \hat{\theta}^{corr} = \hat{\theta} - \frac{1}{I} \hat{B}_{\theta}, \quad
#' \hat{\beta}^{corr} = \hat{\beta} - \frac{1}{I} \hat{B}_{\beta}
#' }
#' where \eqn{\hat{B}} is estimated based on the observed scores \eqn{u_i} and the expected Fisher information \eqn{v_i},
#' calculated as
#' \deqn{
#' \hat{B}_i \approx \frac{E[v_i u_i]}{E[v_i^2]}.
#' }
#' In this implementation, scores and information are calculated based on the logistic form of the Rasch model,
#' where \eqn{p_{ni} = \frac{e^{\theta_n - \beta_i}}{1 + e^{\theta_n - \beta_i}}} is the probability that person \eqn{n}
#' correctly answers item \eqn{i}.
#'
#' Missing responses (NA) in the data matrix \eqn{X} are handled by omitting those entries.
#'
#' @param theta Numeric vector of estimated person parameters \eqn{\hat{\theta}}.
#' @param beta Numeric vector of estimated item parameters \eqn{\hat{\beta}}.
#' @param X Numeric matrix of responses (persons by items), coded as 0/1, with NA allowed for missing responses.
#' @param I Integer scalar representing the number of items (used for bias scaling).
#'
#' @return A named list containing bias-corrected parameter vectors:
#' \describe{
#'   \item{theta}{Numeric vector of bias-corrected person parameters.}
#'   \item{beta}{Numeric vector of bias-corrected item parameters.}
#' }
#'
#' @details
#' The bias correction is inspired by the incidental parameters literature (Neyman & Scott, 1948;
#' Lancaster, 2000; Arellano & Hahn, 2006), which shows that Maximum Likelihood estimators
#' in models with many nuisance parameters (like person parameters in Rasch) are biased when
#' the number of observations per parameter is limited.
#'
#' This function applies a computationally efficient closed-form bias correction using the first and second derivatives
#' of the log-likelihood of the Rasch model likelihood function, evaluated at the JML estimates.
#'
#' The estimator reduces bias by estimating expected score and information terms:
#' \deqn{
#' u_{ni} = X_{ni} - p_{ni}, \quad v_{ni} = p_{ni} (1 - p_{ni})
#' }
#' for person \eqn{n} and item \eqn{i}, and then aggregating these across items or persons.
#'
#' The bias terms for person \eqn{n} and item \eqn{i} are estimated as:
#' \deqn{
#' \hat{B}_{\theta,n} = \frac{\sum_i v_{ni} u_{ni}}{\sum_i v_{ni}^2}, \quad
#' \hat{B}_{\beta,i} = \frac{\sum_n v_{ni} (p_{ni} - X_{ni})}{\sum_n v_{ni}^2}
#' }
#'
#' This method is applicable when the number of items \eqn{I} is fixed and moderate, and
#' the number of persons \eqn{N} is large.
#'
#' @references
#' - Neyman, J., & Scott, E. L. (1948). Consistent Estimates Based on Partially Consistent Observations. Econometrica, 16(1), 1-32.
#' - Lancaster, T. (2000). The incidental parameter problem since 1948. Journal of Econometrics, 95(2), 391-413.
#' - Arellano, M., & Hahn, J. (2006). Understanding Bias in Nonlinear Panel Models. Review of Economic Studies, 73(2), 557-578.
#'
#' @examples
#' \dontrun{
#' # Example data matrix with 2 persons and 4 items:
#' X <- matrix(c(1, 0, 1, NA, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
#' theta <- c(0.5, -0.5)
#' beta <- c(-0.2, 0.1, 0.3, -0.1)
#' I <- ncol(X)
#' corrected <- biasCorrectionJMLE(theta, beta, X, I)
#' print(corrected$theta)
#' print(corrected$beta)
#' }
#' @export
biasCorrection <- function(theta, beta, X, I) {
  stopifnot(is.numeric(theta))
  stopifnot(is.numeric(beta))
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(I))

  N <- nrow(X)
  I <- ncol(X)

  res <- biasCorrectionJMLE(theta = theta, beta = beta, X = X, N = N, I = I)
  # Attach a class for downstream methods
  class(res) <- c("biasCorrection", class(res))
  res
}
biasCorrection.jmleIRT <- function(jmle_obj) {
  if (missing(jmle_obj) || !inherits(jmle_obj, "jmleIRT")) {
    stop("You must provide a valid 'jmleIRT' object.")
  }

  # Extract necessary components from jmle_obj
  theta <- jmle_obj$theta
  beta <- jmle_obj$beta
  X <- jmle_obj$data
  I <- ncol(X)

  N <- nrow(X)
  I <- ncol(X)

  # Call your internal C++ bias correction function
  res <- biasCorrectionJMLE(theta = theta, beta = beta, X = X, N = N, I = I)

  # Attach class for downstream methods
  class(res) <- c("biasCorrection", class(res))

  return(res)
}
