# Analytical First-Order Bias Correction for Joint Maximum Likelihood Estimators in the Rasch Model

This function performs an analytical bias correction for the Joint
Maximum Likelihood (JML) estimators of person (theta) and item (beta)
parameters in the Rasch model.

## Usage

``` r
biasCorrection(theta, beta, X, I)
```

## Arguments

- theta:

  Numeric vector of estimated person parameters \\\hat{\theta}\\.

- beta:

  Numeric vector of estimated item parameters \\\hat{\beta}\\.

- X:

  Numeric matrix of responses (persons by items), coded as 0/1, with NA
  allowed for missing responses.

- I:

  Integer scalar representing the number of items (used for bias
  scaling).

## Value

A named list containing bias-corrected parameter vectors:

- theta:

  Numeric vector of bias-corrected person parameters.

- beta:

  Numeric vector of bias-corrected item parameters.

## Details

The JML estimators are well-known to suffer from incidental parameter
bias when the number of items (\\I\\) is finite, resulting in biased
parameter estimates. This bias is generally of order \\1/I\\.

The bias correction implemented here follows the first-order analytic
correction approach, which approximates the bias term \\B\\ for each
parameter as: \$\$ \hat{\theta}^{corr} = \hat{\theta} - \frac{1}{I}
\hat{B}\_{\theta}, \quad \hat{\beta}^{corr} = \hat{\beta} - \frac{1}{I}
\hat{B}\_{\beta} \$\$ where \\\hat{B}\\ is estimated based on the
observed scores \\u_i\\ and the expected Fisher information \\v_i\\,
calculated as \$\$ \hat{B}\_i \approx \frac{E\[v_i u_i\]}{E\[v_i^2\]}.
\$\$ In this implementation, scores and information are calculated based
on the logistic form of the Rasch model, where \\p\_{ni} =
\frac{e^{\theta_n - \beta_i}}{1 + e^{\theta_n - \beta_i}}\\ is the
probability that person \\n\\ correctly answers item \\i\\.

Missing responses (NA) in the data matrix \\X\\ are handled by omitting
those entries.

The bias correction is inspired by the incidental parameters literature
(Neyman & Scott, 1948; Lancaster, 2000; Arellano & Hahn, 2006), which
shows that Maximum Likelihood estimators in models with many nuisance
parameters (like person parameters in Rasch) are biased when the number
of observations per parameter is limited.

This function applies a computationally efficient closed-form bias
correction using the first and second derivatives of the log-likelihood
of the Rasch model likelihood function, evaluated at the JML estimates.

The estimator reduces bias by estimating expected score and information
terms: \$\$ u\_{ni} = X\_{ni} - p\_{ni}, \quad v\_{ni} = p\_{ni} (1 -
p\_{ni}) \$\$ for person \\n\\ and item \\i\\, and then aggregating
these across items or persons.

The bias terms for person \\n\\ and item \\i\\ are estimated as: \$\$
\hat{B}\_{\theta,n} = \frac{\sum_i v\_{ni} u\_{ni}}{\sum_i v\_{ni}^2},
\quad \hat{B}\_{\beta,i} = \frac{\sum_n v\_{ni} (p\_{ni} -
X\_{ni})}{\sum_n v\_{ni}^2} \$\$

This method is applicable when the number of items \\I\\ is fixed and
moderate, and the number of persons \\N\\ is large.

## References

\- Neyman, J., & Scott, E. L. (1948). Consistent Estimates Based on
Partially Consistent Observations. Econometrica, 16(1), 1-32. -
Lancaster, T. (2000). The incidental parameter problem since 1948.
Journal of Econometrics, 95(2), 391-413. - Arellano, M., & Hahn, J.
(2006). Understanding Bias in Nonlinear Panel Models. Review of Economic
Studies, 73(2), 557-578.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example data matrix with 2 persons and 4 items:
X <- matrix(c(1, 0, 1, NA, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
theta <- c(0.5, -0.5)
beta <- c(-0.2, 0.1, 0.3, -0.1)
I <- ncol(X)
corrected <- biasCorrectionJMLE(theta, beta, X, I)
print(corrected$theta)
print(corrected$beta)
} # }
```
