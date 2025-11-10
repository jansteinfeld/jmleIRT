# PROX Estimation for the Rasch Model

This function implements the PROX routine for Rasch models, ported from
Alexander Robitzsch's R implementation to C++ for speed.

## Usage

``` r
prox_algorithm(dat, dat.resp = NULL, freq = NULL, conv = 0.001, maxiter = 30)
```

## Arguments

- dat:

  Numeric matrix or data frame of binary item responses.

- dat.resp:

  Numeric matrix indicating missing data (1 = observed, 0 = missing). If
  missing, it is generated internally.

- freq:

  Numeric vector of frequencies for each response pattern. Defaults to 1
  for all respondents.

- conv:

  Numeric value specifying the convergence criterion for iteration.
  Defaults to 0.001.

- maxiter:

  Integer specifying the max number of iterations. Defaults to 30.

## Value

A list with elements:

- b:

  Estimated item difficulties (logits).

- theta:

  Estimated person abilities (logits).

- iter:

  Number of iterations performed.

- sigma.i:

  Estimated standard deviations of abilities for persons responding to
  items.

- sigma.n:

  Estimated standard deviations of difficulties for items responded to
  by persons.

## Details

The function starts with initial estimates using standardized raw scores
and applies a bias adjustment to item difficulties based on Robitzsch's
formula: \$\$ d_i = \mu_i - \sqrt{1 + \frac{\sigma_i^2}{2.9}} \times
\text{logit}(\text{item proportion}) \$\$ to correct bias in finite
samples. Then, person and item parameters are updated iteratively using
logistic transformations until convergence or the maximum number of
iterations is reached.

This implementation is primarily for quick and robust approximate Rasch
model estimation, useful as a starting point for joint maximum
likelihood estimation.

## References

Robitzsch, A. (2020). sirt: Supplementary Item Response Theory Models. R
package version 3.10-9. <https://CRAN.R-project.org/package=sirt>

Linacre, J.M. (1994). Many-Facet Rasch Measurement. Chicago: MESA Press.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10, ncol = 10)
result <- prox_algorithm(data)
print(result$b)
} # }
```
