# Iterative PROX Algorithm for Rasch Model Parameter Estimation with Bias Adjustment

This function implements the PROX algorithm for Rasch model parameter
estimation, enhanced with iterative refinement and a bias adjustment for
item difficulties following Alexander Robitzsch's method (as implemented
in the sirt package). It handles missing data in the response matrix by
ignoring missing responses.

## Usage

``` r
prox_robust(responses, max_iter = 10, tol = 1e-06)
```

## Arguments

- responses:

  Numeric matrix of dichotomous item responses (persons in rows, items
  in columns). Missing values (NA) are allowed and ignored in
  calculations.

- max_iter:

  Maximum number of iterations for the refinement step. Default is 10.

- tol:

  Convergence tolerance for iteration updates. Default is 1e-6.

## Value

A list with components:

- person_abilities:

  Numeric vector of estimated person abilities.

- item_difficulties:

  Numeric vector of estimated item difficulties.

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
data(matrix_data) # binary response matrix
out <- prox_robust(response_matrix)
print(out$person_abilities)
print(out$item_difficulties)
} # }
```
