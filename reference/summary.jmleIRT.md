# Summary of JMLE Rasch Model Estimation

Print a summary of the Joint Maximum Likelihood Estimation (JMLE) Rasch
model results, including number of iterations, item difficulties, and
person ability parameters.

## Usage

``` r
# S3 method for class 'jmleIRT'
summary(object, digits = 4, ...)
```

## Arguments

- object:

  An object of class `"jmleIRT"` containing numeric vectors `theta` and
  `beta`. If not a strict fit object, it must still provide `theta` and
  `beta` accessible as `object$theta` and `object$beta`.

- digits:

  Number of digits for numeric rounding (default 4).

- ...:

  Additional arguments passed to or from other methods (currently
  unused).

## Value

Returns the input `object` invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(c(1, 0, 1, NA, 0, 1, 1, 1, 1, 1, 0, 0, NA, NA, 0, 1, 0, 0, 0, NA), nrow = 5, byrow = TRUE)
res <- jmle_estimation(X, max_iter = 200, estimatewle = TRUE, verbose = FALSE)
summary(res)
} # }
```
