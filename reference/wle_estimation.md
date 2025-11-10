# Compute Warm's WLE for given item difficulties

Convenience wrapper to compute person WLE given a response matrix and
item difficulties, using the C++ implementation. NA in X are allowed.

## Usage

``` r
wle_estimation(
  X,
  beta,
  max_iter = 100,
  tol = 1e-10,
  lower_ext = NA_real_,
  upper_ext = NA_real_
)
```

## Arguments

- X:

  Numeric response matrix (0/1) with possible NA.

- beta:

  Numeric vector of item difficulties, length ncol(X).

- max_iter:

  Max iterations for person-level NR (default 100).

- tol:

  Convergence tolerance on NR step size (default 1e-10).

- lower_ext, upper_ext:

  Optional finite bounds for extreme raw scores.

## Value

List with raw_score, wle, standard_error, conv, iterations.
