# Analytical bias correction for jmleIRT objects

Applies the C++ first-order bias correction to an object of class
`jmleIRT` and returns corrected person and item parameters.

## Usage

``` r
# S3 method for class 'jmleIRT'
biasCorrection(jmle_obj, theta = NULL, beta = NULL, X = NULL, I = NULL, ...)
```

## Arguments

- jmle_obj:

  An object of class `"jmleIRT"` as returned by

- theta:

  Person parameters (ignored if jmle_obj is provided)

- beta:

  Item parameters (ignored if jmle_obj is provided)

- X:

  Response matrix (ignored if jmle_obj is provided)

- I:

  Number of items (ignored if jmle_obj is provided)
  [`jmle_estimation()`](jansteinfeld.github.io/jmleIRT/reference/jmle_estimation.md).

- ...:

  further arguments, but currently unused.

## Value

An object of class `"biasCorrection"` containing numeric vectors `theta`
and `beta`.
