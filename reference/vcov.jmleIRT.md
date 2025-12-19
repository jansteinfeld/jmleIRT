# Approximate covariance matrix for jmleIRT objects

Uses diagonal JMLE-based standard errors and returns a block-diagonal
covariance matrix (person block and item block).

## Usage

``` r
# S3 method for class 'jmleIRT'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `"jmleIRT"`.

- ...:

  Currently unused.

## Value

A covariance matrix with row/column names.
