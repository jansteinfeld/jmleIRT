# High-level Rasch model estimation

User-friendly front-end that currently dispatches to JML; kept for
backward-compatibility and future method expansion.

## Usage

``` r
estimate_rasch_model(X, method = "jml", ...)
```

## Arguments

- X:

  Numeric response matrix (0/1, NA allowed).

- method:

  Character, currently only "jml" supported here.

- ...:

  Passed to jmle_estimation().

## Value

List as from jmle_estimation().
