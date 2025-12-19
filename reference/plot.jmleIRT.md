# Plot method for jmleIRT objects

Quick diagnostic plots for Rasch JMLE results.

## Usage

``` r
# S3 method for class 'jmleIRT'
plot(x, type = c("hist", "icc"), items = NULL, ...)
```

## Arguments

- x:

  An object of class `"jmleIRT"`.

- type:

  Character string indicating the type of plot: `"hist"` for histograms
  of person and item parameters, `"icc"` for item characteristic curves
  of selected items.

- items:

  Integer vector of item indices to plot when `type = "icc"`. Defaults
  to the first 5 items (or all if fewer than 5).

- ...:

  Currently unused.
