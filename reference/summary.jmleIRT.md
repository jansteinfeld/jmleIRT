# Summary of JML estimation

Print a brief summary for a JML fit, showing means of person and item
parameters. Expects an object with numeric components \`theta\` and
\`beta\`; designed for objects of class "jmleIRT".

## Usage

``` r
# S3 method for class 'jmleIRT'
summary(object, ...)
```

## Arguments

- object:

  An object of class "jmleIRT" containing numeric vectors \`theta\` and
  \`beta\`. If not a strict fit object, it must still provide \`theta\`
  and \`beta\` accessible as \`object\$theta\` and \`object\$beta\`.

- ...:

  Further arguments passed to or used by methods (currently unused).

- digits:

  Number of digits to display (default is 4).

## Value

The input object, returned invisibly.

\# Minimal example with an adapter object X \<- matrix(c( 1, 0, 1, NA,
0, 1, 1, 1, 1, 1, 0, 0, NA, NA, 0, 1, 0, 0, 0, NA ), nrow=5, byrow=TRUE)

res \<- jmle_estimation(X, max_iter=200, estimatewle=TRUE,
verbose=FALSE) summary(res)

## See also

\[jmle_estimation\]
