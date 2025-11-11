# Print method for JMLE Rasch Model Objects

This function prints a brief summary of a Joint Maximum Likelihood
Estimation (JMLE) Rasch model object of class `jmleIRT`.

## Usage

``` r
# S3 method for class 'jmleIRT'
print(x, ...)
```

## Arguments

- x:

  An object of class `jmleIRT` returned from a JMLE Rasch model
  estimation function.

- ...:

  Additional arguments passed to or from other methods (currently
  ignored).

## Value

Invisibly returns the input object `x`.

## Details

The print method displays basic information including: - The type of the
model object - Number of iterations used in the estimation procedure -
Prompts the user to use
[`summary()`](https://rdrr.io/r/base/summary.html) for more detailed
output.
