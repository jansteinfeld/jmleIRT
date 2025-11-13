# First steps with the Package jmleIRT

## Introduction

The `jmleIRT` package provides an implementation of Joint Maximum
Likelihood Estimation (JMLE) for the Rasch (1PL) model. It allows
estimation of item difficulties and person abilities from dichotomous
response data with bias correction and weighted likelihood estimation
options.

## Rasch Model Specification

For person $p = 1,\ldots,P$ and item $i = 1,\ldots,I$, the Rasch model
posits the probability of a correct response $X_{pi}$ as:

$$P\left( X_{pi} = 1 \mid \theta_{p},\beta_{i} \right) = \frac{\exp\left( \theta_{p} - \beta_{i} \right)}{1 + \exp\left( \theta_{p} - \beta_{i} \right)}$$

where $\theta_{p}$ represents the ability of person $p$ and $\beta_{i}$
the difficulty of item $i$.

## Joint Maximum Likelihood Estimation

The joint likelihood across all persons and items is

$$L(\theta,\beta) = \prod\limits_{p = 1}^{P}\prod\limits_{i = 1}^{I}P\left( X_{pi} \mid \theta_{p},\beta_{i} \right)$$

JMLE alternates updates for $\theta_{p}$ and $\beta_{i}$ until
convergence, providing simultaneous estimates of person and item
parameters.

## Bias Correction

Finite samples and extreme scores induce bias in JMLE estimates. The
package implements:

- **Epsilon bias stabilization:** Adding a small $\varepsilon > 0$ to
  avoid infinite logits and improve numerical stability.
- **Analytical bias correction:** Drawing on methods from psychometric
  literature to reduce bias in item difficulty estimates.

## Using the jmleIRT Package

Key features: - Missing values by design are allowed (`NA`). - Optional
epsilon adjustment (`eps`) to reduce estimation bias from extreme
scores. - Centering based on either “items” or “persons” to ensure
identifiability. - Post-hoc bias correction
(`bias_correction = TRUE/FALSE`). - Optional computation of Warm’s
weighted likelihood estimates (WLE) for person parameters
(`estimatewle = TRUE`). - Fast Newton-Raphson optimization in C++.

### Simulating Data

We first simulate a dataset of person responses.

``` r
N <- 100  # persons
I <- 10   # items
X <- matrix(rbinom(N * I, 1, 0.5), nrow = N)
# randomly set 10% of entries to missing
X[sample(length(X), size = 0.1 * length(X))] <- NA
head(X)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    1    0    1    0    0    0    0    1     1
#> [2,]    1    1    1    0    1    0    1    1    0     0
#> [3,]    1    0    0    1    1    0    0    0    0     1
#> [4,]    0    0    0    1    0    1    1    1    1     1
#> [5,]   NA    0   NA    0    0    0   NA    1    0     0
#> [6,]    0    1    1    0    0    1    1    0    0     0
```

### Basic JML Estimation

Run the JML estimation with centering on “items” (default):

``` r
fit <- jmle_estimation(
  X, max_iter = 500, conv = 1e-5,
  center = "items", bias_correction = FALSE,
  estimatewle = FALSE, verbose = FALSE
)

str(fit)
#> List of 8
#>  $ data           : num [1:100, 1:10] 0 1 1 0 NA 0 0 0 1 NA ...
#>  $ theta          : num [1:100] -0.415 0.415 -0.415 0.415 -1.81 ...
#>  $ beta           : num [1:10] 0.472 -0.495 -0.216 -0.269 0.124 ...
#>  $ iterations     : int 5
#>  $ converged      : chr "TRUE"
#>  $ bias_correction: logi FALSE
#>  $ center         : chr "items"
#>  $ wle_estimate   : num [1:100] NA NA NA NA NA NA NA NA NA NA ...
#>  - attr(*, "class")= chr [1:2] "jmleIRT" "list"
```

The output is a list with:

- `theta`: person ability estimates
- `beta`: item difficulty estimates
- `iterations`: number of NR iterations
- `converged`: convergence flag
- `bias_correction`: whether applied
- `center`: ‘items’ or ‘persons’  
- `wle_estimate`: person WLEs (if requested)

### Bias Correction

Bias correction rescales item difficulty estimates.

``` r
fit_bc <- jmle_estimation(X, center = "items", bias_correction = TRUE)
summary(fit_bc$beta)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -0.445178 -0.230402  0.009101  0.000000  0.234625  0.425021
```

Compare to the uncorrected version:

``` r
summary(fit$beta)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -0.49464 -0.25600  0.01011  0.00000  0.26069  0.47225
```

### Warm’s WLE Estimates

Warm’s weighted likelihood estimates can be obtained in addition to JML
estimates.

``` r
fit_wle <- jmle_estimation(X, center = "items", estimatewle = TRUE)
head(fit_wle$wle_estimate)
#> [1] -0.3779385  0.3781464 -0.3779385  0.3781464 -1.4876787 -0.3779385
```

### Handling Missing Data

The function handles missing responses (`NA`) without failure:

``` r
X_miss <- X
X_miss[1:5, 1:2] <- NA
fit_miss <- jmle_estimation(X_miss, center = "items")

length(fit_miss$theta)
#> [1] 100
length(fit_miss$beta)
#> [1] 10
```

### Conclusion

The **jmleIRT** package provides a simple and efficient workflow for
Rasch JML estimation with useful options such as centering, bias
correction, and WLEs, while being robust to missing data.
