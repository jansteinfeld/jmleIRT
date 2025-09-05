
# jmleIRT: Rasch Model Joint Maximum Likelihood Estimation in R <img src="man/figures/jmleIRT.png" width="120" align="right" alt=""/>

[![CRAN
status](https://www.r-pkg.org/badges/version/jmleIRT)](https://cran.r-project.org/package=jmleIRT)
[![R-CMD-check](https://github.com/jansteinfeld/jmleIRT/actions/workflows/check-full.yaml/badge.svg)](https://github.com/jansteinfeld/jmleIRT/actions/workflows/check-full.yaml)

[![GitHub
version](https://img.shields.io/github/r-package/v/jansteinfeld/jmleIRT?label=version&logo=github)](https://github.com/jansteinfeld/jmleIRT/)
[![GitHub
release](https://img.shields.io/github/v/release/jansteinfeld/jmleIRT?label=release&logo=github)](https://github.com/jansteinfeld/jmleIRT/)
[![GitHub
issues](https://img.shields.io/github/issues-raw/jansteinfeld/jmleIRT?label=issues&logo=github)](https://github.com/jansteinfeld/jmleIRT/issues)

[![codecov](https://codecov.io/gh/jansteinfeld/jmleIRT/branch/master/graph/badge.svg?token=11lw4stBoI)](https://app.codecov.io/gh/jansteinfeld/jmleIRT)
[![CRAN
version](https://img.shields.io/cran/v/jmleIRT?label=CRAN%20version)](https://cran.r-project.org/package=jmleIRT)
[![CRAN
checks](https://badges.cranchecks.info/summary/jmleIRT.svg)](https://cran.r-project.org/web/checks/check_results_jmleIRT.html)
[![downloads](http://cranlogs.r-pkg.org/badges/last-month/jmleIRT?color=blue)](https://cran.r-project.org/package=jmleIRT)
[![License](https://img.shields.io/cran/l/jmleIRT)](https://opensource.org/license/GPL-3.0)

## Overview

The `jmleIRT` package provides tools to estimate Rasch (1PL) model
parameters using **Joint Maximum Likelihood Estimation (JMLE)** in R.
The package supports bias correction and Warm’s Weighted Likelihood
Estimates suitable for psychometric research.

------------------------------------------------------------------------

## Features

- Joint Maximum Likelihood Estimation of person abilities $\theta_p$ and
  item difficulties $\beta_i$.
- Bias correction to reduce estimation bias, including epsilon bias
  stabilization.
- Warm’s Weighted Likelihood Estimates (WLE) for person abilities.
- Support for missing data (missing responses coded as `NA`).
- Verbose mode to trace iteration progress.
- Comprehensive unit tests for reliability.

------------------------------------------------------------------------

## Installation

Install from CRAN:

``` r
install.packages("jmleIRT")
```

Or from GitHub development version:

``` r
# install.packages("remotes")
remotes::install_github("jansteinfeld/jmleIRT")
```

------------------------------------------------------------------------

## Usage Example

``` r
library(jmleIRT)

# Simulate response  40 persons (p) × 5 items (i)
set.seed(123)
X <- matrix(rbinom(40 * 5, 1, 0.5), nrow = 40)

# Estimate Rasch model using JMLE with bias correction
fit <- jmle_estimation(X, max_iter = 200, conv = 1e-5,
                      center = "items", bias_correction = TRUE)

# Item difficulties (centered)
print(fit$beta)

# Compute Weighted Likelihood Estimates for persons
wle_res <- estimate_wle(X, fit$beta)

# Person ability estimates via WLE
print(wle_res$wle)
```

------------------------------------------------------------------------

## Mathematical Model

## The Rasch Model

The Rasch model is a fundamental Item Response Theory (IRT) model used
for dichotomous item data. The probability that person $p$ with ability
$\theta_p$ answers item $i$ with difficulty $\beta_i$ correctly is
modeled as:

$$P(X_{pi} = 1 \mid \theta_p, \beta_i) = \frac{e^{\theta_p - \beta_i}}{1 + e^{\theta_p - \beta_i}}$$

where $X_{pi} \in \{0,1\}$ indicates the correctness of response.

## Joint Maximum Likelihood Estimation (JMLE)

JMLE simultaneously estimates person abilities
$\theta = (\theta_1, \ldots, \theta_P)$ and item difficulties
$\beta = (\beta_1, \ldots, \beta_I)$ by maximizing the joint likelihood
of observed responses:

$$L(\theta, \beta) = \prod_{p=1}^P \prod_{i=1}^I P(X_{pi} \mid \theta_p, \beta_i)$$

This is typically solved via an iterative algorithm alternating updates
for $\theta_p$ and $\beta_i$, until convergence criteria are met.

Although easy to implement and computationally efficient, JMLE estimates
suffer from bias, especially in small samples or with extreme response
patterns.

Here is the updated version of the Bias Correction Methods section,
formatted for integration into your README or vignette:

------------------------------------------------------------------------

## Bias Correction Methods

The Joint Maximum Likelihood Estimation (JMLE) for the Rasch model is
known to produce biased parameter estimates, especially in finite
samples and with extreme response patterns. The `jmleIRT` package
implements two key bias correction techniques:

- **Epsilon adjustment:** Proposed by Bertoli Barsotti, Punzo, and
  Lando, this method adds a small positive constant $$\varepsilon > 0$$
  during the estimation process. This adjustment prevents infinite
  log-odds for perfect or zero scores, ensuring finite item difficulty
  and ability estimates. It also helps alleviate bias related to extreme
  response patterns by smoothing the estimation.
- **Analytical bias correction:** This approach involves modifying the
  item difficulty estimates based on theoretical bias expressions
  derived from asymptotic expansions or minimum divergence estimators.
  It aims to reduce the structural bias inherent in the JMLE estimates
  as described in psychometric literature (e.g., Bertoli Barsotti &
  Punzo, 2012; Lando & Bertoli Barsotti, 2014).

In combination, the epsilon adjustment promotes numerical stability and
improved bias behavior at the estimation level, while analytical bias
correction provides targeted, theoretically grounded refinement of item
parameter estimates to enhance overall accuracy and robustness of the
JMLE algorithm.

## Package Features

- Estimation of person and item parameters using bias-corrected JMLE
- optional ability estimation based on Warm’s Weighted Likelihood
  Estimates (WLE)
- Handling of missing data
- documentation and examples
- Open-source development on GitHub for transparency and collaboration
- Optimized C++ backend for performance

## References

- Bertoli Barsotti, L., & Punzo, A. (2012). Comparison of two bias
  reduction techniques for the Rasch model. *Electronic Journal of
  Applied Statistical Analysis*, 5(3), 360-366.
- Lando, T., & Bertoli Barsotti, L. (2014). A modified minimum
  divergence estimator: some preliminary results for the Rasch model.
  *Electronic Journal of Applied Statistical Analysis*, 7(1), 37-57.
- Wright, B. D., & Panchapakesan, N. (1969). A Procedure for Sample-Free
  Item Analysis. *Educational and Psychological Measurement*, 29(1),
  23-48.

------------------------------------------------------------------------

## License

GPL-3

------------------------------------------------------------------------

Contributions welcome via GitHub: \[jansteinfeld/jmleIRT\].
