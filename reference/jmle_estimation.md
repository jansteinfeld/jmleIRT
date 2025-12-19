# Fit Rasch model by JML with optional WLE

Fits the Rasch (1PL) model using joint maximum likelihood via an
efficient C++ Newton-Raphson implementation. Missing responses (NA) are
allowed. Optional epsilon adjustment helps avoid biased estimates due to
extreme-score degeneracy, and centering ensures identifiability.
Optionally returns Warm's WLE person estimates computed with the
estimated item difficulties.

## Usage

``` r
jmle_estimation(
  X,
  max_iter = 1000,
  conv = 1e-06,
  eps = 0,
  bias_correction = c("none", "simple", "analytic"),
  center = "items",
  max_update = 1.5,
  verbose = FALSE,
  estimatewle = FALSE,
  wle_adj = 1e-08
)
```

## Arguments

- X:

  Numeric matrix of responses (0/1) with possible NA for missing.

- max_iter:

  Maximum NR iterations in C++ (default 1000).

- conv:

  Convergence threshold on maximum absolute parameter change.

- eps:

  Small positive value used to adjust extreme person scores. When
  `eps = 0` (default), persons with all 0 or all 1 responses are
  excluded from the iterative JMLE updates and assigned -Inf / +Inf
  afterward. When `eps > 0`, no deletion is performed and extreme scores
  are regularized internally.

- bias_correction:

  Character, one of `"none"`, `"simple"`, or `"analytic"`. `"none"`: no
  bias correction; `"simple"`: fast finite-sample scaling of item
  difficulties by factor (I - 1) / I; `"analytic"`: post-hoc first-order
  analytical bias correction using
  [`biasCorrection()`](jansteinfeld.github.io/jmleIRT/reference/biasCorrection.md).

- center:

  either "items" or "persons"; enforce mean(beta)=0 or mean(theta)=0.

- max_update:

  Step-size clip per update to stabilize NR (default 1.5).

- verbose:

  Logical, periodic progress from C++ (every 10 iters).

- estimatewle:

  Logical, compute Warm's WLE per person using estimated beta.

- wle_adj:

  Small adjustment to avoid extreme raw scores in WLE (default 1e-8).

## Value

A list with components: - theta: numeric vector of person parameters -
beta: numeric vector of item difficulties (mean 0 if center = "items") -
iterations: integer iterations used - converged: logical convergence
indicator - bias_correction_mode: selected bias correction mode -
center: centering used ("items" or "persons") - wle_estimate: numeric
vector of WLE (if estimatewle = TRUE) or NA - se_theta, se_beta:
approximate JMLE-based standard errors (if implemented) -
theta_analytic, beta_analytic: analytically bias-corrected parameters if
`bias_correction = "analytic"`, otherwise `NULL`
