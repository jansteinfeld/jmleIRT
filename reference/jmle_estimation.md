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
  bias_correction = FALSE,
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

  Epsilon adjustment for extreme person scores (default 0).

- bias_correction:

  Logical, simple post-hoc scaling of item betas.

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
beta: numeric vector of item difficulties (mean 0 if center = TRUE) -
iterations: integer iterations used - converged: logical convergence
indicator - bias_correction, centered: echoes of inputs - wle_estimate:
numeric vector of WLE (if estimatewle = TRUE) or NA
