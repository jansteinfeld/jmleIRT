# jml_estimation.R
# Rasch model (1PL) estimation via JML (C++) with optional WLE; NA tolerated.
# Public R API aligned with C++ implementations in estimate_jmle.cpp

#' Fit Rasch model by JML with optional WLE
#'
#' Fits the Rasch (1PL) model using joint maximum likelihood via an efficient
#' C++ Newton-Raphson implementation. Missing responses (NA) are allowed.
#' Optional epsilon adjustment helps avoid biased estimates due to extreme-score degeneracy, and
#' centering ensures identifiability. Optionally returns Warm's WLE person
#' estimates computed with the estimated item difficulties.
#'
#' @param X Numeric matrix of responses (0/1) with possible NA for missing.
#' @param max_iter Maximum NR iterations in C++ (default 1000).
#' @param conv Convergence threshold on maximum absolute parameter change.
#' @param eps Epsilon adjustment for extreme person scores (default 0).
#' @param bias_correction Logical, simple post-hoc scaling of item betas.
#' @param center either "items" or "persons"; enforce mean(beta)=0 or mean(theta)=0.
#' @param max_update Step-size clip per update to stabilize NR (default 1.5).
#' @param verbose Logical, periodic progress from C++ (every 10 iters).
#' @param estimatewle Logical, compute Warm's WLE per person using estimated beta.
#' @param wle_adj Small adjustment to avoid extreme raw scores in WLE (default 1e-8).
#'
#' @return A list with components:
#'   - theta: numeric vector of person parameters
#'   - beta: numeric vector of item difficulties (mean 0 if center = TRUE)
#'   - iterations: integer iterations used
#'   - converged: logical convergence indicator
#'   - bias_correction, centered: echoes of inputs
#'   - wle_estimate: numeric vector of WLE (if estimatewle = TRUE) or NA
#' @export
jmle_estimation <- function(
  X,
  max_iter = 1000,
  conv = 1e-6,
  eps = 0,
  bias_correction = FALSE,
  center = "items",
  max_update = 1.5,
  verbose = FALSE,
  estimatewle = FALSE,
  wle_adj = 1e-8
) {
  # Basic validation and coercion
  X <- .coerce_binary_matrix(X)
  stopifnot(is.numeric(max_iter), length(max_iter) == 1L)
  stopifnot(is.numeric(conv), length(conv) == 1L)
  stopifnot(is.numeric(eps), length(eps) == 1L)
  stopifnot(is.logical(bias_correction), length(bias_correction) == 1L)
  stopifnot(is.character(center), length(center) == 1L)
  stopifnot(is.numeric(max_update), length(max_update) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)
  stopifnot(is.logical(estimatewle), length(estimatewle) == 1L)

  if (ncol(X) == 0L) {
    stop("X must have at least one item (column).")
  }
  if (nrow(X) == 0L) {
    stop("X must have at least one person (row).")
  }

  # Delegate to compiled routine (Rcpp)
  res <- estimate_jmle(
    X_ = X,
    max_iter = as.integer(max_iter),
    conv = conv,
    eps = eps,
    bias_correction = bias_correction,
    center = center,
    max_update = max_update,
    verbose = verbose,
    estimatewle = estimatewle,
    wle_adj = wle_adj
  )

  # Attach a class for downstream methods
  class(res) <- c("jmleIRT", class(res))
  res
}

#' Compute Warm's WLE for given item difficulties
#'
#' Convenience wrapper to compute person WLE given a response matrix and
#' item difficulties, using the C++ implementation. NA in X are allowed.
#'
#' @param X Numeric response matrix (0/1) with possible NA.
#' @param beta Numeric vector of item difficulties, length ncol(X).
#' @param max_iter Max iterations for person-level NR (default 100).
#' @param tol Convergence tolerance on NR step size (default 1e-10).
#' @param lower_ext,upper_ext Optional finite bounds for extreme raw scores.
#'
#' @return List with raw_score, wle, standard_error, conv, iterations.
#' @export
wle_estimation <- function(
  X,
  beta,
  max_iter = 100,
  tol = 1e-10,
  lower_ext = NA_real_,
  upper_ext = NA_real_
) {
  X <- .coerce_binary_matrix(X)
  beta <- .coerce_numeric_vector(beta, n = ncol(X), name = "beta")
  res <- estimate_wle(
    X_ = X,
    beta = beta,
    max_iter = as.integer(max_iter),
    tol = tol,
    lower_ext = lower_ext,
    upper_ext = upper_ext
  )
  class(res) <- c("wle_estimation", class(res))
  res
}

#' High-level Rasch model estimation
#'
#' User-friendly front-end that currently dispatches to JML; kept for
#' backward-compatibility and future method expansion.
#'
#' @param X Numeric response matrix (0/1, NA allowed).
#' @param method Character, currently only "jml" supported here.
#' @param ... Passed to jmle_estimation().
#' @return List as from jmle_estimation().
#' @export
estimate_rasch_model <- function(X, method = "jml", ...) {
  method <- match.arg(method, choices = c("jml"))
  switch(method,
         jml = jmle_estimation(X, ...))
}

# Internal helpers ----------------------------------------------------------

# Coerce to numeric 0/1 matrix with NA tolerated; accept logical/integer.
.coerce_binary_matrix <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) stop("X must be a matrix.")
  # Require a numeric-like matrix (numeric, integer, or logical)
  if (!is.numeric(X) && !is.integer(X) && !is.logical(X)) {
    stop("X must be a numeric matrix (0/1 with NA allowed).")
  }
  storage.mode(X) <- "double"
  # Allow values in {0,1,NA}
  bad <- !(is.na(X) | X == 0 | X == 1)
  if (any(bad, na.rm = TRUE)) {
    stop("X must contain only 0, 1, or NA.")
  }
  X
}

.coerce_numeric_vector <- function(x, n = NULL, name = "x") {
  if (!is.numeric(x)) stop(sprintf("%s must be numeric.", name))
  if (!is.null(n) && length(x) != n) {
    stop(sprintf("Length of %s must be %d.", name, n))
  }
  as.numeric(x)
}

