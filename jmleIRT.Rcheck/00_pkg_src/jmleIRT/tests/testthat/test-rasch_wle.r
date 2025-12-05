test_that("estimate_wle runs and returns fields", {
  set.seed(1)
  X <- matrix(rbinom(50, 1, 0.5), nrow = 10, ncol = 5)
  beta <- rep(0, ncol(X))
  res <- estimate_wle(X, beta, max_iter = 50, tol = 1e-8)
  expect_named(res, c("raw_score","wle","standard_error","conv","iterations"))
  expect_equal(length(res$raw_score), nrow(X))
  expect_equal(length(res$wle), nrow(X))
  expect_equal(length(res$standard_error), nrow(X))
})

test_that("estimate_wle handles NA-only rows and returns NA where appropriate", {
  X <- matrix(NA_real_, nrow = 3, ncol = 4)
  beta <- rep(0, 4)
  res <- estimate_wle(X, beta)
  expect_true(all(is.na(res$wle)))
  expect_true(all(res$conv == 0))
  expect_true(all(is.na(res$standard_error)))
}) 

test_that("estimate_wle handles extremes via x*-stabilization (wle_adj)", {
  # First two rows are extremes; third is mixed
  X <- rbind(rep(0, 4), rep(1, 4), c(1,0,1,0))
  beta <- c(-1, -0.5, 0.5, 1)

  # With Warm-style stabilization, estimates should be finite and SE computed
  res <- estimate_wle(X, beta, wle_adj = 1e-8)

  # raw scores
  expect_equal(res$raw_score[1], 0L)
  expect_equal(res$raw_score[2], 4L)
  expect_equal(res$raw_score[3], 2L)

  # Extremes: finite WLE due to fractional x*; should converge
  expect_true(is.finite(res$wle[1]))
  expect_true(is.finite(res$wle[2]))
  expect_true(res$conv[1] == 1)
  expect_true(res$conv[2] == 1)

  # SE for stabilized extremes should be finite and positive
  expect_true(is.finite(res$standard_error[1]) && res$standard_error[1] > 0)
  expect_true(is.finite(res$standard_error[2]) && res$standard_error[2] > 0)

  # Mixed row should also be finite with positive SE
  expect_true(is.finite(res$wle[3]))
  expect_true(is.finite(res$standard_error[3]) && res$standard_error[3] > 0)
}) 

test_that("estimate_wle respects input validation (matrix and beta length)", {
  X <- matrix(c(0,1,0,1), nrow = 2)
  expect_error(estimate_wle(list(1,2), beta = c(0,0)), "must be a matrix")
  expect_error(estimate_wle(X, beta = 1:3), "Length of 'beta' must equal number of columns")
}) 

test_that("estimate_jmle optionally returns WLE using final beta", {
  set.seed(123)
  X <- matrix(rbinom(80, 1, 0.5), nrow = 20, ncol = 4)

  # Run JMLE with WLE enabled; pass a small wle_adj for stability
  fit <- estimate_jmle(
    X,
    max_iter = 200,
    conv = 1e-6,
    eps = 0.0,
    bias_correction = FALSE,
    center = "items",
    max_update = 1.5,
    verbose = FALSE,
    estimatewle = TRUE,
    wle_adj = 1e-8
  )

  expect_true(all(c("theta","beta","iterations","converged","bias_correction","center","wle_estimate") %in% names(fit)))
  # wle_estimate present and length equals nrow(X)
  expect_equal(length(fit$wle_estimate), nrow(X))
  # beta length equals ncol(X)
  expect_equal(length(fit$beta), ncol(X))
}) 
