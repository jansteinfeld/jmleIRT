set.seed(123)

make_X <- function(N = 10, I = 6, na_frac = 0.1) {
  X <- matrix(rbinom(N * I, 1, 0.5), nrow = N)
  if (na_frac > 0) {
    nNA <- ceiling(length(X) * na_frac)
    X[sample(length(X), nNA)] <- NA
  }
  X
}


# --- Basic functionality ---
test_that("jmle_estimation runs and returns expected structure", {
  X <- make_X(20, 8, na_frac = 0.15)
  fit <- jmle_estimation(
    X,
    max_iter = 200, conv = 1e-5, eps = 0,
    bias_correction = "none", center = "items",
    max_update = 1.5, verbose = FALSE, estimatewle = FALSE
  )
  expect_type(fit$theta, "double")
  expect_type(fit$beta, "double")
  expect_equal(length(fit$theta), nrow(X))
  expect_equal(length(fit$beta), ncol(X))
  expect_true(abs(mean(fit$beta)) < 1e-6) # centered
  # expect_type(fit$converged, "logical")
  expect_type(fit$iterations, "integer")
  expect_type(fit$bias_correction, "logical")
  expect_type(fit$center, "character")
  expect_true(all(is.na(fit$wle_estimate))) # not requested
})

# --- WLE estimation ---
test_that("jmle_estimation computes WLE when requested", {
  X <- make_X(30, 6, na_frac = 0.05)
  fit <- jmle_estimation(X, estimatewle = TRUE, max_iter = 300, center = "items")
  expect_equal(length(fit$theta), nrow(X))
  expect_equal(length(fit$wle_estimate), nrow(X))
  expect_false(all(is.na(fit$wle_estimate)))
})

# --- Bias correction ---
test_that("jmle_estimation applies bias correction option", {
  X <- make_X(25, 7)
  fit1 <- jmle_estimation(X, bias_correction = "none", center = "items")
  fit2 <- jmle_estimation(X, bias_correction = "simple", center = "items")
  fit3 <- jmle_estimation(X, bias_correction = "analytic", center = "items")
  expect_true("theta_analytic" %in% names(fit3))
  expect_true("beta_analytic" %in% names(fit3))
  expect_equal(length(fit1$beta), length(fit2$beta))
  expect_false(isTRUE(all.equal(fit1$beta, fit2$beta)))
})

# --- Handling missing data ---
test_that("jmle_estimation tolerates missing responses", {
  X <- make_X(20, 6, na_frac = 0.2)
  fit <- jmle_estimation(X, center = "items")
  expect_equal(length(fit$theta), nrow(X))
  expect_equal(length(fit$beta), ncol(X))
})

# --- Error handling ---
test_that("jmle_estimation errors on invalid input", {
  expect_error(jmle_estimation(matrix(2, 5, 5)),
    regexp = "must contain only 0, 1, or NA"
  )
  expect_error(jmle_estimation(matrix("a", 5, 5)),
    regexp = "numeric matrix"
  )
})


test_that("jmle_estimation applies analytic bias correction", {
  X <- make_X(20, 6)
  fit <- jmle_estimation(X, bias_correction = "analytic", center = "items")
  expect_true("theta_analytic" %in% names(fit))
  expect_true("beta_analytic" %in% names(fit))
  expect_equal(length(fit$theta_analytic), nrow(X))
  expect_equal(length(fit$beta_analytic), ncol(X))
})
