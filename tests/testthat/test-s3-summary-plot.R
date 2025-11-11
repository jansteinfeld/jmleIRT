library(testthat)

X <- matrix(c(
  1, 0, 1, NA,
  0, 1, 1, 1,
  1, 1, 0, 0,
  NA, NA, 0, 1,
  0, 0, 0, NA
), nrow = 5, byrow = TRUE)

test_that("jmle_estimation runs and returns expected structure", {
  res <- jmle_estimation(X, max_iter = 200, estimatewle = TRUE, verbose = FALSE, center = "items")
  expect_s3_class(res, "jmleIRT")
  expect_true(is.numeric(res$theta))
  expect_true(is.numeric(res$beta))
  expect_true(is.integer(res$iterations) || is.numeric(res$iterations))
  expect_true(is.logical(res$converged))
  expect_true(length(res$theta) == nrow(X))
  expect_true(length(res$beta) == ncol(X))
  expect_true(length(res$wle_estimate) == nrow(X))
})

test_that("jmle_estimation convergence flag behaves correctly", {
  res <- jmle_estimation(X, max_iter = 1, center = "items")
  expect_false(res$converged)
})

test_that("jmle_estimation throws errors on invalid input", {
  expect_error(jmle_estimation(matrix(nrow = 0, ncol = 5), center = "items"), "at least one person")
  expect_error(jmle_estimation(matrix(nrow = 5, ncol = 0), center = "items"), "at least one item")
  expect_error(jmle_estimation(as.matrix(data.frame(a = c(1, 2)))), "X must contain only 0, 1, or NA")
})

# test_that(" summary works and prints correct output", {
#   res <- jmle_estimation(X, center="items")
#   # assuming a summary.jmleIRT method exists
#   expect_output(summary(res), "Summary of JMLE Rasch Model Estimation")
# })

test_that("wle_estimation returns sensible output", {
  res_jmle <- jmle_estimation(X, center = "items")
  res_wle <- wle_estimation(X, beta = res_jmle$beta)
  expect_true(is.numeric(res_wle$wle))
  expect_true(length(res_wle$wle) == nrow(X))
  expect_true(all(names(res_wle) %in% c("raw_score", "wle", "standard_error", "conv", "iterations")))
})

test_that("center argument centers beta and shifts theta accordingly", {
  res_cent <- jmle_estimation(X, center = "items")
  expect_equal(round(mean(res_cent$beta), 6), 0)
  res_no_cent <- jmle_estimation(X, center = "none")
  expect_true(abs(mean(res_no_cent$beta)) > 1e-4)
})

test_that("bias_correction scales betas appropriately", {
  res_bc <- jmle_estimation(X, bias_correction = TRUE, center = "items")
  expect_true(!all(abs(res_bc$beta) == 0))
})

test_that("eps parameter affects initialization and avoids infinite logits", {
  X_extreme <- matrix(1, nrow = 2, ncol = 4)
  expect_silent(jmle_estimation(X_extreme, eps = 0.1, center = "items"))
})

test_that("max_update clips parameter update sizes", {
  expect_silent(jmle_estimation(X, max_update = 0.01, center = "items"))
})

test_that("verbose option prints progress messages", {
  expect_output(jmle_estimation(X, verbose = TRUE, center = "items"), "Iter")
})
