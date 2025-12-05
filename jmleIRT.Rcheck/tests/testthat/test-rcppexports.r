test_that("RcppExports wrappers callable", {
  X <- matrix(rbinom(20, 1, 0.5), nrow = 5)
  beta <- rep(0, ncol(X))
  out_wle <- estimate_wle(X, beta, max_iter = 10L)
  expect_true(is.list(out_wle))
  out_jml <- estimate_jmle(
    X, max_iter = 10L, conv = 1e-4, eps = 0, bias_correction = FALSE,
    center = "items", max_update = 1.0, verbose = FALSE, estimatewle = FALSE
  )
  expect_true(is.list(out_jml))
})
