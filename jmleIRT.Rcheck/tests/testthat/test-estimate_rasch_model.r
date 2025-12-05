test_that("estimate_rasch_model dispatches to JML", {
  X <- matrix(rbinom(30, 1, 0.4), nrow = 6)
  fit <- estimate_rasch_model(X, method = "jml", max_iter = 80)
  expect_true(is.list(fit))
  expect_true(all(c("theta","beta","iterations","converged","wle_estimate") %in% names(fit)))
})