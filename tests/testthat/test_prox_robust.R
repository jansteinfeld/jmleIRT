test_that("prox_algorithm produces reasonable estimates", {
  set.seed(123)
  dat <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10, ncol = 10)

  res <- prox_algorithm(dat)

  # Check that output is a list with specific names
  expect_named(res, c("b", "theta", "iter", "sigma.i", "sigma.n"))

  # Item difficulties and abilities should be numeric
  expect_true(is.numeric(res$b))
  expect_true(is.numeric(res$theta))

  # Iterations should be positive integer
  expect_true(res$iter > 0)

  # Sigma estimates should be positive
  expect_true(all(res$sigma.i >= 0))
  expect_true(all(res$sigma.n >= 0))
})
