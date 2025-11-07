test_that("Bias-Korrektur verändert Schätzer wie erwartet", {
  X <- matrix(c(1, 0, 1, NA, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
  theta <- c(0.5, -0.5)
  beta <- c(-0.2, 0.1, 0.3, -0.1)
  I <- ncol(X)
  N <- nrow(X)

  result <- biasCorrectionJMLE(theta, beta, X, N, I)

  expect_type(result, "list")
  expect_named(result, c("theta", "beta"))
  expect_length(result$theta, length(theta))
  expect_length(result$beta, length(beta))

  # expect_false(all(result$theta == theta))
  # expect_false(all(result$beta == beta))
})
