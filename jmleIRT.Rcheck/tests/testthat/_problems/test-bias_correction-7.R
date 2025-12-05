# Extracted from test-bias_correction.r:7

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "jmleIRT", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
X <- matrix(c(1, 0, 1, NA, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
theta <- c(0.5, -0.5)
beta <- c(-0.2, 0.1, 0.3, -0.1)
I <- ncol(X)
result <- biasCorrectionJMLE(theta, beta, X, I)
