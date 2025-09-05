test_that(".coerce_binary_matrix accepts matrix/data.frame and 0/1/NA", {
  Xm <- matrix(c(0,1,NA,1,0,0), nrow = 2)
  res <- jmleIRT:::`.coerce_binary_matrix`(Xm)
  expect_true(is.matrix(res))
  expect_type(res, "double")

  Xdf <- data.frame(a = c(0,1), b = c(NA,0))
  res2 <- jmleIRT:::`.coerce_binary_matrix`(Xdf)
  expect_true(is.matrix(res2))
  expect_equal(storage.mode(res2), "double")
})

test_that(".coerce_binary_matrix errors on non-matrix-like", {
  expect_error(jmleIRT:::`.coerce_binary_matrix`(list(1,2)), "must be a matrix") 
})

test_that(".coerce_binary_matrix errors on values outside {0,1,NA}", {
  Xbad <- matrix(c(-1,2,0,1), nrow = 2)
  expect_error(jmleIRT:::`.coerce_binary_matrix`(Xbad), "contain only 0, 1, or NA")
})

test_that(".coerce_numeric_vector enforces numeric and length", {
  expect_error(jmleIRT:::`.coerce_numeric_vector`(c("a","b")), "must be numeric")
  expect_error(jmleIRT:::`.coerce_numeric_vector`(1:2, n = 3, name = "beta"),
               "Length of beta must be 3")
  b <- jmleIRT:::`.coerce_numeric_vector`(1:3, n = 3, name = "beta")
  expect_equal(as.numeric(1:3), b)
})
