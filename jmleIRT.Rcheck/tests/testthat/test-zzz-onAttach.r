test_that(".onAttach emits the intended startup message (captured)", {
  res <- testthat::evaluate_promise(jmleIRT:::.onAttach(tempdir(), "jmleIRT"))
  expect_true(length(res$messages) >= 1)
  expect_true(any(grepl("jmleIRT", res$messages)))
  expect_true(any(grepl("Rasch \\(1PL\\) JML", res$messages)))
  expect_equal(res$warnings, character())
  expect_null(res$result)
})