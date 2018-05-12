context("test-rfast99.R")

test_that("rfast99 class", {
  output <- rfast99(factors = 20, n = 100, q = "qunif", q.arg = list(min = 0, max = 1), replicate = 20, conf = 0.95)
  expect_that(output, is_a("rfast99"))
  expect_equal(output$M, 4)
  expect_equal(output$replicate, 20)
  expect_equal(output$conf, 0.95)
  expect_true(class(output$a) == "array")
  expect_true(class(output$omega) == "numeric")
})

test_that("rfast99 factor", {
  factors <- c("A","B","C")
  output <- rfast99(factors = factors, n = 100, q = "qunif", q.arg = list(min = 0, max = 1), replicate = 20, conf = 0.95)
  expect_true(class(output$factors) == "character")
})
