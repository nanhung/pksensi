context("test-rfast99.R")

test_that("rfast99 single_time", {
  FFPK <- function(parameters, times, dose = 320){
    A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
    CONC <- A*exp(- parameters[2] * times) - A*exp(- parameters[1] * times)
    return(CONC)
  }
  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))
  set.seed(1234)
  x<-rfast99(factors=c("KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 20, conf = 0.95)
  times <- 0.5
  y<-solve_fun(x, model = FFPK, times = times, output = "output")
  tell2(x,y)

  expect_equal(dimnames(y)[[3]], "0.5")
  expect_equal(dim(y)[3], 1)
})

test_that("rfast99 class", {
  FFPK <- function(parameters, times, dose = 320){
    A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
    CONC <- A*exp(- parameters[2] * times) - A*exp(- parameters[1] * times)
    return(CONC)
  }
  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))
  set.seed(1234)
  x<-rfast99(factors=c("KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 10, conf = 0.95)
  times <- seq(from = 0.25, to = 12.25, by = 0.5)
  y<-solve_fun(x, model = FFPK, times = times, output = "output")
  tell2(x,y)

  expect_silent(plot(x))
  expect_that(x, is_a("rfast99"))
  expect_equal(x$M, 4)
  expect_equal(x$replicate, 10)
  expect_equal(x$conf, 0.95)
  expect_true(class(x$a) == "array")
  expect_true(class(x$omega) == "numeric")
})

test_that("rfast99 error", {
  expect_error(rfast99(factors = 20, n = 100, q.arg = list(min = 0, max = 1)),
               "Please assign the distribution(s) of quantile function", fix = TRUE)
  expect_error(rfast99(factors = 20, n = 100, q = "qunif"),
               "Please assign the arguments for defined distribution", fix = TRUE)
})

test_that("rfast99 factor", {
  factors <- c("A","B","C")
  output <- rfast99(factors = factors, n = 100, q = "qunif", q.arg = list(min = 0, max = 1))
  expect_true(class(output$factors) == "character")
})

test_that("rfast99 M", {
  output <- rfast99(factors = 3, n = 100, M = 1, q = "qunif", q.arg = list(min = 0, max = 1))
  expect_equal(output$M, 1)
})

test_that("rfast99 rep1", {
  FFPK <- function(parameters, times, dose = 320){
    A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
    CONC <- A*exp(- parameters[2] * times) - A*exp(- parameters[1] * times)
    return(CONC)
  }
  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))
  set.seed(1234)
  x<-rfast99(factors=c("KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 1, conf = 0.95)
  times <- 0.5
  y<-solve_fun(x, model = FFPK, times = times, output = "output")
  tell2(x,y)

  expect_equal(x$replicate, 1)
})

test_that("rfast99 rep2", {
  FFPK <- function(parameters, times, dose = 320){
    A <- (dose * parameters[1])/( parameters[3]*( parameters[1]- parameters[2]))
    CONC <- A*exp(- parameters[2] * times) - A*exp(- parameters[1] * times)
    return(CONC)
  }
  q = "qunif"
  q.arg = list(list(min = log(0.5), max = log(1.5)),
               list(min = log(0.02), max = log(0.3)),
               list(min = log(20), max = log(60)))
  set.seed(1234)
  x<-rfast99(factors=c("KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 1, conf = 0.95)
  times <- c(0.5, 1)
  y<-solve_fun(x, model = FFPK, times = times, output = "output", lnparam = T)
  tell2(x,y)

  expect_equal(x$replicate, 1)
})

