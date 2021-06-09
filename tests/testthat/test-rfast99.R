context("test-rfast99.R")

test_that("rfast99 single_timepoint_no_replication", {

  FFPK <- pksensi::FFPK

  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1),
               list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))

  set.seed(1234)
  x<-rfast99(params=c("F","KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 1, conf = 0.95)
  time <- 0.5
  out<-solve_fun(x, model = FFPK, time = time, vars = "output")
  check(out)
  print(out)

  expect_equal(dim(out$y)[2], 1)
  expect_error(pksim(out), "The time point must greater than 1", fix = TRUE)
})


test_that("rfast99 single_timepoint_20_replicaiton", {

  FFPK <- pksensi::FFPK

  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1),
               list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))

  set.seed(1234)
  x<-rfast99(params=c("F","KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 20, conf = 0.95)
  time <- 0.5
  out<-solve_fun(x, model = FFPK, time = time, vars = "output")
  check(out)
  print(out)

  expect_equal(dimnames(out$y)[[3]], "0.5")
  expect_equal(dim(out$y)[3], 1)
})

test_that("rfast99 class", {

  FFPK <- pksensi::FFPK

  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1),
               list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))

  set.seed(1234)
  x<-rfast99(params=c("F","KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 20, conf = 0.95)

  time <- seq(from = 0.25, to = 12.25, by = 0.5)
  out<-solve_fun(x, model = FFPK, time = time, vars = "output")
  check(out)
  print(out)

  #heat_check(out, text = TRUE)
  #heat_check(out, level = F)
  #cheat_check(out, index = "CI")

  expect_silent(pksim(out))
  expect_silent(pksim(out, log = TRUE))
  expect_silent(plot(out))
  expect_that(x, is_a("rfast99"))
  expect_equal(x$M, 4)
  expect_equal(x$replicate, 20)
  expect_equal(x$conf, 0.95)
  expect_true(class(x$a) == "array")
  expect_true(class(x$omega) == "numeric")
})

test_that("rfast99 error", {
  expect_error(rfast99(params = 20, n = 100, q.arg = list(min = 0, max = 1)),
               "Please assign the distribution(s) of quantile function", fix = TRUE)
  expect_error(rfast99(params = 20, n = 100, q = "qunif"),
               "Please assign the arguments for defined distribution", fix = TRUE)
})

test_that("rfast99 factor", {
  params <- c("A","B","C")
  output <- rfast99(params = params, n = 100, q = "qunif", q.arg = list(min = 0, max = 1))
  expect_true(class(output$params) == "character")
})

test_that("rfast99 M", {
  output <- rfast99(params = 3, n = 100, M = 1, q = "qunif", q.arg = list(min = 0, max = 1))
  expect_equal(output$M, 1)
})

test_that("rfast99 rep1", {

  FFPK <- pksensi::FFPK

  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1),
               list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))

  set.seed(1234)
  x<-rfast99(params=c("F","KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 20, conf = 0.95)

  time <- 0.5
  out<-solve_fun(x, model = FFPK, time = time, vars = "output")
  check(out)
  print(out)

  expect_equal(x$replicate, 20)
})

test_that("rfast99 rep2", {

  FFPK <- pksensi::FFPK

  q = "qunif"
  q.arg = list(list(min = 0.5, max = 1),
               list(min = 0.5, max = 1.5),
               list(min = 0.02, max = 0.3),
               list(min = 20, max = 60))

  set.seed(1234)
  x<-rfast99(params=c("F","KA","KE","V"), n = 100, q = q, q.arg = q.arg,
             rep = 20, conf = 0.95)

  time <- c(0.5, 1)
  out<-solve_fun(x, model = FFPK, time = time, vars = "output", lnparam = TRUE)
  check(out, time = 1)
  print(out)

  expect_equal(x$replicate, 20)
})

