context("test-models.R")

test_that("models pbtk1cpt_model", {

  pbtk1cpt_model()
  mName <- "pbtk1cpt"
  #mcsim_pkg()
  #pksensi:::compile_model_pkg(mName)
  #q <- c("qunif", "qunif", "qunif", "qnorm")
  #q.arg = list(list(min = 0.5, max = 1),
  #             list(min = 0.5, max = 1.5),
  #             list(min = 0.02, max = 0.3),
  #             list(mean = 70, sd = 5))

  #set.seed(1234)
  #params <- c("vdist", "ke", "kgutabs", "BW")
  #x <- rfast99(params, n = 200, q = q, q.arg = q.arg, replicate = 10)
  #t <- seq(from = 0.01, to = 24.01, by = 1)
  #conditions <- c("Agutlument = 10")
  #outputs <- c("Ccompartment")
  #out <- solve_mcsim(x, mName = mName, params = params,
  #                   time = t, vars = outputs, condition = conditions)

  expect_true(file.exists("pbtk1cpt.model"))
  #expect_true(ifelse(.Platform$OS.type == "unix",
  #                   file.exists("mcsim.pbtk1cpt"),
  #                   file.exists("mcsim.pbtk1cpt.exe")))
})

test_that("models pbtk1cpt_model", {

  pbpk_apap_model()
  expect_true(file.exists("pbpk_apap.model"))

})
