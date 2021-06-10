context("test-models.R")

test_that("models pbtk1cpt_model", {

  pbtk1cpt_model()
  compile_model("pbtk1cpt")
  compile_model("pbtk1cpt", application = 'R')
  expect_true(file.exists("pbtk1cpt.model"))

})

test_that("models pbtk1cpt_model", {

  pbpk_apap_model()
  expect_true(file.exists("pbpk_apap.model"))

})
