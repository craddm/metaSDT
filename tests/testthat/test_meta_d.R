library(metaSDT)

context("Check meta-d MLE")

test_that("meta_d MLE gives expected result", {
  nR_S1 <- c(32, 23, 14, 9)
  nR_S2 <- c(6, 13, 26, 33)
  expect_equal_to_reference(fit_meta_d_MLE(nR_S1, nR_S2), file = "met_mleConst.rds")
  expect_equal_to_reference(fit_meta_d_MLE(nR_S1, nR_S2, add_constant = FALSE), file = "met_mleNoConst.rds")
})

