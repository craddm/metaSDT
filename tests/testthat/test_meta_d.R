
context("Check meta-d measures")

test_that("meta_d MLE gives expected result", {
  nR_S1 <- c(32, 23, 14, 9)
  nR_S2 <- c(6, 13, 26, 33)
  expect_equal_to_reference(fit_meta_d_MLE(nR_S1, nR_S2),
                            file = "met_mleConst.rds")

  expect_equal_to_reference(fit_meta_d_MLE(nR_S1, nR_S2,
                                           add_constant = FALSE),
                            file = "met_mleNoConst.rds")

})

test_that("meta_d_SSE gives expected results", {

  nR_S1 <- c(29, 33, 12, 7)
  nR_S2 <- c(3, 15, 21, 43)

  expect_equal_to_reference(fit_meta_d_SSE(nR_S1, nR_S2),
                            file = "met_sseConst.rds")

  expect_equal_to_reference(fit_meta_d_SSE(nR_S1, nR_S2,
                                           add_constant = FALSE),
                            file = "met_sseNoConst.rds")

})
