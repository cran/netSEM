test_that("Check if pathwayRMSE output is data frame", {
  x <- netSEMp1(acrylic)
  expect_s3_class(pathwayRMSE(x, maxlen=2), "data.frame")
})