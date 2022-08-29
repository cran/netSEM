# plot.netSEMp1
test_that("Check if plot.netSEMp1 has correct output", {
  x <- netSEMp1(acrylic)
  expect_s3_class(plot.netSEMp1(x), c("grViz","htmlwidget"))
})

# plot.netSEMp2
test_that("Check if plot.netSEMp2 has correct output", {
  x <- netSEMp2(acrylic)
  expect_s3_class(plot.netSEMp2(x), c("grViz","htmlwidget"))
})