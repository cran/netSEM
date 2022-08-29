test_that("summary.netSEMp1 is showing correct output", {
  p1_ans <- netSEMp1(acrylic)
  expect_length(summary(p1_ans), nrow(p1_ans$table))
})

test_that("summary.netSEMp2 is showing correct output", {
  p2_ans <- netSEMp2(acrylic)
  expect_length(summary(p2_ans), nrow(p2_ans$res.print))
})