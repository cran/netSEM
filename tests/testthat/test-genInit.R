test_that("genInit function works", {
  df <-
    genInit(bounds = list(
      a1 = c(-4, -3),
      a2 = c(-1, 1),
      a3 = c(3, 4)
    ), k = 5)
  expect_equal(nrow(df), 5)
  expect_length(df, 3)
  expect_true(range(df)[1] >= -4)
})
