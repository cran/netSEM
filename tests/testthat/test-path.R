test_that("path() gives correct output", {
  netSEM_ans <- netSEMp1(acrylic)
  expect_type(path(netSEM_ans, from = "IrradTot", to = "IAD2"), "list")
  expect_match(path(netSEM_ans, from = "IrradTot", to = "IAD2")$model.print, "IAD2 =")
})

test_that("path() errors if input class is not netSEMp1", {
  expect_error(path(acrylic, from = "IrradTot", to = "IAD2"),
               "x is not of class 'netSEM'")
})

test_that("path() errors if argument input is from model types", {
  netSEM_ans <- netSEMp1(acrylic)
  expect_error(path(netSEM_ans, from = "SL", to = "IAD2"),
               "'from' variable is not found.")
  expect_error(path(netSEM_ans, from = "IrradTot", to = "SL"),
               "'to' variable is not found.")
})

test_that("path() shows correct output if path doesn't exist", {
  netSEM_ans <- netSEMp1(acrylic)
  expect_type(path(netSEM_ans, from = "IAD2", to = "IrradTot"),
              "logical")
})