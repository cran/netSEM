test_that("netSEMp1 gives correct output", {
  expect_s3_class(netSEMp1(pet), c("netSEMp1", "list"))
})

test_that("netSEMp1 shows error if exogenous variable doesn't exist", {
  expect_error(netSEMp1(pet, exogenous = 100),
               "exogenous location out of range!")
  expect_error(netSEMp1(pet, exogenous = 0.5),
               "object 'exogenous.loc' not found")
  expect_error(netSEMp1(pet, exogenous = "test"),
               "exogenous 'test' does not exist!")
})

test_that("netSEMp1 shows error if endogenous variable doesn't exist", {
  expect_error(netSEMp1(pet, endogenous = 100),
               "endogenous location out of range!")
  expect_error(netSEMp1(pet, endogenous = 0.5),
               "replacement has length zero")
  expect_error(netSEMp1(pet, endogenous = "test"),
               "endogenous 'test' does not exist!")
})