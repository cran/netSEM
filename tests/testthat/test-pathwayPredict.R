test_that("Check if pathwayPredict output is list", {
    x <- netSEMp1(acrylic)
    paths <- pathwayRMSE(x, maxlen=3)
    expect_type(pathwayPredict(x, paths[10,2]), "list")
  })