#fichier test_getInformation.R
library(testthat)

test_that("test", {

  lim = 500000
  name = seq(0,250000,length.out = 6)
  output = c("0-50kb","50kb-100kb","100kb-150kb","150kb-200kb","200kb-250kb",">250kb")

  names(output) = name

  ## Test for the output conformity
  expect_identical(getLengthVector(lim),output)
})
