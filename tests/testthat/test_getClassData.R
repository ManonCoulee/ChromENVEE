#fichier test_getClassData.R
library(testthat)

test_that("Test the conformity of initial data", {
  data_table = c("table","expresion")

  expect_error(getInformation(data_table),
    regexp = "'data_table' must be a GRanges object or a list of GRanges object")

  expect_error(plotChromatinState(genomic_region = data_table),
    regexp = "'genomic_region' must be a data frame or a list of data frame")

})
