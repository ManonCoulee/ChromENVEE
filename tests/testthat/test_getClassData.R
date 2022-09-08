#fichier test_getClassData.R
library(testthat)

test_that("Test the conformity of initial data", {
  dataTable = c("table","expresion")

  expect_error(getInformation(dataTable),
    regexp = "'dataTable' must be a GRanges object or a list of GRanges object")

  expect_error(plotChromatinState(tableChromatinState = dataTable),
    regexp = "'tableChromatinState' must be a data frame or a list of data frame")

})
