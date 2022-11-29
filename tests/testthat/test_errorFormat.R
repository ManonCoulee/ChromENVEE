#fichier test_getClassData.R
library(testthat)

test_that("Test the conformity of initial data", {
  data_table = c("table","expresion")

  data_frame = data.frame(c(1:10),c(1:10))

  ## Test for plotChromatinState
  expect_error(plotChromatinState(tableChromatinState = data_table),
    regexp = "'tableChromatinState' must be a dataframe or a list of dataframe")

  ## Test for enhancerAnnotation
  expect_error(enhancerAnnotation(enhancerTable = data_frame),
    regexp = "'enhancerTable' is not a GRanges table")

  ## Test for plotEnhancerExpression
  expect_error(plotEnhancerExpression(scale = "exp", dataTable = data_frame),
    regexp = "'scale' must be 'none', 'log10','log2' value. Default is 'none'")

  expect_error(plotEnhancerExpression(distance = "0"),
    regexp = "'distance' must be a numeric object")
})
