#fichier test_getClassData.R
library(testthat)

test_that("Test the conformity of initial data", {
  data_table = c("table","expresion")

  stateName = c("TSS","TSSFlnk","TSSFlnkD")
  stateNumber = c("U1","U2")
  color = c("#B71C1C","#E65100","#E65100","#43A047")

  data_frame = data.frame(c(1:10),c(1:10))

  ## Test for getStateColor
  expect_error(getStateColor(stateName,stateNumber,color),
    regexp = "Not same length between each parameters !")

  ## Test for getInformation
  expect_error(getInformation(dataTable = data_table),
    regexp = "'dataTable' must be a GRanges object or a list of GRanges object")
  expect_error(getInformation(dataTable = data_frame),
    regexp = "'dataTable' must be a GRanges object or a list of GRanges object")

  ## Test for plotChromatinState
  expect_error(plotChromatinState(tableChromatinState = data_table),
    regexp = "'tableChromatinState' must be a data frame or a list of data frame")

  ## Test for enhancerAnnotation
  expect_error(enhancerAnnotation(enhancerTable = data_frame),
    regexp = "'enhancerTable' is not a GRanges table")

  ## Test for plotEnhancerExpression
  expect_error(plotEnhancerExpression(scale = "exp", dataTable = data_frame),
    regexp = "'scale' must be 'none', 'log10','log2' value. Default is 'none'")

  expect_error(plotEnhancerExpression(distance = "0"),
    regexp = "'distance' must be a numeric object")
})
