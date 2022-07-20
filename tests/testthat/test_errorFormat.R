#fichier test_getClassData.R
library(testthat)

test_that("Test the conformity of initial data", {
  data_table = c("table","expresion")

  state_name = c("TSS","TSSFlnk","TSSFlnkD")
  state_number = c("U1","U2")
  color_value = c("#B71C1C","#E65100","#E65100","#43A047")

  data_frame = data.frame(c(1:10),c(1:10))

  ## Test for stateColor
  expect_error(stateColor(state_name,state_number,color_value),
    regexp = "Not same length between each parameters !")

  ## Test for getInformation
  expect_error(getInformation(data_table = data_table),
    regexp = "'data_table' must be a GRanges object or a list of GRanges object")
  expect_error(getInformation(data_table = data_frame),
    regexp = "'data_table' must be a GRanges object or a list of GRanges object")

  ## Test for plotChromatinState
  expect_error(plotChromatinState(genomic_region = data_table),
    regexp = "'genomic_region' must be a data frame or a list of data frame")

  ## Test for enhancerAnnotation
  expect_error(enhancerAnnotation(enhancer_table = data_frame),
    regexp = "'enhancer_table' is not a GRanges table")

  ## Test for plotEnhancerExpression
  expect_error(plotEnhancerExpression(distance = "a",data_table = data_frame),
    regexp = "'distance' must be a numeric object")
  expect_error(plotEnhancerExpression(scale = "exp" ,data_table = data_frame),
    regexp = "'scale' must be 'none', 'log10','log2' value. Default is 'none'")

})
