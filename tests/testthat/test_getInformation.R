#fichier test_getInformation.R
library(testthat)
library(GenomicRanges)

test_that("test", {

  data_table = list("table","expresion")
  data_frame = data.frame(c(1:10),c(1:10))

  data(listTableEnhancer)
  df = lapply(listTableEnhancer,function(x){x[,1]})

  ## Test for error message
  expect_error(getInformation(dataTable = data_table),
    regexp = "'dataTable' must be a list of GRanges object")

  expect_error(getInformation(dataTable = data_frame),
    regexp = "'dataTable' must be a GRanges object or a list of GRanges object")

  expect_error(getInformation(dataTable = listTableEnhancer[[1]][,1]),
    regexp = "GRanges object need a 'sample_name' column")

  expect_error(getInformation(dataTable = df),
    regexp = "GRanges object need a 'sample_name' column")

})
