#fichier test_stateColor.R
library(testthat)

test_that("test", {
  stateName = c("TSS","TSSFlnk","TSSFlnkD")
  stateNumber = c("U1","U2")
  color = c("#B71C1C","#E65100","#E65100","#43A047")

  expect_error(getStateColor(stateName,stateNumber,color),
    regexp = "Not same length between each parameters !")
})
