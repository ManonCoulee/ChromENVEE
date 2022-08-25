#fichier test_stateColor.R
library(testthat)

test_that("test", {
  state_name = c("TSS","TSSFlnk","TSSFlnkD")
  state_number = c("U1","U2")
  color_value = c("#B71C1C","#E65100","#E65100","#43A047")

  expect_error(getStateColor(state_name,state_number,color_value),
    regexp = "Not same length between each parameters !")
})
