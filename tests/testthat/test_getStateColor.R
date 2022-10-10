#fichier test_getStateColor.R
library(testthat)

test_that("test", {

  df1 = data.frame(stateName = c("TSS","TSSFlnk"),
    stateNumber = c("U1","U2"),
    color = c("#B71C1C","#E65100"))

  df2 = data.frame(stateName = c("TSS","TSSFlnk"),
    stateNumber = c("U1","U2"),
    colorValue = c("#B71C1C","#E65100"))

  output = list("stateName" = c("TSS" = "#B71C1C","TSSFlnk" = "#E65100"),
    "stateNumber" = c("U1" = "#B71C1C","U2" = "#E65100"))

  ## Test for error message
  expect_error(getStateColor(df1),
    regexp = "colorTable colnames must be 'stateName','stateNumber' and 'colorValue'")

  ## Test for the output conformity
  expect_identical(getStateColor(df2),output)
})
