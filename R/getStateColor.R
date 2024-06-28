#' Function to return chromatin state color table
#'
#' @title getStateColor
#' @param colorTable dataframe color information
#'
#' @examples
#' colorTable = system.file("extdata", colorTable, package = "ChromENVEE")
#' data(colorTable)
#' colorTable
#'
#' @return list of vector with each chromatin state is associated to color
#' @export
getStateColor <- function(colorTable) {

  lapply(c("stateName","stateNumber","colorValue"), function(colName) {
    if(!(colName %in% colnames(colorTable))) {
      stop("colorTable colnames must be 'stateName','stateNumber' and 'colorValue'")
    }
  })

  colorName <- colorTable$colorValue
  names(colorName) <- colorTable$stateName

  colorNumber <- colorTable$colorValue
  names(colorNumber) <- colorTable$stateNumber

  col <- list("stateName" = colorName, "stateNumber" = colorNumber)

  return(col)
}
