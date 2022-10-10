#' Function to return chromatin state color table
#'
#' @param colorTable a data frame which contains color information
#'
#' @return list of vector with each chromatin state is associated to color
#' @export
getStateColor = function(colorTable) {

  lapply(c("stateName","stateNumber","colorValue"), function(x){
    if(!(x %in% colnames(colorTable))) {
      stop("colorTable colnames must be 'stateName','stateNumber' and 'colorValue'")
    }
  })

  colorName = colorTable$colorValue
  names(colorName) = colorTable$stateName

  colorNumber = colorTable$colorValue
  names(colorNumber) = colorTable$stateNumber

  col = list("stateName" = colorName, "stateNumber" = colorNumber)

  return(col)
}
