#' Function to return chromatin state color table
#'
#' @param stateName a vector of chromatin state name
#' @param stateNumber a vector of chromatin state number
#' @param color a vector of color to associated to each chromatin state
#'
#' @return list of vector with each chromatin state is associated to color
#' @export
getStateColor = function(stateName, stateNumber, color) {

  ## Test si les 3 variable ont la même taille
  if(length(stateName) != length(stateNumber)) {
    stop("Not same length between each parameters !")
  }
  if(length(stateName) != length(color)) {
    stop("Not same length between each parameters !")
  }
  if(length(stateNumber) != length(color)) {
    stop("Not same length between each parameters !")
  }

  ## Creer le tableau
  colorName = color
  names(colorName) = stateName

  colorNumber = color
  names(colorNumber) = stateNumber

  col = list("stateName" = colorName, "stateNumber" = colorNumber)

  return(col)
}
