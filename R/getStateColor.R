#' Function to return chromatin state color table
#'
#' @param state_name a vector of chromatin state name
#' @param state_number a vector of chromatin state number
#' @param color a vector of color to associated to each chromatin state
#'
#' @return list of vector with each chromatin state associate to color
#' @export
getStateColor = function(state_name, state_number, color) {

  ## Test si les 3 variable ont la même taille
  if(length(state_name) != length(state_number)) {
    stop("Not same length between each parameters !")
  }
  if(length(state_name) != length(color)) {
    stop("Not same length between each parameters !")
  }
  if(length(state_number) != length(color)) {
    stop("Not same length between each parameters !")
  }

  ## Creer le tableau
  colorName = color
  names(colorName) = state_name

  colorNumber = color
  names(colorNumber) = state_number

  col = list("state_name" = colorName, "state_number" = colorNumber)

  return(col)
}
