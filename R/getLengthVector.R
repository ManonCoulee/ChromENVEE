#' Function to create a vector of limit value
#'
#' @title getLengthVector
#' @param lim limit value for distance
#'
#' @import stringr
#'
#' @return a vector of limit
getLengthVector <- function(lim) {

  lim <- lim/2
  limitVector <- seq(0, lim, length.out = 6)

  limitLabel <- unlist(lapply(seq_len(length(limitVector)), function(limValue) {

    label <- paste0(str_replace(as.character(as.integer(limitVector[limValue])), "000$", "kb"), "-",
      str_replace(as.character(as.integer(limitVector[limValue + 1])), "000$", "kb"))

    if((limValue + 1) > length(limitVector)) {
      label <- paste0(">",str_replace(as.character(limitVector[limValue]), "000$", "kb"))
    }
    return(label)
  }))

  names(limitLabel) <- limitVector
  return(limitLabel)
}
