#' Function to create a vector of limit value
#'
#' @title getLengthVector
#' @param lim limit value for distance
#'
#' @import stringr
#'
#' @return a vector of limit
getLengthVector = function(lim) {

  lim = lim/2
  limit = seq(0,lim,length.out = 6)

  limit_label = unlist(lapply(seq_len(length(limit)), function(l) {
    lab = paste0(str_replace(as.character(as.integer(limit[l])),"000$","kb"),"-",
      str_replace(as.character(as.integer(limit[l+1])),"000$","kb"))
    if((l+1) > length(limit)) {
      lab = paste0(">",str_replace(as.character(limit[l]),"000$","kb"))
    }
    return(lab)
  }))

  names(limit_label) =  limit
  return(limit_label)
}
