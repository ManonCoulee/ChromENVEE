#' Function to return information (expression and distance) for each enhancer
#'
#' @title getInformation
#' @param enhancerTable GRanges object or list of GRanges output of enhancerExpression
#'
#' @import ggplot2
#' @importFrom methods is
#'
#' @return distance, expression, gene name information
#' @export
getInformation <- function(enhancerTable) {
  if(is(enhancerTable, "list")) {
    enhancerInformationTable <- lapply(enhancerTable, function(enhancerTableFile) {
      if(!is(enhancerTableFile, "GRanges")) {
        stop("'enhancerTable' must be a list of GRanges object")
      }
      if(length(enhancerTableFile$sampleName) == 0) {
        stop("GRanges object need a 'sampleName' column")
      }
      enhancerTableFile <- enhancerTableFile[enhancerTableFile$geneAssociation != 0,]
      geneName <- unlist(strsplit(unlist(enhancerTableFile$geneList), ";"))
      expression <- unlist(strsplit(unlist(enhancerTableFile$geneExpression), ";"))
      distance <- unlist(strsplit(unlist(enhancerTableFile$distance), ";"))

      # geneName <- enhancerTableFile$geneList ##RADA
      # expression <- enhancerTableFile$geneExpression ##RADA
      # distance <- enhancerTableFile$distance ##RADA
      chromatinState <- enhancerTableFile$chromatinState

      dataFrame <- data.frame(geneName, expression, distance)
      dataFrame <- dataFrame[dataFrame$geneName != "NA", ]
      dataFrame$sampleName <- unique(enhancerTableFile$sampleName)
      group <- unique(enhancerTableFile$chromatinState)
      dataFrame$chromatinState <- unlist(lapply(group, function(state) {
        count <- sum(enhancerTableFile[enhancerTableFile$chromatinState == state]$geneAssociation)
        # count <- length(enhancerTableFile[enhancerTableFile$chromatinState == state]) ##RADA
        rep(state, count)
      }))
      dataFrame <- dataFrame[dataFrame$expression != "NA", ]
      dataFrame$expression <- as.numeric(dataFrame$expression)
      dataFrame$distance <- as.numeric(dataFrame$distance)

      return(dataFrame)
    })

    enhancerInformationTableMerge <- do.call(rbind, enhancerInformationTable)
    return(enhancerInformationTableMerge)

  } else if(is(enhancerTable, "GRanges")) {
    if(length(enhancerTable$sampleName) == 0) {
      stop("GRanges object need a 'sampleName' column")
    }

    enhancerTable <- enhancerTable[enhancerTable$geneAssociation != 0,]
    geneName <- unlist(strsplit(unlist(enhancerTable$geneList), ";"))
    expression <- unlist(strsplit(unlist(enhancerTable$geneExpression), ";"))
    distance <- unlist(strsplit(unlist(enhancerTable$distance), ";"))

    # geneName <- enhancerTable$geneList ##RADA
    # expression <- enhancerTable$geneExpression ##RADA
    # distance <- enhancerTable$distance ##RADA

    dataFrame <- data.frame(geneName, expression, distance)
    group <- unique(enhancerTable$chromatinState)
    dataFrame$chromatinState <- unlist(lapply(group, function(state) {
      count <- sum(enhancerTable[enhancerTable$chromatinState == state]$geneAssociation)
      # count <- length(enhancerTable[enhancerTable$chromatinState == state]) ##RADA
      rep(state, count)
    }))
    dataFrame <- dataFrame[dataFrame$expression != "NA", ]
    dataFrame$distance <- as.numeric(dataFrame$distance)
    dataFrame$expression <- as.numeric(dataFrame$expression)

    return(dataFrame)
  } else {
    stop("'enhancerTable' must be a GRanges object or a list of GRanges object")
  }
}
