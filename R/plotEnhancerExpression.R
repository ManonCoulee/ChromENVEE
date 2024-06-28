#' Function to create a plot which represent the distribution of gene expression
#'
#' @title plotEnhancerExpression
#' @param dataTable GRanges object or list of GRanges output of enhancerExpression
#' @param colorTable dataframe with color information
#' @param scale value (log10, log2 or none) to rescale expression (default = "none")
#' @param distance numeric value to focus on region (default = 0)
#' @param xlab x-axis label (default = "")
#' @param ylab y-axis label (default = "gene expression (CPM)")
#'
#' @import ggplot2
#'
#' @examples
#' data(listTableEnhancer)
#' data(genomeFile)
#' data(geneExpression)
#' data(colorTable)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' expression = enhancerExpression(anno, geneExpression)
#' plotEnhancerExpression(expression, colorTable = colorTable)
#'
#' @return ggplot2 figure corresponding to the distribution of gene expression
#' @export
plotEnhancerExpression = function(enhancerTable,
  xlab = "",
  ylab = "gene expression (CPM)",
  scale = "none",
  colorTable,
  distance = 0) {

  if(!is(distance, "numeric")) {
    stop("'distance' must be a numeric object")
  }

  if (!(scale %in% c("none","log10","log2"))) {
    stop("'scale' must be 'none', 'log10','log2' value. Default is 'none'")
  }

  informationTable <- getInformation(enhancerTable)
  informationTable$expression <- as.numeric(informationTable$expression)
  informationTable$distance <- as.numeric(informationTable$distance)

  if(distance != 0) {
    informationTable <- informationTable[informationTable$distance <= distance,]
  }

  col <- getStateColor(colorTable = colorTable)

  if(scale == "log10") {
    informationTable$expression <- log10(informationTable$expression + 0.01)
  } else if(scale == "log2") {
    informationTable$expression <- log2(informationTable$expression + 0.01)
  }

  if(is(enhancerTable, "list")) {
      plot = ggplot(informationTable, aes(x = sampleName, y = expression)) +
        geom_violin(aes(fill = chromatinState), color = "black") +
        geom_boxplot(aes(x = sampleName, fill = chromatinState), width = 0.1) +
        scale_fill_manual(values = col$stateNumber) +
        xlab(xlab) + ylab(ylab) +
        themePlot() +
        theme(axis.text.x = element_text(),
          axis.text.y = element_text(),
          legend.position = "none")
  } else {
      plot = ggplot(informationTable, aes(x = "", y = expression)) +
        geom_violin(aes(fill = chromatinState),color = "black") +
        geom_boxplot(width = 0.1) +
        scale_fill_manual(values = col$stateNumber) +
        xlab(xlab) + ylab(ylab) +
        themePlot() +
        theme(axis.text.x = element_text(angle = 45),
          axis.text.y = element_text(),
          legend.position = "none")
  }
  return(plot)
}
