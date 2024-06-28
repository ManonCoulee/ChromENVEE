#' Function to create a plot which represent the distribution of distance gene-enhancer
#'
#' @title plotGeneDistance
#' @param enhancerTable GRanges object or list of GRanges output of enhancerExpression
#' @param xlab x-axis label (default = "")
#' @param ylab y-axis label (default = "distance enhancer-gene (bp))
#' @param limit limit for distance analysis (default = 500000 (500kb))
#'
#' @import ggplot2
#'
#' @examples
#' data(listTableEnhancer)
#' data(genomeFile)
#' data(geneExpression)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' expression = enhancerExpression(anno, geneExpression)
#' plotGeneDistance(expression)
#'
#' @return distribution of distance gene-enhancer
#' @export
plotGeneDistance <- function(enhancerTable,
					limit = 500000,
					xlab = "",
					ylab = "distance enhancer-gene (bp)") {

  lim <- getLengthVector(limit)

  informationTable = getInformation(enhancerTable)
  informationTable$distance = as.numeric(informationTable$distance)
  informationTable$expression = as.numeric(informationTable$expression)

  for(limit in names(lim)) {
    pos = which(informationTable$distance > as.numeric(limit))
    informationTable[pos,"distanceRed"] = lim[limit]
  }

  if(is(enhancerTable, "list")) {
    plot = ggplot(informationTable,aes(x = sampleName, fill = factor(distanceRed, levels = lim))) +
      geom_bar(stat = "count", position = "fill") +
      coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
      theme_void() +
      theme(axis.text.y = element_text(),
        legend.position = "bottom") +
      scale_fill_manual(values = c("#f1948a", "#c39bd3", "#85c1e9", "#76d7c4", "#f7dc6f", "#f0b27a"))
  } else {
    plot = ggplot(informationTable,aes(x = "", fill = factor(distanceRed, levels = lim))) +
      geom_bar(stat = "count", position = "fill") +
      coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
      theme_void() +
      theme(axis.text.y = element_text(),
        legend.position = "bottom") +
      scale_fill_manual(values = c("#f1948a", "#c39bd3", "#85c1e9", "#76d7c4", "#f7dc6f", "#f0b27a"))
  }
  return(plot)
}
