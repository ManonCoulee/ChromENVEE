#' Function to create a plot which represent the distribution of distance gene-enhancer
#'
#' @title plotGeneDistance
#' @param dataTable GRanges object or list of GRanges output of enhancerExpression
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
plotGeneDistance = function(dataTable,
					limit = 500000,
					xlab = "",
					ylab = "distance enhancer-gene (bp)") {

  lim = getLengthVector(limit)

  information_table = getInformation(dataTable)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in names(lim)) {
    pos = which(information_table$distance > as.numeric(l))
    information_table[pos,"distance_red"] = lim[l]
  }

  if(is(dataTable, "list")) {
    p = ggplot(information_table,aes(x = sample_name, fill = factor(distance_red,levels = lim))) +
      geom_bar(stat = "count", position = "fill") +
      coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
      theme_void() +
      theme(axis.text.y = element_text(),
        legend.position = "bottom") +
      scale_fill_manual(values = c("#f1948a","#c39bd3","#85c1e9","#76d7c4","#f7dc6f","#f0b27a"))
  } else {
    p = ggplot(information_table,aes(x = "", fill = factor(distance_red,levels = lim))) +
      geom_bar(stat = "count", position = "fill") +
      coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
      theme_void() +
      theme(axis.text.y = element_text(),
        legend.position = "bottom") +
      scale_fill_manual(values = c("#f1948a","#c39bd3","#85c1e9","#76d7c4","#f7dc6f","#f0b27a"))
  }
  return(p)
}
