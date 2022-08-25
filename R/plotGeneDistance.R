#' Function to create a plot which represent the distribution of distance gene-enhancer
#'
#' @param data_table a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param xlab a string
#' @param ylab a string
#' @param limit a value of limit for distance analysis
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the distribution of distance gene-enhancer
#' @export
plotGeneDistance = function(data_table, limit = 500000,
  xlab = "",
  ylab = "distance enhancer-gene (bp)") {

  lim = getLengthVector(limit)

  information_table = getInformation(data_table)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in names(lim)) {
    pos = which(information_table$distance > as.numeric(l))
    information_table[pos,"distance_red"] = lim[l]
  }

  p = ggplot(information_table,aes(x = sample_name, fill = factor(distance_red,levels = lim))) +
    geom_bar(stat = "count", position = "fill") +
    coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
    themePlot() +
    theme(axis.text.y = element_text(),
      legend.position = "bottom") +
    scale_fill_manual(values = c("#f1948a","#c39bd3","#85c1e9","#76d7c4","#f7dc6f","#f0b27a"))

  return(p)
}
