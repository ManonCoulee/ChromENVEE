#' Function to create a plot which represent the distribution of distance gene-enhancer
#'
#' @param data_table a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param xlab a string
#' @param ylab a string
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the distribution of distance gene-enhancer
#' @export
plotGeneDistance = function(data_table, xlab = "", ylab = "distance enhancer-gene (bp)") {

  limit = c(0,50000,100000,150000,200000,250000)
  limit_label = c("0-50kb","50-100kb","100-150kb","150-200kb","200-250kb",">250kb")

  information_table = getInformation(data_table)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in limit) {
    pos = which(information_table$distance > l)
    information_table[pos,"distance_red"] = limit_label[limit == l]
  }

  p = ggplot(information_table,aes(x = sample_name, fill = factor(distance_red,levels = limit_label))) +
    geom_bar(stat = "count", position = "fill") +
    coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
    themePlot() +
    theme(axis.text.y = element_text(),
      legend.position = "bottom")

    # theme_bw() + theme(strip.background  = element_blank(),
    #   text = element_text(size=35, angle = 0),
    #   panel.grid.major = element_line(colour = "grey80"),
    #   panel.border = element_blank(),
    #   axis.ticks = element_blank(),
    #   axis.text.x= element_blank(),
    #   panel.grid.minor.x=element_blank(),
    #   panel.grid.major.x=element_blank(),
    #   legend.position = "bottom") +
    scale_fill_manual(values = c("#f1948a","#c39bd3","#85c1e9","#76d7c4","#f7dc6f","#f0b27a"))

  return(p)
}
