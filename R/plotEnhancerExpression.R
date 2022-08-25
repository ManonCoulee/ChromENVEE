#' Function to create a plot which represent the distribution of gene expression
#'
#' @param data_table a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param color a list of color value
#' @param state_number a list of chromatin state number
#' @param state_name a list of chromatin state name
#' @param scale a value (log10, log2 or none) to rescale expression (default = "none")
#' @param distance a value to zoom the expression focus (default = 0)
#' @param xlab a string
#' @param ylab a string
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the distribution of gene expression
#' @export
plotEnhancerExpression = function(data_table,
  xlab = "",
  ylab = "gene expression log10(CPM)",
  scale = "none",
  color,
  state_number,
  state_name, distance = 0) {

  if(class(distance) != "numeric") {
    stop("'distance' must be a numeric object")
  }

  if (!(scale %in% c("none","log10","log2"))) {
    stop("'scale' must be 'none', 'log10','log2' value. Default is 'none'")
  }

  information_table = getInformation(data_table)
  information_table$expression = as.numeric(information_table$expression)
  information_table$distance = as.numeric(information_table$distance)

  if(distance != 0) {
    information_table = information_table[information_table$distance <= distance,]
  }

  col = getStateColor(state_name = state_name, state_number = state_number, color = color)

  if(scale == "log10") {
    information_table$expression = log10(information_table$expression+0.01)
  } else if(scale == "log2") {
    information_table$expression = log2(information_table$expression+0.01)
  }

  p = ggplot(information_table,aes(x = sample_name, y = expression)) +
    geom_violin(aes(fill = chromatin_state),color = "black") +
    geom_boxplot(width = 0.1) +
    scale_fill_manual(values = col$state_number) +
    xlab(xlab) + ylab(ylab) +
    themePlot() +
    theme(axis.text.x = element_text(),
      axis.text.y = element_text(),
      legend.position = "none")

  return(p)
}
