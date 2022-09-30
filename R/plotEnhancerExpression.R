#' Function to create a plot which represent the distribution of gene expression
#'
#' @param dataTable a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param color a list of color value
#' @param stateNumber a list of chromatin state number
#' @param stateName a list of chromatin state name
#' @param scale a value (log10, log2 or none) to rescale expression (default = "none")
#' @param distance a value to zoom the expression focus (default = 0)
#' @param xlab a string (default = "")
#' @param ylab a string (default = "gene expression (CPM)")
#'
#' @import ggplot2
#'
#' @return ggplot2 figure corresponding to the distribution of gene expression
#' @export
plotEnhancerExpression = function(dataTable,
  xlab = "",
  ylab = "gene expression (CPM)",
  scale = "none",
  color,
  stateNumber,
  stateName, distance = 0) {

  if(class(distance) != "numeric") {
    stop("'distance' must be a numeric object")
  }

  if (!(scale %in% c("none","log10","log2"))) {
    stop("'scale' must be 'none', 'log10','log2' value. Default is 'none'")
  }

  information_table = getInformation(dataTable)
  information_table$expression = as.numeric(information_table$expression)
  information_table$distance = as.numeric(information_table$distance)

  if(distance != 0) {
    information_table = information_table[information_table$distance <= distance,]
  }

  col = getStateColor(stateName = stateName, stateNumber = stateNumber, color = color)

  if(scale == "log10") {
    information_table$expression = log10(information_table$expression+0.01)
  } else if(scale == "log2") {
    information_table$expression = log2(information_table$expression+0.01)
  }

  if(class(dataTable) == "list") {
      p = ggplot(information_table,aes(x = sample_name, y = expression)) +
        geom_violin(aes(fill = chromatin_state),color = "black") +
        geom_boxplot(width = 0.1) +
        scale_fill_manual(values = col$stateNumber) +
        xlab(xlab) + ylab(ylab) +
        themePlot() +
        theme(axis.text.x = element_text(),
          axis.text.y = element_text(),
          legend.position = "none")
  } else {
      p = ggplot(information_table,aes(x = "", y = expression)) +
        geom_violin(aes(fill = chromatin_state),color = "black") +
        geom_boxplot(width = 0.1) +
        scale_fill_manual(values = col$stateNumber) +
        xlab(xlab) + ylab(ylab) +
        themePlot() +
        theme(axis.text.x = element_text(),
          axis.text.y = element_text(),
          legend.position = "none")
  }

  return(p)
}
