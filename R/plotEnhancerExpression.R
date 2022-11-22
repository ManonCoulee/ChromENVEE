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
#' listTableEnhancer = system.file("extdata", listTableEnhancer, package = "ChromENVEE")
#' data(listTableEnhancer)
#' genomeFile = system.file("extdata", genomeFile, package = "ChromENVEE")
#' data(genomeFile)
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' colorTable = system.file("extdata", "colorTable.rda", package = "ChromENVEE")
#' data(colorTable)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' expression = enhancerExpression(anno, geneExpression)
#' plotEnhancerExpression(expression, colorTable = colorTable)
#'
#' @return ggplot2 figure corresponding to the distribution of gene expression
#' @export
plotEnhancerExpression = function(dataTable,
  xlab = "",
  ylab = "gene expression (CPM)",
  scale = "none",
  colorTable, distance = 0) {

  if(!is(distance, "numeric")) {
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

  col = getStateColor(colorTable = colorTable)

  if(scale == "log10") {
    information_table$expression = log10(information_table$expression+0.01)
  } else if(scale == "log2") {
    information_table$expression = log2(information_table$expression+0.01)
  }

  if(is(dataTable, "list")) {
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
