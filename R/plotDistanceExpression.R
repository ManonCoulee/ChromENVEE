#' Function to create a plot which represent the gene expression in function distance gene-enhancer
#'
#' @title plotDistanceExpression
#' @param dataTable GRanges object or list of GRanges output of enhancerExpression
#' @param colorTable dataframe which contains color information (colorTable data)
#' @param xlab x-axis label (default = "")
#' @param ylab y-axis label (default = "log(CPM)")
#' @param limit limit for distance analysis (default = 500000 (500kb))
#'
#' @import ggplot2
#' @importFrom stats aggregate
#' @importFrom methods is
#'
#' @examples
#' listTableEnhancer = system.file("extdata", listTableEnhancer, package = "ChromENVEE")
#' data(listTableEnhancer)
#' genomeFile = system.file("extdata", genomeFile, package = "ChromENVEE")
#' data(genomeFile)
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' colorTable = system.file("extdata", "colorTable", package = "ChromENVEE")
#' data(colorTable)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' expression = enhancerExpression(anno, geneExpression)
#' plotDistanceExpression(expression, colorTable = colorTable)
#'
#' @return gene expression in function distance gene-enhancer
#' @export
plotDistanceExpression = function(dataTable,
					xlab = "",
          ylab = "log(CPM)",
          colorTable,
          limit = 500000) {

  col = getStateColor(colorTable = colorTable)

  lim = getLengthVector(limit)

  information_table = getInformation(dataTable)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in names(lim)) {
    pos = which(information_table$distance > as.numeric(l))
    information_table[pos,"limit"] = lim[l]
  }

  if(is(dataTable, "list")) {
    ylab = "mean(CPM)"
    list_mean = lapply(lim, function(x) {
      tt = information_table[information_table$limit == x,]
      df = aggregate(tt$expression,list(tt$sample_name), mean)
      df$grp = x
      colnames(df) = c("sample_name","expression","grp")
      return(df)
    })
    mean = do.call(rbind,list_mean)

    p = ggplot(mean,aes(x = factor(grp, levels = lim), y = expression,
      color = sample_name,
      group = sample_name, linetype = sample_name)) +
      geom_point(size = 3) +
      geom_line(size = 2) +
      scale_color_manual(values = col$stateName) +
      scale_linetype_manual(values = rep(c("solid","dashed"),length(unique(mean$sample_name))/2)) +
      xlab(xlab) + ylab(ylab) +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "none")
  } else {
    p = ggplot(information_table,aes(x = factor(limit, levels = lim), y = log(expression+0.01))) +
      geom_boxplot(aes(fill = chromatin_state),width = 0.1) +
      scale_fill_manual(values = col$stateNumber) +
      xlab(xlab) + ylab(ylab) +
      themePlot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(),
        legend.position = "none")
  }

  return(p)
}
