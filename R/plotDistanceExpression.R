#' Function to create a plot which represent the gene expression in function distance gene-enhancer
#'
#' @param dataTable a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param colorTable a data frame which contains color information
#' @param xlab a string (default = "")
#' @param ylab a string (default = "log(CPM)")
#' @param limit a value of limit for distance analysis (default = 500000 (500kb))
#'
#' @import ggplot2
#'
#' @return ggplot2 figure corresponding to the mean of gene expression in function distance gene-enhancer
#' @export
plotDistanceExpression = function(dataTable,
  xlab = "",
  ylab = "log(CPM)",
  colorTable, limit = 500000) {

  col = getStateColor(colorTable)

  lim = getLengthVector(limit)

  information_table = getInformation(dataTable)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in names(lim)) {
    pos = which(information_table$distance > as.numeric(l))
    information_table[pos,"limit"] = lim[l]
  }

  if(class(dataTable) == "list") {
    ylab = "mean(CPM)"
    list_mean = lapply(lim, function(x) {
      tt = information_table[information_table$limit == x,]
      df = aggregate(tt$expression,list(tt$sample_name,tt$chromatin_state), mean)
      df$grp = x
      colnames(df) = c("sample_name","chromatin_state","expression","grp")
      return(df)
    })
    mean = do.call(rbind,list_mean)

    p = ggplot(mean,aes(x = factor(grp, levels = lim), y = expression,
      color = chromatin_state,
      group = sample_name, linetype = sample_name)) +
      geom_point(size = 3) +
      geom_line(size = 2) +
      scale_color_manual(values = col$stateNumber) +
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
