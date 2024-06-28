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
#' data(listTableEnhancer)
#' data(genomeFile)
#' data(geneExpression)
#' data(colorTable)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' expression = enhancerExpression(anno, geneExpression)
#' plotDistanceExpression(expression, colorTable = colorTable)
#'
#' @return gene expression in function distance gene-enhancer
#' @export
plotDistanceExpression <- function(enhancerTable,
					xlab = "",
          ylab = "log(CPM)",
          colorTable,
          limit = 500000) {

  col <- getStateColor(colorTable = colorTable)

  lim <- getLengthVector(limit)

  informationTable <- getInformation(enhancerTable)
  informationTable$distance <- as.numeric(informationTable$distance)
  informationTable$expression <- as.numeric(informationTable$expression)

  for(limit in names(lim)) {
    pos <- which(informationTable$distance > as.numeric(limit))
    informationTable[pos,"limit"] <- lim[limit]
  }

  if(is(enhancerTable, "list")) {
    ylab <- "mean(CPM)"
    listMean <- lapply(lim, function(limit) {
      limitTable <- informationTable[informationTable$limit == limit,]
      grpExpressionTable <- aggregate(limitTable$expression, list(limitTable$sampleName), mean)
      grpExpressionTable$grp <- limit
      colnames(grpExpressionTable) <- c("sampleName","expression","grp")
      return(grpExpressionTable)
    })
    meanTable <- do.call(rbind,listMean)

    plot <- ggplot(meanTable, aes(x = factor(grp, levels = lim), y = expression,
      color = sampleName,
      group = sampleName, linetype = sampleName)) +
      geom_point(size = 3) +
      geom_line(size = 2) +
      scale_color_manual(values = col$stateName) +
      scale_linetype_manual(values = rep(c("solid", "dashed"),
				round(length(unique(meanTable$sampleName))/2))) +
      xlab(xlab) + ylab(ylab) +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "none")
  } else {
    plot <- ggplot(informationTable, aes(x = factor(limit, levels = lim), y = log(expression + 0.01))) +
      geom_boxplot(aes(fill = chromatinState), width = 0.1) +
      scale_fill_manual(values = col$stateNumber) +
      xlab(xlab) + ylab(ylab) +
      themePlot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(),
        legend.position = "none")
  }

  return(plot)
}
