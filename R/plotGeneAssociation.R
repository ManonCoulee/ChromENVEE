#' Function which plot the number of gene associated to enhancer
#'
#' @title plotGeneAssociation
#' @param table GRanges object or list of GRanges output of enhancerAnnotation
#' @param merge boolean, if TRUE, list of dataframe is merge (default = FALSE)
#'
#' @import ggplot2
#'
#' @examples
#' data(listTableEnhancer)
#' data(genomeFile)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' plotGeneAssociation(anno)
#'
#' @return ggplot2 figure of the number of gene associated to enhancer
#' @export
plotGeneAssociation <- function(enhancerTable, merge = FALSE) {

  if (merge && is(enhancerTable, "list")) {
    listTableAssociation <- lapply(names(enhancerTable), function(name) {
      subEnhancerTable <- enhancerTable[[name]]
      countAssociationTable <- data.frame(table(subEnhancerTable$geneAssociation))
      colnames(countAssociationTable) <- c("geneAssociation", "count")
      countAssociationTable$geneAssociation <- as.numeric(countAssociationTable$geneAssociation)
      countAssociationTable$name <- name
      return(countAssociationTable)
    })
    tableAssociation = do.call(rbind, listTableAssociation)

    ## Interpolation polynomiale
    plot <- ggplot(tableAssociation, aes(x = as.numeric(geneAssociation), y = count, color = name)) +
      geom_smooth(method = "lm", formula = y ~ poly(x, 11), se = FALSE) +
      xlab("number of gene associated") + ylab("frequence") + labs(color = "") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")
  } else {

    tableAssociation <- data.frame(table(enhancerTable$geneAssociation))
    colnames(tableAssociation) <- c("geneAssociation","count")

    plot <- ggplot(tableAssociation, aes(x = as.numeric(geneAssociation), y = count)) +
      geom_point() +
      geom_smooth(method = "lm", formula = y ~ poly(x, 10), se = FALSE, color = "red") +
      xlab("number of gene associated") + ylab("frequence") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")
  }
  return(plot)
}
