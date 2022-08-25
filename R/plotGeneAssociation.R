#' Function which plot the number of gene associated to enhancer
#'
#' @param table GRanges table contains results obtains by enhancerAnnotation()
#' @param all a boolean, if TRUE, list of dataframe is merge
#'
#' @import ggplot2
#'
#' @return a plot of the number of gene associated to enhancer
#' @export
plotGeneAssociation = function(table, all = F) {

  if (all == T && class(table) == "list") {
    list_table = lapply(names(table), function(name) {
      tt = table[[name]]
      df = data.frame(table(tt$gene_association))
      colnames(df) = c("gene_association","count")
      df$gene_association = as.numeric(df$gene_association)
      df$name = name
      return(df)
    })
    df = do.call(rbind,list_table)

    ## Interpolation polynomiale
    p = ggplot(df,aes(x = as.numeric(gene_association), y = count, color = name)) +
      geom_smooth(method = "lm",formula = y~poly(x,11),se = F) +
      xlab("number of gene associated") + ylab("frequence") + labs(color = "") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")
  } else {

    df = data.frame(table(table$gene_association))
    colnames(df) = c("gene_association","count")

    p = ggplot(df,aes(x = as.numeric(gene_association), y = count)) +
      geom_point() +
      geom_smooth(method = "lm",formula = y~poly(x,10),se = F,color = "red") +
      xlab("number of gene associated") + ylab("frequence") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")
  }
  return(p)
}
