#' Function with associate at each enhancer the nearly gene
#'
#' @param enhancer_table GRanges table contains genomic position return by ChromHMM tool
#' @param genome a bed genome annotation file
#' @param interval a numeric value corresponding to the distance from enhancer where gene is looking for
#'
#' @import ggplot2
#'
#' @return the table with infromation concerning gene associate enhancer
#' @export
plotGeneAssociation = function(table, all = F, state_name, state_number, state_color) {

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
      xlab("number of gene associate") + ylab("frequence") + labs(color = "") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")

      # theme_bw() + theme(strip.background  = element_blank(),
      #   text = element_text(size=35, angle = 0),
      #   panel.grid.major = element_line(colour = "grey80"),
      #   panel.border = element_blank(),
      #   axis.ticks = element_blank(),
      #   panel.grid.minor.x=element_blank(),
      #   panel.grid.major.x=element_blank(),
      #   legend.text = element_text(size = 15),
      #   legend.position = "bottom")
  } else {

    df = data.frame(table(table$gene_association))
    colnames(df) = c("gene_association","count")

    p = ggplot(df,aes(x = as.numeric(gene_association), y = count)) +
      geom_point() +
      geom_smooth(method = "lm",formula = y~poly(x,10),se = F,color = "red") +
      xlab("number of gene associate") + ylab("frequence") +
      themePlot() +
      theme(axis.text.x = element_text(),
        axis.text.y = element_text(),
        legend.position = "bottom")

      # theme_bw() + theme(strip.background  = element_blank(),
      #   text = element_text(size=35, angle = 0),
      #   panel.grid.major = element_line(colour = "grey80"),
      #   panel.border = element_blank(),
      #   axis.ticks = element_blank(),
      #   panel.grid.minor.x=element_blank(),
      #   panel.grid.major.x=element_blank(),
      #   legend.position = "bottom")
  }
  return(p)
}
