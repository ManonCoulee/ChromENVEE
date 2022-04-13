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
plotEnhancerExpression = function(data_table, xlab = "", ylab = "gene expression log10(cpm)",
  scale = "none", model = 1) {

  information_table = get_information(data_table)
  information_table$value = as.numeric(information_table$value)

  if(scale == "log10") {
    information_table$value = log10(information_table$value+0.01)
  } else if(scale == "log") {
    information_table$value = log(information_table$value)
  }

  if(model == 1) {
    color_value = colorTable$states15
  } else if(model == 2) {
    color_value = colorTable$states18
  } else if(model == 3) {
    color_value = colorTable$states25
  }

  p = ggplot(information_table,aes(x = sample_name, y = value)) +
    geom_violin(aes(fill = chromatin_state),color = "black") +
    geom_boxplot(width = 0.1) +
    scale_fill_manual(values = color_value) +
    xlab(xlab) + ylab(ylab) +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=35, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      legend.position = "none")

  return(p)
}

get_information = function(data_table) {
  if(class(data_table) == "list") {
    df = lapply(data_table,function(x) {
      value = unlist(strsplit(unlist(x$gene_expression),";"))
      df = data.frame(value)
      df$sample_name = unique(x$sample_name)
      group = unique(x$chromatin_state)
      df$chromatin_state =  unlist(lapply(group, function(state) {
        count = sum(x[x$chromatin_state == state]$gene_association)
        rep(state,count)
      }))
      return(df)
    })
    data_frame = do.call(rbind,df)
    data_frame = data_frame[data_frame$value != "NA",]
    return(data_frame)
  } else if(class(data_table) == "GRanges") {
    value = unlist(strsplit(unlist(data_table$gene_expression),";"))
    data_frame = data.frame(value)
    data_frame$sample_name = unique(data_table$sample_name)
    group = unique(data_table$chromatin_state)
    data_frame$chromatin_state =  unlist(lapply(group, function(state) {
      count = sum(data_table[data_table$chromatin_state == state]$gene_association)
      rep(state,count)
    }))
    data_frame = data_frame[data_frame$value != "NA",]
    return(data_frame)
  } else {
    print("'data_table' must be a GRanges object or a list of GRanges object")
  }
}
