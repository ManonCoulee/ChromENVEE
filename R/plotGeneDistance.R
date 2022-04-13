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
plotGeneDistance = function(data_table, xlab = "", ylab = "distance enhancer-gene (bp)") {

  limit = c(0,50000,100000,150000,200000,250000)
  limit_label = c("0-50kb","50-100kb","100-150kb","150-200kb","200-250kb",">250kb")

  information_table = get_information(data_table)
  information_table$value = as.numeric(information_table$value)

  for(l in limit) {
    pos = which(information_table$value > l)
    information_table[pos,"distance"] = limit_label[limit == l]
  }

  p = ggplot(information_table,aes(x = sample_name, fill = factor(distance,levels = limit_label))) +
    geom_bar(stat = "count", position = "fill") +
    coord_flip() + labs(fill = "") + xlab(xlab) + ylab(ylab) +
    theme_bw() + theme(strip.background  = element_blank(),
      text = element_text(size=35, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x= element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      legend.position = "bottom") +
    scale_fill_manual(values = c("#f1948a","#c39bd3","#85c1e9","#76d7c4","#f7dc6f","#f0b27a"))

  return(p)
}

get_information = function(data_table) {
  if(class(data_table) == "list") {
    df = lapply(data_table,function(x) {
      value = unlist(strsplit(unlist(x$distance),";"))
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
    value = unlist(strsplit(unlist(data_table$distance),";"))
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
