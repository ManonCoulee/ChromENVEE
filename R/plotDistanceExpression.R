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
plotDistanceExpression = function(data_table, xlab = "", ylab = "log(CPM)", model = 1) {

  limit = c(0,50000,100000,150000,200000,250000)
  limit_label = c("0-50kb","50-100kb","100-150kb","150-200kb","200-250kb",">250kb")

  information_table = get_information(data_table)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  tt = information_table[information_table$distance < 100000,]

  p = ggplot(tt,aes(x = sample_name, y = log(expression+0.01))) +
    geom_boxplot(aes(fill = sample_name),color = "black") +
    scale_fill_manual(values = rep(c("#F5B041","#666633","#99FF66","#FFEB3B"),each = 2)) +
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


  # for(l in limit) {
  #   pos = which(information_table$distance > l)
  #   information_table[pos,"limit"] = limit_label[limit == l]
  # }
  #
  # if(model == 1) {
  #   color_value = colorTable$states15
  # } else if(model == 2) {
  #   # color_value = colorTable$states18
  #   color_value = rep(c("#F5B041","#666633","#99FF66","#FFEB3B"),each = 2)
  # } else if(model == 3) {
  #   color_value = colorTable$states25
  # }
  #
  # if(class(data_table) == "list") {
  #   ylab = "mean(CPM)"
  #   list_mean = lapply(limit_label, function(x) {
  #     tt = information_table[information_table$limit == x,]
  #     df = aggregate(tt$expression,list(tt$sample_name), mean)
  #     df$grp = x
  #     colnames(df) = c("sample_name","expression","grp")
  #     # df$name = paste0(df$state,df$sample_name)
  #     return(df)
  #   })
  #   mean = do.call(rbind,list_mean)
  #
  #   p = ggplot(mean,aes(x = factor(grp, levels = limit_label), y = expression, color = sample_name,
  #     group = sample_name, linetype = sample_name)) +
  #     geom_point(size = 3) +
  #     geom_line(size = 2) +
  #     scale_color_manual(values = color_value) +
  #     scale_linetype_manual(values = rep(c("solid","dashed"),4)) +
  #     xlab(xlab) + ylab(ylab) +
  #     theme_bw() + theme(strip.background  = element_blank(),
  #       text = element_text(size=35, angle = 0),
  #       panel.grid.major = element_line(colour = "grey80"),
  #       panel.border = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.text.x = element_text(size = 30, angle = 0, hjust = 1),
  #       panel.grid.minor.x=element_blank(),
  #       panel.grid.major.x=element_blank(),
  #       legend.position = "none")
  # } else {
  #   p = ggplot(information_table,aes(x = factor(limit, levels = limit_label), y = log(expression+0.01))) +
  #     geom_boxplot(aes(fill = chromatin_state),width = 0.1) +
  #     scale_fill_manual(values = color_value) +
  #     xlab(xlab) + ylab(ylab) +
  #     theme_bw() + theme(strip.background  = element_blank(),
  #       text = element_text(size=35, angle = 0),
  #       panel.grid.major = element_line(colour = "grey80"),
  #       panel.border = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
  #       panel.grid.minor.x=element_blank(),
  #       panel.grid.major.x=element_blank(),
  #       legend.position = "none")
  # }

  return(p)
}

get_information = function(data_table) {
  if(class(data_table) == "list") {
    df = lapply(data_table,function(x) {
      expression = unlist(strsplit(unlist(x$gene_expression),";"))
      distance = unlist(strsplit(unlist(x$distance),";"))
      df = data.frame(expression,distance)
      df$sample_name = unique(x$sample_name)
      group = unique(x$chromatin_state)
      df$chromatin_state =  unlist(lapply(group, function(state) {
        count = sum(x[x$chromatin_state == state]$gene_association)
        rep(state,count)
      }))

      df = df[df$distance < 100000,]
      df$expression = as.numeric(df$expression)

      return(df)
    })

    data_frame = do.call(rbind,df)
    data_frame = data_frame[data_frame$expression != "NA",]
    return(data_frame)
  } else if(class(data_table) == "GRanges") {
    expression = unlist(strsplit(unlist(data_table$gene_expression),";"))
    distance = unlist(strsplit(unlist(data_table$distance),";"))
    data_frame = data.frame(expression,distance)
    data_frame$sample_name = unique(data_table$sample_name)
    group = unique(data_table$chromatin_state)
    data_frame$chromatin_state =  unlist(lapply(group, function(state) {
      count = sum(data_table[data_table$chromatin_state == state]$gene_association)
      rep(state,count)
    }))
    data_frame = data_frame[data_frame$expression != "NA",]
    return(data_frame)
  } else {
    print("'data_table' must be a GRanges object or a list of GRanges object")
  }
}
