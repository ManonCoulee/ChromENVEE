#' Function to create a plot which represent the gene expression in function distance gene-enhancer
#'
#' @param data_table a GRanges table or list of GRanges obtains by enhancerExpression function
#' @param color a list of color value
#' @param state_number a list of chromatin state number
#' @param state_name a list of chromatin state name
#' @param xlab a string
#' @param ylab a string
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the mean of gene expression in function distance gene-enhancer
#' @export
plotDistanceExpression = function(data_table, xlab = "", ylab = "log(CPM)",color, state_number,
  state_name) {

  col = stateColor(state_name = state_name, state_number = state_number, color = color)

  limit = c(0,50000,100000,150000,200000,250000)
  limit_label = c("0-50kb","50-100kb","100-150kb","150-200kb","200-250kb",">250kb")

  information_table = getInformation(data_table)
  information_table$distance = as.numeric(information_table$distance)
  information_table$expression = as.numeric(information_table$expression)

  for(l in limit) {
    pos = which(information_table$distance > l)
    information_table[pos,"limit"] = limit_label[limit == l]
  }

  if(class(data_table) == "list") {
    ylab = "mean(CPM)"
    list_mean = lapply(limit_label, function(x) {
      tt = information_table[information_table$limit == x,]
      df = aggregate(tt$expression,list(tt$sample_name,tt$chromatin_state), mean)
      df$grp = x
      colnames(df) = c("sample_name","chromatin_state","expression","grp")
      return(df)
    })
    mean = do.call(rbind,list_mean)

    p = ggplot(mean,aes(x = factor(grp, levels = limit_label), y = expression,
      color = chromatin_state,
      group = sample_name, linetype = sample_name)) +
      geom_point(size = 3) +
      geom_line(size = 2) +
      scale_color_manual(values = col$state_number) +
      scale_linetype_manual(values = rep(c("solid","dashed"),length(unique(mean$sample_name))/2)) +
      xlab(xlab) + ylab(ylab) +
      theme_bw() + theme(strip.background  = element_blank(),
        text = element_text(size=35, angle = 0),
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 30, angle = 0, hjust = 1),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        legend.position = "none")
  } else {
    p = ggplot(information_table,aes(x = factor(limit, levels = limit_label), y = log(expression+0.01))) +
      geom_boxplot(aes(fill = chromatin_state),width = 0.1) +
      scale_fill_manual(values = col$state_number) +
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
  }

  return(p)
}
