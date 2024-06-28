#' Function to return the theme usable for each plot
#'
#' @return the theme usable with ggplot
themePlot <- function() {
  theme <- theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=20, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())
  return(theme)
}
