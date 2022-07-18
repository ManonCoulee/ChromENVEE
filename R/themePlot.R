#' Function to return summary table from ChromHMM data and create associated plot
#'
#' @param genomic_region a table contains genomic region associate each chromatin state (output from ChromHH)
#' @param state_name a vector of chromatin state name
#' @param state_number a vector of chromatin state number (name from ChromHMM)
#' @param color a vector of color to colored plot
#' @param plot a boolean to create plot
#' @param merge a boolean to merge data if it's a list of table
#' @param filename a string to name the plot create
#' @param ylab a string
#' @param xlab a string
#'
#' @return the table contains distribution of different chromatin state
#' @export
themePlot = function() {
  theme = theme_minimal() + theme(strip.background  = element_blank(),
      text = element_text(size=35, angle = 0),
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank())

  return(theme)
}
