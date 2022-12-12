#' Function with estimate the predominant state associated at each gene
#'
#' @title predominantState
#' @param table GRanges with the output of geneEnvironment
#' @param state list of chromatin state
#' @param header list of column in table to analysed with UMAP
#' @param neighbors numeric value (see umap package)
#' @param metric chracter value (see umap package)
#' @param dist numeric value (see umap package)
#'
#' @import umap
#'
#' @examples
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' chromatinState = system.file("extdata", chromatinState, package = "ChromENVEE")
#' data(chromatinState)
#' colorTable = system.file("extdata", colorTable, package = "ChromENVEE")
#' data(colorTable)
#' enviro = geneEnvironment(geneExpression[1:100,], chromatinState, unique(colorTable$stateName))
#' predominantState(enviro,state = unique(colorTable$stateName),
#'   header = unique(colorTable$stateName))
#'
#' @return table with 2-dimensional value and predominant chromatin state
#' @export
predominantState = function(table,
    state,
    header,
    neighbors = 32,
    metric = "euclidean",
    dist = 0.5) {

  message("==> It will be take few minutes to process")
  table = data.frame(table)
  result = invisible(umap(table[,header], n_neighbors = neighbors,
      metric = metric, min_dist = dist))
  colnames(result$layout) = c("UMAP1","UMAP2")

  result = cbind(result$data,result$layout)
  result = as.data.frame(result)

  result$state = sapply(rownames(result), function(gene) {
  	pos = which(table[gene,state] == max(table[gene,state]))
  	x = colnames(table[,state])[pos[1]]
  	return(x)
  })

  return(result)
}
