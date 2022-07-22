#' Function with estimate the predominant state associated at each gene
#'
#' @param table a table corresponding to the result of geneEnvironment function
#' @param state a list of chromatin state
#' @param header a list of column in table to analysed with UMAP
#' @param neighbors a numeric value (see umapr package)
#' @param metric a chracter value (see umapr package)
#' @param dist a numeric value (see umapr package)
#'
#' @import umapr
#'
#' @return the table with 2-dimensional value and predominent chromatin state
#' @export
predominentState = function(table, state, header, neighbors = 32, metric = "euclidean", dist = 0.5) {

  result = umap(table[,header], n_neighbors = neighbors, metric = metric, min_dist = dist)

  result$state = sapply(rownames(result), function(gene) {
  	pos = which(table[gene,state] == max(table[gene,state]))
  	x = colnames(table[,state])[pos[1]]
  	return(x)
  })

  return(result)
}
