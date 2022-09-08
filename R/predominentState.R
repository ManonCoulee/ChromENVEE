#' Function with estimate the predominant state associated at each gene
#'
#' @param table a table corresponding to the result of geneEnvironment function
#' @param state a list of chromatin state
#' @param header a list of column in table to analysed with UMAP
#' @param neighbors a numeric value (see umapr package)
#' @param metric a chracter value (see umapr package)
#' @param dist a numeric value (see umapr package)
#'
#' @import umap
#'
#' @return table with 2-dimensional value and predominent chromatin state
#' @export
predominentState = function(table, state, header, neighbors = 32, metric = "euclidean", dist = 0.5) {

  message("\n==> It will be take few minutes to process\n")
  result = invisible(umap(table[,header], n_neighbors = neighbors, metric = metric, min_dist = dist))
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
