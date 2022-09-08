#' Function which associated at each enhancer the gene expression level
#'
#' @param enhancerTable GRanges table contains genomic position (ex. return by ChromHMM tool)
#' @param geneExpressionTable a table contains gene expression and gene Ensembl name (RNAseq data)
#'
#' @return GRanges table with for each genomic position all gene associated and the expression of these genes
#' @export
enhancerExpression = function(enhancerTable, geneExpressionTable) {

  enhancerTable$gene_expression = unlist(lapply(1:length(enhancerTable), function(enhancer){
    list_gene = enhancerTable[enhancer]$gene_list
    list_gene = unlist(strsplit(list_gene,";"))

    expression = unlist(lapply(list_gene, function(gene) {
      x = geneExpressionTable[geneExpressionTable$gene_ENS == gene,"gene_expression"]
      if(length(x) != 1) {
        x = NA
      }
      return(x)
    }))
    return(paste0(expression, collapse = ";"))
  }))
  return(enhancerTable)
}
