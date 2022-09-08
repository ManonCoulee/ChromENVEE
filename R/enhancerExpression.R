#' Function which associated at each enhancer the gene expression level
#'
#' @param enhancer_table GRanges table contains genomic position (ex. return by ChromHMM tool)
#' @param gene_expression_table a table contains gene expression and gene Ensembl name (RNAseq data)
#'
#' @return GRanges table with for each genomic position all gene associated and the expression of these genes
#' @export
enhancerExpression = function(enhancer_table, gene_expression_table) {

  enhancer_table$gene_expression = unlist(lapply(1:length(enhancer_table), function(enhancer){
    list_gene = enhancer_table[enhancer]$gene_list
    list_gene = unlist(strsplit(list_gene,";"))

    expression = unlist(lapply(list_gene, function(gene) {
      x = gene_expression_table[gene_expression_table$gene_ENS == gene,"gene_expression"]
      if(length(x) != 1) {
        x = NA
      }
      return(x)
    }))
    return(paste0(expression, collapse = ";"))
  }))
  return(enhancer_table)
}
