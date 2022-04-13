#' Function with associate at each enhancer the gene expression level
#'
#' @param enhancer_table GRanges table contains genomic position return by ChromHMM tool
#' @param gene_expression_table a table contains gene expression and gene Ensembl gene (RNAseq)
#'
#' @return GRanges table with for each genomic position gene associate and the expression of gene
#' @export
enhancerExpression = function(enhancer_table, gene_expression_table) {

  print("data")
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
