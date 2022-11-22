#' Function which associated at each enhancer the gene expression level
#'
#' @title enhancerExpression
#' @param enhancerTable GRanges object contains genomic position
#' @param geneExpressionTable gene expression and gene Ensembl name (RNAseq data)
#'
#' @examples
#' listTableEnhancer = system.file("extdata", listTableEnhancer, package = "ChromENVEE")
#' data(listTableEnhancer)
#' genomeFile = system.file("extdata", genomeFile, package = "ChromENVEE")
#' data(genomeFile)
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' enhancerExpression(anno, geneExpression)
#'
#' @return GRanges table with for each genomic position all gene associated and the expression of these genes
#' @export
enhancerExpression = function(enhancerTable, geneExpressionTable) {

  enhancerTable$gene_expression = unlist(lapply(seq_len(length(enhancerTable)), function(enhancer){
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
