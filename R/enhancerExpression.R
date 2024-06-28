#' Function which associated at each enhancer the gene expression level
#'
#' @title enhancerExpression
#' @param enhancerTable GRanges object contains genomic position
#' @param geneExpressionTable gene expression and gene Ensembl name (RNAseq data)
#'
#' @examples
#' data(listTableEnhancer)
#' data(genomeFile)
#' data(geneExpression)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' enhancerExpression(anno, geneExpression)
#'
#' @return GRanges table with for each genomic position all gene associated and the expression of these genes
#' @export
enhancerExpression <- function(enhancerTable, geneExpressionTable) {

  enhancerTable$geneExpression <- unlist(lapply(seq_len(length(enhancerTable)), function(enhancer) {
    listGene <- enhancerTable[enhancer]$geneList
    listGene <- unlist(strsplit(listGene, ";"))

    expression <- unlist(lapply(listGene, function(gene) {
      subGeneExpressionTable <- geneExpressionTable[geneExpressionTable$geneENS == gene, "geneExpression"]
      if(length(subGeneExpressionTable) != 1) {
        subGeneExpressionTable <- NA
      }
      return(subGeneExpressionTable)
    }))
    return(paste0(expression, collapse = ";"))
  }))
  return(enhancerTable)
}
