#' Function with associate at each enhancer the gene expression level
#'
#' @param enhancer_table GRanges table contains genomic position return by ChromHMM tool
#' @param gene_expression_table a table contains gene expression and gene Ensembl gene (RNAseq)
#'
#' @return GRanges table with for each genomic position gene associate and the expression of gene
#' @export
geneExpression = function(enhancer_table, gene_expression_table) {

  ## Return gene common between gene expression and enhancer_table
  common_gene = intersect(gene_expression_table$gene_ENS,
    unlist(strsplit(enhancer_table$gene_ENS,";")))

  ## Return position common gene in enhancer table
  pos_common_gene_enhancer = unlist(lapply(common_gene,list_pos_gene_common_enhancer,
    enhancer_table$gene_ENS))

  ## Reduce enhancer_table with only enhancer associate gene in gene expression
  common_gene_enhancer_table = enhancer_table[unique(pos_common_gene_enhancer),]
  common_gene_enhancer_table$gene_expression = "NA"

  ## Return expression of each common gene
  expression = sapply(common_gene, function(gene,gene_table,enhancer_table) {
    expression = gene_table[gene_table$gene_ENS == gene,"gene_expression"]
  },gene_table = gene_expression_table,enhancer_table = common_gene_enhancer_table)

  ## Return position of each common gene in enhancer table reduce
  pos_common_gene_enhancer_table = unlist(lapply(common_gene_enhancer_table$gene_ENS,
    function(gene, enhancer_gene) {
      gene2 = unlist(strsplit(gene,";"))
      pos = which(enhancer_gene %in% gene2)
      return(pos)
  }, enhancer_gene = names(expression)))


  pos_common_gene_enhancer_table = unlist(lapply(common_gene_enhancer_table$gene_ENS,
    function(enhancer_gene, gene) {
      list_enhancer_gene = unlist(strsplit(enhancer_gene,";"))
      # print(list_enhancer_gene)
      pos = which(gene %in% list_enhancer_gene)

      ## Take in count enhancer associate many genes and randomly select gene expression
      if(length(pos) > 1) {
        pos = sample(pos,1)
      }
      return(pos)
  }, gene = names(expression)))


  common_gene_enhancer_table$gene_expression = expression[pos_common_gene_enhancer_table]

  return(common_gene_enhancer_table)
}

list_pos_gene_common_enhancer = function(gene,enhancer_gene) {
  pos = which(grepl(gene,enhancer_gene))
  return(pos)
}
