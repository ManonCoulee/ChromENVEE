#' A list of 11 GRanges data frame, each data frame contains enhancer chromatin state information and gene associated information
#'
#' @format A list of 8 GRanges data frame:
#' \describe{
#'  \item{seqnames}{chromosome number}
#'  \item{ranges}{start and end position, in bp}
#'  \item{strand}{strand value}
#'  \item{chromatin_state}{chromatin state number}
#'  \item{sample_name}{name of the sample}
#'  \item{start_500kb}{position 500kb before start position}
#'  \item{end_500kb}{position 500kb after end position}
#'  \item{gene_association}{number of gene associated enhancer}
#'  \item{distance}{list of distance gene-enhancer in bp}
#'  \item{gene_list}{list of gene name associated enhancer}
#'  \item{gene_expression}{list of gene expression associated enhancer}
#'}
#'
"listTableEnhancerGeneExpression"
