#' Function which associated at each enhancer all gene in the interval arount the enhancer
#'
#' @title enhancerAnnotation
#' @param enhancerTable GRanges object or list of GRanges object with enhancer position
#' @param genome GRanges of genome annotation file
#' @param interval distance from enhancer where gene is looking for (default = 500000 (500kb))
#' @param nCore number of core used (default nCore = 1)
#'
#' @import GenomicRanges
#' @import parallel
#' @importFrom methods as
#' @importFrom methods is
#'
#' @examples
#' listTableEnhancer = system.file("extdata", listTableEnhancer, package = "ChromENVEE")
#' data(listTableEnhancer)
#' genomeFile = system.file("extdata", genomeFile, package = "ChromENVEE")
#' data(genomeFile)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' anno
#'
#' @return GRanges with infromation about gene associated with enhancer
#' @export
enhancerAnnotation = function(enhancerTable, genome,
		interval = 500000,
		nCore = 1){

  if(!is(enhancerTable,"GRanges")){
    stop("'enhancerTable' is not a GRanges table")
  }

  ## Add interval value at each side enhancer
  enhancerTable$start_500kb = GenomicRanges::start(enhancerTable) - interval
  enhancerTable$end_500kb = GenomicRanges::end(enhancerTable) + interval

  ## Add TSS position in function strand
  genome$TSS = GenomicRanges::start(genome)
	genome[GenomicRanges::strand(genome)@values == "-",]$TSS = end(
			genome[GenomicRanges::strand(genome)@values == "-",])

  ## Divide gene table in function of the chromosome
  list_genome_table = lapply(unique(seqnames(genome)),function(chr){
		tt = genome[seqnames(genome) == chr,]
		tt = as.data.frame(tt)[,c("TSS","gene_ENS")]
    return(tt)
  })
  names(list_genome_table) = unique(seqnames(genome))

  # Return gene associate enhancer
  list_sub_enhancerTable = split(seq_len(length(enhancerTable)),
    rep_len(seq_len(nCore), length(enhancerTable)))

  list_enhancerTable = parallel::mclapply(list_sub_enhancerTable, function(pos,table){
    tt = table[pos,]

    results = lapply(seq_along(tt), comparisonPositionGeneEnhancer, enhancer_table = tt,
      list_genome_table = list_genome_table)

    tt$gene_association = unlist(lapply(results,function(x){return(x[1])}))
    tt$distance = unlist(lapply(results,function(x){return(x[2])}))
    tt$gene_list = unlist(lapply(results,function(x){return(x[3])}))
    return(tt)
  }, table = enhancerTable, mc.cores = nCore)

  table_final = unlist(methods::as(list_enhancerTable,"GRangesList"))
  return(table_final)
}

comparisonPositionGeneEnhancer = function(enhancer,enhancer_table,list_genome_table){

  ## Chromosome associate position
  chr_value_enhancer = GenomicRanges::seqnames(enhancer_table[enhancer,])@values

  genome = list_genome_table[[as.character(chr_value_enhancer)]]

  ## Position of enhancer
  start = GenomicRanges::start(enhancer_table[enhancer,])
  end = GenomicRanges::end(enhancer_table[enhancer,])
  start_500kb =  enhancer_table[enhancer,]$start_500kb
  end_500kb = enhancer_table[enhancer,]$end_500kb

  tt = genome[genome$TSS < start,]
  sub_genome = tt[tt$TSS > start_500kb,]
  sub_genome$distance = abs(sub_genome$TSS - start)
  tt = genome[genome$TSS > end,]
  sub_genome2 = tt[tt$TSS < end_500kb,]
  sub_genome2$distance = abs(sub_genome2$TSS - end)
  genome_res = rbind(sub_genome,sub_genome2)

  nb_gene = nrow(genome_res)
  distance = paste0(genome_res$distance, collapse = ";")
  gene = paste0(genome_res$gene_ENS, collapse = ";")

  return(list(nb_gene,distance,gene))
}
