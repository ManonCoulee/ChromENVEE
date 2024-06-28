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
#' data(listTableEnhancer)
#' data(genomeFile)
#' anno = enhancerAnnotation(listTableEnhancer[[1]], genomeFile)
#' anno
#'
#' @return GRanges with infromation about gene associated with enhancer
#' @export
enhancerAnnotation <- function(enhancerTable,
		genome,
		interval = 500000,
		nCore = 1){

  if(!is(enhancerTable, "GRanges")) {
    stop("'enhancerTable' is not a GRanges table")
  }

  ## Add interval value at each side enhancer
  enhancerTable$start500kb <- GenomicRanges::start(enhancerTable) - interval
  enhancerTable$end500kb <-GenomicRanges::end(enhancerTable) + interval

  ## Add TSS position in function strand
  genome$TSS <- GenomicRanges::start(genome)
	genome[GenomicRanges::strand(genome)@values == "-", ]$TSS <- end(
			genome[GenomicRanges::strand(genome)@values == "-", ])

  ## Divide gene table in function of the chromosome
  listGenomeTable <- lapply(unique(seqnames(genome)), function(chr) {
		genomeTable <- genome[seqnames(genome) == chr, ]
		genomeTable <- as.data.frame(genomeTable)[,c("TSS", "geneENS")]
    return(genomeTable)
  })
  names(listGenomeTable) <- unique(seqnames(genome))

  # Return gene associate enhancer
  listSubEnhancerPosition <- split(seq_len(length(enhancerTable)),
    rep_len(seq_len(nCore), length(enhancerTable)))

  listEnhancerTable <- parallel::mclapply(listSubEnhancerPosition,
		function(position, enhancerTable) {

    subEnhancerTable <- enhancerTable[position, ]
    results <- lapply(seq_along(subEnhancerTable), comparisonPositionGeneEnhancer,
			subEnhancerTable = subEnhancerTable, listGenomeTable = listGenomeTable)

    subEnhancerTable$geneAssociation = unlist(lapply(results, function(x) {return(x[1])}))
    subEnhancerTable$distance = unlist(lapply(results, function(x) {return(x[2])}))
    subEnhancerTable$geneList = unlist(lapply(results, function(x) {return(x[3])}))

    return(subEnhancerTable)
  }, enhancerTable = enhancerTable, mc.cores = nCore)

  tableFinal <- unlist(methods::as(listEnhancerTable, "GRangesList"))

  return(tableFinal)
}

comparisonPositionGeneEnhancer <- function(seqLength, subEnhancerTable, listGenomeTable){

  ## Chromosome associate position
  chrValueEnhancer <- GenomicRanges::seqnames(subEnhancerTable[seqLength,])@values

  genome <- listGenomeTable[[as.character(chrValueEnhancer)]]
	if(is.null(genome)) {
		return(list(nbGene = 0,distance = NA, gene = NA))
	}

  ## Position of enhancer
  start <- GenomicRanges::start(subEnhancerTable[seqLength, ])
  end <- GenomicRanges::end(subEnhancerTable[seqLength, ])
  start500kb <-  subEnhancerTable[seqLength, ]$start500kb
  end500kb <- subEnhancerTable[seqLength, ]$end500kb

  genomeRed <- genome[genome$TSS < start,]
  genomeRed1 <- genomeRed[genomeRed$TSS > start500kb,]
  genomeRed1$distance <- abs(genomeRed1$TSS - start)

  genomeRed <- genome[genome$TSS > end,]
  genomeRed2 <- genomeRed[genomeRed$TSS < end500kb,]
  genomeRed2$distance <- abs(genomeRed2$TSS - end)
  genomeRedFinal <- rbind(genomeRed1,genomeRed2)

  nbGene <- nrow(genomeRedFinal)
	if(nbGene == 0){
		distance <- NA
		gene <- NA
	} else {
		distance <- paste0(genomeRedFinal$distance, collapse = ";")
		gene <- paste0(genomeRedFinal$geneENS, collapse = ";")
		# genomeMin <- genomeRedFinal[genomeRedFinal$distance == min(genomeRedFinal$distance),] ## RADA
		# if(nrow(genomeMin) == 2){ ## RADA
		# 	distance <- NA
		# 	gene <- NA
		# } else { ## RADA
		# 	distance <- genomeMin$distance
		# 	gene <- genomeMin$geneENS
		# }
	}
  return(list(nbGene, distance, gene))
}
