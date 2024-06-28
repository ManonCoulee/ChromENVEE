#' Function which estimated the chromatin state in the environment around gene (estimated with interval around gene TSS)
#'
#' @title geneEnvironment
#' @param geneExpressionTable table with gene expression and gene Ensembl name (ex. geneExpression)
#' @param tableChromatinState GRanges with information on the chromatin state (ex. chromatinState)
#' @param stateOrder list of chromatin state
#' @param interval numeric value corresponding to environment distance (default = 3000 (3kb))
#'
#' @import GenomicRanges
#'
#' @examples
#' geneFile = system.file("extdata", "gene.tsv", package = "ChromENVEE")
#' geneExpression = read.table(geneFile, h = TRUE)
#' data(chromatinState)
#' data(colorTable)
#' geneEnvironment(geneExpression, chromatinState, unique(colorTable$stateName))
#'
#' @return GRanges with the distribution of each chromatin state in the environment
#' @export
geneEnvironment <- function(geneExpressionTable,
		tableChromatinState,
		stateOrder,
	 	interval = 3000) {

	geneExpressionTable$TSS <- geneExpressionTable$start
	geneExpressionTable[geneExpressionTable$strand == "-", "TSS"] <-
		geneExpressionTable[geneExpressionTable$strand == "-", "end"]
	# as.numeric(apply(geneExpressionTable,1,function(line) {
	#   if(line["strand"] == "+") {
	#     return(line["start"])
	#   } else {
	#     return(line["end"])
	#   }
	# }))

	geneExpressionTable$TSS3kbBefore <- geneExpressionTable$TSS - interval
	geneExpressionTable$TSS3kbAfter <- geneExpressionTable$TSS + interval

	geneExpressionTable[,stateOrder] <- 0

	listTableChromatinState <- lapply(unique(seqnames(tableChromatinState)), function(chr) {
		subTableChromatinState = tableChromatinState[seqnames(tableChromatinState) == chr,]
		return(subTableChromatinState)
	})
	names(listTableChromatinState) = unique(seqnames(tableChromatinState))

	geneExpressionTable[,stateOrder] = do.call("rbind", lapply(rownames(geneExpressionTable),
			function(geneName, geneExpressionTable, listTableChromatinState) {

		BeforePosGene <- geneExpressionTable[geneName, "TSS3kbBefore"]
		AfterPosGene <- geneExpressionTable[geneName, "TSS3kbAfter"]
		chr <- geneExpressionTable[geneName, "seqnames"]

		tableChromatinStateChr <- listTableChromatinState[[as.character(chr)]]

		## Remove chromatin state before interval
		tableRedChromatinState <- tableChromatinStateChr[end(tableChromatinStateChr) > BeforePosGene,]
		## Remove chromatin state after interval
		tableRedChromatinState <- tableRedChromatinState[start(tableRedChromatinState) < AfterPosGene,]

		s = start(tableRedChromatinState) < BeforePosGene
		start(tableRedChromatinState[s]) = BeforePosGene
		e = end(tableRedChromatinState) > AfterPosGene
		end(tableRedChromatinState[e]) = AfterPosGene

		tableRedChromatinState$coverage <- abs(GenomicRanges::start(tableRedChromatinState) - GenomicRanges::end(tableRedChromatinState)) / (interval*2)

		endStateFirst <- GenomicRanges::end(tableRedChromatinState[1])
		startStateLast <- GenomicRanges::start(tableRedChromatinState[length(tableRedChromatinState)])
		tableRedChromatinState[1]$coverage <- abs(BeforePosGene - endStateFirst) / (interval*2)
		tableRedChromatinState[length(tableRedChromatinState)]$coverage <- abs(startStateLast - AfterPosGene) / (interval*2)

		value <- geneExpressionTable[geneName, stateOrder]
		# stateCoverage <- unlist(lapply(unique(stateOrder), function(state){
		# 	sum(tableRedChromatinState[tableRedChromatinState$stateName == state,]$coverage)
		# }))
		tableCoverage <- aggregate(tableRedChromatinState$coverage,
				by = list(tableRedChromatinState$stateName), sum)
		value[,as.character(tableCoverage$`Group.1`)] <- tableCoverage$x

		# names(stateCoverage) = unique(stateOrder)
		# value[,names(stateCoverage)] = stateCoverage
		return(value)

	}, geneExpressionTable = geneExpressionTable, listTableChromatinState = listTableChromatinState))

	return(geneExpressionTable)
}
