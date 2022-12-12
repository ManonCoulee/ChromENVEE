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
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' chromatinState = system.file("extdata", chromatinState, package = "ChromENVEE")
#' data(chromatinState)
#' colorTable = system.file("extdata", colorTable, package = "ChromENVEE")
#' data(colorTable)
#' geneEnvironment(geneExpression[1:100,], chromatinState, unique(colorTable$stateName))
#'
#' @return GRanges with the distribution of each chromatin state in the environment
#' @export
geneEnvironment = function(geneExpressionTable,tableChromatinState, stateOrder,
	 	interval = 3000) {

	geneExpressionTable$TSS = as.numeric(apply(geneExpressionTable,1,function(line) {
	  if(line["strand"] == "+") {
	    return(line["start"])
	  } else {
	    return(line["end"])
	  }
	}))

	geneExpressionTable$TSS_moins_3kb = geneExpressionTable$TSS - interval
	geneExpressionTable$TSS_plus_3kb = geneExpressionTable$TSS + interval

	geneExpressionTable[,stateOrder] = 0

	list_tableChromatinState = lapply(unique(seqnames(tableChromatinState)),function(chr) {
		tt = tableChromatinState[seqnames(tableChromatinState) == chr,]
		return(tt)
	})
	names(list_tableChromatinState) = unique(seqnames(tableChromatinState))

	geneExpressionTable[,stateOrder] = do.call("rbind",lapply(rownames(geneExpressionTable),
			function(gene, geneExpressionTable,list_tableChromatinState) {

		m3kb = geneExpressionTable[gene,"TSS_moins_3kb"]
		p3kb = geneExpressionTable[gene,"TSS_plus_3kb"]
		chr = geneExpressionTable[gene,"chr"]

		tableChromatinState_chr = list_tableChromatinState[[chr]]

		## Remove chromatin state before interval
		tt = tableChromatinState_chr[end(tableChromatinState_chr) > m3kb,]
		## Remove chromatin state after interval
		tt = tt[start(tt) < p3kb,]

		tt$coverage = unlist(lapply(seq_along(tt),function(state) {
			tt_state = tt[state,]
			S = GenomicRanges::start(tt_state)
			E = GenomicRanges::end(tt_state)

			if(E > p3kb) {
				length_overlapping = abs(S - p3kb)
			} else if(S < m3kb){
				length_overlapping = abs(m3kb - E)
			} else {
				length_overlapping = abs(S - E)
			}
			return(length_overlapping)
		}))/(interval*2)

		value = geneExpressionTable[gene,stateOrder]
		t = unlist(lapply(unique(tt$state_name), function(x){
			sum(tt[tt$state_name == x,]$coverage)
		}))
		names(t) = unique(tt$state_name)
		value[,names(t)] = t
		return(value)

	}, geneExpressionTable = geneExpressionTable,
		list_tableChromatinState = list_tableChromatinState))
	return(geneExpressionTable)
}
