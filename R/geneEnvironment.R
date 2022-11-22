#' Function which estimated the chromatin state in the environment around gene (estimated with interval around gene TSS)
#'
#' @title geneEnvironment
#' @param geneExpressionTable table with gene expression and gene Ensembl name (ex. geneExpression)
#' @param tableChromatinState table with information on the chromatin state (ex. chromatinState)
#' @param stateOrder list of chromatin state
#' @param interval numeric value corresponding to environment distance (default = 3000 (3kb))
#'
#' @examples
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' chromatinState = system.file("extdata", chromatinState, package = "ChromENVEE")
#' data(chromatinState)
#' geneExpression = system.file("extdata", geneExpression, package = "ChromENVEE")
#' data(geneExpression)
#' colorTable = system.file("extdata", colorTable, package = "ChromENVEE")
#' data(colorTable)
#' geneEnvironment(geneExpression, chromatinState, unique(colorTable$stateName))
#'
#' @return table with the distribution of each chromatin state in the environment
#' @export
geneEnvironment = function(geneExpressionTable,tableChromatinState, stateOrder, interval = 3000) {

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

	list_tableChromatinState = lapply(unique(tableChromatinState$chr),function(chr) {
		tt = tableChromatinState[tableChromatinState$chr == chr,]
		return(tt)
	})
	names(list_tableChromatinState) = unique(tableChromatinState$chr)

	geneExpressionTable[,stateOrder] = do.call(rbind,lapply(rownames(geneExpressionTable), function(gene) {
		gene_value = geneExpressionTable[gene,]
		chr = gene_value$chr

		tableChromatinState_chr = list_tableChromatinState[[chr]]

		## Remove chromatin state before interval
		tt = tableChromatinState_chr[!((tableChromatinState_chr$start < gene_value$TSS_moins_3kb)
		 	& (tableChromatinState_chr$end < gene_value$TSS_moins_3kb)),]
		## Remove chromatin state after interval
		tt = tt[!((tt$start > gene_value$TSS_plus_3kb)
			& (tt$end > gene_value$TSS_plus_3kb)),]

		tt$coverage = unlist(lapply(rownames(tt),function(state) {
			tt_state = tt[state,]
			if((tt_state$start > gene_value$TSS_moins_3kb) & (tt_state$end < gene_value$TSS_plus_3k)){
				length_overlapping = abs(tt_state$start - tt_state$end)
			} else if((tt_state$start < gene_value$TSS_moins) & (tt_state$end < gene_value$TSS_plus)){
				length_overlapping = abs(gene_value$TSS_moins - tt_state$end)
			} else {
				length_overlapping = abs(tt_state$start - gene_value$TSS_plus)
			}
			return(length_overlapping)
		}))

		return(t(data.frame(unlist(lapply(stateOrder, function(x){
			sum(tt[tt$state_name == x,"coverage"])/(interval*2)
		})))))
	}))
	return(geneExpressionTable)
}
