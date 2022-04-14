#' Function with associate at each enhancer the nearly gene
#'
#' @param enhancer_table GRanges table contains genomic position return by ChromHMM tool
#' @param genome a bed genome annotation file
#' @param interval a numeric value corresponding to the distance from enhancer where gene is looking for
#'
#' @import ggplot2
#'
#' @return the table with infromation concerning gene associate enhancer
#' @export
state_overlapping = function(table,table_chromHMM, state_order, interval = 3000) {

	## Return TSS position (start if + strand or end if - strand)
	table$TSS = as.numeric(apply(table,1,function(line) {
	  if(line["strand"] == "+") {
	    return(line["start"])
	  } else {
	    return(line["end"])
	  }
	}))

	table$TSS_moins_3kb = table$TSS - interval
	table$TSS_plus_3kb = table$TSS + interval

	table[,state_order] = 0

	list_chromHMM_table = lapply(unique(table_chromHMM$chr),function(chr) {
		tt = table_chromHMM[table_chromHMM$chr == chr,]
		return(tt)
	})
	names(list_chromHMM_table) = unique(table_chromHMM$chr)

	table[,state_order] = do.call(rbind,lapply(rownames(table), function(gene) {
		gene_value = table[gene,]
		chr = gene_value$chr

		table_chromHMM_chr = list_chromHMM_table[[chr]]

		## Remove chromatin state before interval
		tt = table_chromHMM_chr[!((table_chromHMM_chr$start < gene_value$TSS_moins_3kb)
		 	& (table_chromHMM_chr$end < gene_value$TSS_moins_3kb)),]
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

		return(t(data.frame(unlist(lapply(state_order, function(x){
			sum(tt[tt$state_name == x,"coverage"])/(interval*2)
		})))))
	}))
	return(table)
}
