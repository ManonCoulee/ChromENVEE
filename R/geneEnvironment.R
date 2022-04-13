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
state_overlapping = function(gene,state_order) {
  # print(gene[,"chr"])
  chr = gene["chr"]
  table_chromHMM_chr = table_chromHMM[table_chromHMM$chr == chr,]

	## Return the distance overlapped by the chromatin state
  list_overlap = unlist(apply(table_chromHMM_chr,1, function(state) {
		TSS_moins = as.numeric(gene["TSS_moins_3kb"])
		TSS_plus = as.numeric(gene["TSS_plus_3kb"])
	 	start= as.numeric(state["start"])
		end = as.numeric(state["end"])

  	if((start < TSS_moins) & (end < TSS_moins)) {
    } else if((start > TSS_plus) & (end > TSS_plus)) {
    } else {
			if((start > TSS_moins) & (end < TSS_plus)){
				length_overlapping = abs(start - end)
			} else if((start < TSS_moins) & (end < TSS_plus)){
				length_overlapping = abs(TSS_moins - end)
			} else {
				length_overlapping = abs(start - TSS_plus)
			}
      return(list("state" = state["state_name"], "overlap" = length_overlapping))
    }
	}))

	## Create a table with length associate each chromatin state in the environment
	table_list_overlap = data.frame(list_overlap)
	table_overlap = data.frame(table_list_overlap[seq(1,nrow(table_list_overlap),2),1])
	table_overlap$length = table_list_overlap[seq(2,nrow(table_list_overlap),2),1]
	colnames(table_overlap) = c("state","length")

	## Create count table
	table_overlap_reduce = data.frame(state_order)
	table_overlap_reduce$length = 0
	x = lapply(table_overlap_reduce$state_order, function(state) {
		return(sum(as.numeric(table_overlap[table_overlap$state == state,"length"])))
	})

  return(unlist(x))
}
