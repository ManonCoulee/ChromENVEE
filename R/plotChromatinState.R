#' Function to return summary table from ChromHMM data and create associated plot
#'
#' @title plotChromatinState
#' @param tableChromatinState GRanges object containing information about the chromatin state (ex. chromatinState data)
#' @param colorTable dataframe which contains color information (ex. colorTable data)
#' @param plot boolean to create plot (default = TRUE)
#' @param merge boolean to merge data if it's a list of table (default = F), if TRUE, list of dataframe is merge
#' @param filename string to name the plot create (default = "chromatin_state_distribution")
#' @param ylab x-axis label (default = "chromatin state (%)")
#' @param xlab y-axis label (default = "")
#'
#' @importFrom methods is
#' @import GenomicRanges
#'
#' @examples
#' data(chromatinState)
#' data(colorTable)
#' state = plotChromatinState(chromatinState, colorTable)
#' state
#'
#' @return table contains distribution of different chromatin state
#' @export
plotChromatinState = function(tableChromatinState, colorTable,
					plot = TRUE,
					merge = FALSE,
					filename = "chromatin_state_distribution",
					ylab = "chromatin state (%)",
					xlab = "") {

	if(is(tableChromatinState,"GRanges")) {
		table = distributionChromatinState(tableChromatinState,colorTable)
	} else if(is(tableChromatinState, "list")) {
		table = lapply(tableChromatinState,distributionChromatinState,colorTable)
	} else {
		stop("'tableChromatinState' must be a GRanges or a list of GRanges")
	}

	if(merge && is(table, "list")) {
		table = do.call("rbind",lapply(table,"[",,))
	}

	if(plot) {
		if(is(table, "list")) {
			lapply(table,plotDistributionChromatinState,filename = filename,
				colorTable = colorTable,
				merge = merge,ylab = ylab,xlab = xlab)
		} else {
			plotDistributionChromatinState(table,filename = filename,
				colorTable = colorTable,
				merge = merge,ylab = ylab,xlab = xlab)
		}
	}
	return(table)
}

distributionChromatinState = function(tableChromatinState,colorTable) {

	if(is(tableChromatinState, "GRanges")) {
		resume = data.frame("state" = unique(colorTable$stateName))
		resume$state = factor(resume$state, levels = unique(colorTable$stateName))
		rownames(resume) = unique(colorTable$stateName)

		tableChromatinState$size = abs(GenomicRanges::start(tableChromatinState) -
			GenomicRanges::end(tableChromatinState))
		tableChromatinState$state_name = factor(tableChromatinState$state,
				levels = colorTable$stateNumber, labels = colorTable$stateName)
		genome_length = sum(tableChromatinState$size)

		resume$coverage = unlist(lapply(rownames(resume), function(state) {
				(sum(tableChromatinState[tableChromatinState$state_name == state,]$size)/genome_length)*100
		}))
		resume$sample_name = unique(tableChromatinState$name)

		return(resume)
	} else {
		stop("'tableChromatinState' must be a GRanges or a list of GRanges")
	}
}

plotDistributionChromatinState = function(table, filename, colorTable,
	merge,ylab,xlab) {

	col = getStateColor(colorTable = colorTable)

	p = ggplot(table, aes(y = coverage, x = factor(sample_name))) +
		geom_bar(aes(fill = factor(state),
			alpha = factor(sample_name)),
			position = "dodge",stat = "identity") +
		facet_grid(.~factor(state), switch = "x") +
		scale_alpha_manual(values = seq(0.5,1, (0.5/(length(unique(table$sample_name))-1)))) +
		scale_fill_manual(values = col$stateName) +
		xlab(xlab) + ylab(ylab) +
		labs(fill = "chomatin state", alpha = "cell type") +
		themePlot() +
		theme(strip.text.x = element_text(size=35, angle = 90, hjust = 1),
			axis.text.y = element_text(size=35, angle = 90, hjust = 1),
			legend.position="bottom",
			legend.text = element_text(size=30, angle = 0),
	    legend.title = element_text(size=30, angle = 0)) +
		guides(fill = "none")

	ggsave(filename = filename,plot = p, width = 20, height = 12, device = 'png', dpi = 300)
}
