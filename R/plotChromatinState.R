#' Function to return summary table from ChromHMM data and create associated plot
#'
#' @param tableChromatinState a bed table containing information about the chromatin state (ex. chromatin_state data)
#' @param colorTable a data frame which contains color information
#' @param plot a boolean to create plot (default = T)
#' @param merge a boolean to merge data if it's a list of table (default = F), if TRUE, list of dataframe is merge
#' @param filename a string to name the plot create (defualt = "chromatin_state_distribution")
#' @param ylab a string (default = "chromatin state (%)")
#' @param xlab a string (default = "")
#'
#' @return table contains distribution of different chromatin state
#' @export
plotChromatinState = function(tableChromatinState, colorTable,
					plot = T,
					merge = F,
					filename = "chromatin_state_distribution",
					ylab = "chromatin state (%)",
					xlab = "") {

	if(class(tableChromatinState) == "data.frame") {
		table = distributionChromatinState(tableChromatinState,stateName,stateNumber)
	} else if(class(tableChromatinState) == "list") {
		table = lapply(tableChromatinState,distribution_chromatin_state,stateName,stateNumber)
	} else {
		stop("'tableChromatinState' must be a data frame or a list of data frame")
	}

	if(merge == T & class(table) == "list") {
		table = do.call("rbind",lapply(table,"[",,))
	}

	if(plot == T) {
		if(class(table) == "list") {
			lapply(table,plot_distribution_chromatin_state,filename = filename, stateName = stateName,
				color = color, stateNumber = stateNumber,
				merge = merge,ylab = ylab,xlab = xlab)
		} else {
			plotDistributionChromatinState(table,filename = filename, stateName = stateName,
				color = color, stateNumber = stateNumber,
				merge = merge,ylab = ylab,xlab = xlab)
		}
	}
	return(table)
}

distributionChromatinState = function(tableChromatinState,stateName, stateNumber) {

	if(class(tableChromatinState) == "data.frame") {
		resume = data.frame("state" = unique(stateName))
		resume$state = factor(resume$state, levels = unique(stateName))
		rownames(resume) = unique(stateName)

		tableChromatinState$size = abs(tableChromatinState$start - tableChromatinState$end)
		tableChromatinState$state_name = factor(tableChromatinState$state, levels = stateNumber, labels = stateName)
		genome_length = sum(tableChromatinState$size)

		resume$coverage = unlist(lapply(rownames(resume), function(state) {
				(sum(tableChromatinState[tableChromatinState$state_name == state,"size"])/genome_length)*100
		}))
		resume$sample_name = unique(tableChromatinState$name)

		return(resume)
	} else {
		stop("'tableChromatinState' must be a data frame or a list of data frame")
	}
}

plotDistributionChromatinState = function(table, filename, color, stateName, stateNumber,
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
