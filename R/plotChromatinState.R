#' Function to return summary table from ChromHMM data and create associated plot
#'
#' @param genomic_region a table contains genomic region associate each chromatin state (output from ChromHH)
#' @param state_name a vector of chromatin state name
#' @param state_number a vector of chromatin state number (name from ChromHMM)
#' @param color a vector of color to colored plot
#' @param plot a boolean to create plot
#' @param merge a boolean to merge data if it's a list of table
#' @param filename a string to name the plot create
#' @param ylab a string
#' @param xlab a string
#'
#' @return the table contains distribution of different chromatin state
#' @export
plotChromatinState = function(genomic_region, state_name, state_number, color,
					plot = T,
					merge = F,
					filename = "chromatin_state_distribution",
					ylab = "chromatin state (%)",
					xlab = "") {

	if(class(genomic_region) == "data.frame") {
		table = distribution_chromatin_state(genomic_region,state_name,state_number)
	} else if(class(genomic_region) == "list") {
		table = lapply(genomic_region,distribution_chromatin_state,state_name,state_number)
	} else {
		stop("'genomic_region' must be a data frame or a list of data frame")
	}

	if(merge == T & class(table) == "list") {
		table = do.call("rbind",lapply(table,"[",,))
		# table = table[table$state != "Het",]
		# table = table[table$state != "ZNF/Rpts",]
	}

	if(plot == T) {
		if(class(table) == "list") {
			lapply(table,plot_distribution_chromatin_state,filename = filename, state_name = state_name,
				color = color, state_number = state_number,
				merge = merge,ylab = ylab,xlab = xlab)
		} else {
			plot_distribution_chromatin_state(table,filename = filename, state_name = state_name,
				color = color, state_number = state_number,
				merge = merge,ylab = ylab,xlab = xlab)
		}
	}
	return(table)
}

distribution_chromatin_state = function(genomic_region,state_name, state_number) {

	if(class(genomic_region) == "data.frame") {
		resume = data.frame("state" = unique(state_name))
		resume$state = factor(resume$state, levels = unique(state_name))
		rownames(resume) = unique(state_name)

		genomic_region$size = abs(genomic_region$start - genomic_region$end)
		genomic_region$state_name = factor(genomic_region$state, levels = state_number, labels = state_name)
		genome_length = sum(genomic_region$size)

		resume$coverage = unlist(lapply(rownames(resume), function(state) {
				(sum(genomic_region[genomic_region$state_name == state,"size"])/genome_length)*100
		}))
		resume$sample_name = unique(genomic_region$name)

		return(resume)
	} else {
		stop("'genomic_region' must be a data frame or a list of data frame")
	}
}

plot_distribution_chromatin_state = function(table, filename, color, state_name, state_number,
	merge,ylab,xlab) {

	col = stateColor(state_name = state_name, state_number = state_number, color = color)

	p = ggplot(table, aes(y = coverage, x = factor(sample_name, levels = c("Kit_m","Kit_p","SC","RS")))) +
		geom_bar(aes(fill = factor(state),
			alpha = factor(sample_name, levels = c("Kit_m","Kit_p","SC","RS"))),
			position = "dodge",stat = "identity") +
		facet_grid(.~factor(state), switch = "x") +
		scale_alpha_manual(values = seq(0.5,1, (0.5/(length(unique(table$sample_name))-1)))) +
		scale_fill_manual(values = col$state_name) +
		xlab(xlab) + ylab(ylab) +
		labs(fill = "chomatin state", alpha = "cell type") +
		theme_minimal() + theme(strip.background  = element_blank(),
				text = element_text(size=35, angle = 90),
				panel.grid.major = element_line(colour = "grey80"),
				panel.border = element_blank(),
				axis.ticks = element_blank(),
				axis.text.x = element_blank(),
				axis.text.y = element_text(size=35, angle = 90, hjust = 1),
				panel.grid.minor.x=element_blank(),
				panel.grid.major.x=element_blank()) +
		theme(legend.position="bottom", legend.text = element_text(size=30, angle = 0),
	    legend.title = element_text(size=30, angle = 0)) +
		guides(fill = "none")

	ggsave(filename = filename,plot = p, width = 20, height = 12, device = 'png', dpi = 300)
}
