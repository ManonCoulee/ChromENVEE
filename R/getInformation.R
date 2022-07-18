#' Function to return information (expression and distance) for each enhancer
#'
#' @param data_table a GRanges table or list of GRanges obtains by enhancerExpression function
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the distribution of distance gene-enhancer
#' @export
getInformation = function(data_table) {
  if(class(data_table) == "list") {
    df = lapply(data_table,function(x) {
      gene_name = unlist(strsplit(unlist(x$gene_list),";"))
      expression = unlist(strsplit(unlist(x$gene_expression),";"))
      distance = unlist(strsplit(unlist(x$distance),";"))
      df = data.frame(gene_name,expression,distance)
      df$sample_name = unique(x$sample_name)
      group = unique(x$chromatin_state)
      df$chromatin_state =  unlist(lapply(group, function(state) {
        count = sum(x[x$chromatin_state == state]$gene_association)
        rep(state,count)
      }))
      df = df[df$expression != "NA",]
      df$expression = as.numeric(df$expression)

      return(df)
    })

    data_frame = do.call(rbind,df)
    return(data_frame)

  } else if(class(data_table) == "GRanges") {
    gene_name = unlist(strsplit(unlist(data_table$gene_list),";"))
    expression = unlist(strsplit(unlist(data_table$gene_expression),";"))
    distance = unlist(strsplit(unlist(data_table$distance),";"))
    data_frame = data.frame(gene_name,expression,distance)
    data_frame$sample_name = unique(data_table$sample_name)
    group = unique(data_table$chromatin_state)
    data_frame$chromatin_state =  unlist(lapply(group, function(state) {
      count = sum(data_table[data_table$chromatin_state == state]$gene_association)
      rep(state,count)
    }))
    data_frame = data_frame[data_frame$expression != "NA",]
    return(data_frame)
  } else {
    stop("'data_table' must be a GRanges object or a list of GRanges object")
  }
}
