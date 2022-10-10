#' Function to return information (expression and distance) for each enhancer
#'
#' @param dataTable a GRanges table or list of GRanges obtains by enhancerExpression function
#'
#' @import ggplot2
#'
#' @return the ggplot2 figure corresponding to the distribution of distance gene-enhancer
#' @export
getInformation = function(dataTable) {
  if(class(dataTable) == "list") {
    df = lapply(dataTable,function(x) {
      if(class(x) != "GRanges") {
        stop("'dataTable' must be a list of GRanges object")
      }
      if(length(x$sample_name) == 0) {
        stop("GRanges object need a 'sample_name' column")
      }
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

  } else if(class(dataTable) == "GRanges") {
    if(length(dataTable$sample_name) == 0) {
      stop("GRanges object need a 'sample_name' column")
    }
    gene_name = unlist(strsplit(unlist(dataTable$gene_list),";"))
    expression = unlist(strsplit(unlist(dataTable$gene_expression),";"))
    distance = unlist(strsplit(unlist(dataTable$distance),";"))
    data_frame = data.frame(gene_name,expression,distance)
    group = unique(dataTable$chromatin_state)
    data_frame$chromatin_state =  unlist(lapply(group, function(state) {
      count = sum(dataTable[dataTable$chromatin_state == state]$gene_association)
      rep(state,count)
    }))
    data_frame = data_frame[data_frame$expression != "NA",]
    return(data_frame)
  } else {
    stop("'dataTable' must be a GRanges object or a list of GRanges object")
  }
}
