#' Function with associate at each enhancer the nearly gene
#'
#' @param enhancer_table GRanges table contains genomic position return by ChromHMM tool
#' @param genome a bed genome annotation file
#' @param interval a numeric value corresponding to the distance from enhancer where gene is looking for
#'
#' @import GenomicFeatures
#'
#' @return the table with infromation concerning gene associate enhancer
#' @export
enhancerAnnotation = function(enhancer_table,genome, interval = 500000, nCore = 1) {

  set.seed(10)
  if(class(enhancer_table) != "GRanges"){
    print("'enhancer_table' is not a GRanges table")
  }
  ## Add interval value at each side enhancer
  enhancer_table$start_500kb = start(enhancer_table) - interval
  enhancer_table$end_500kb = end(enhancer_table) + interval

  ## Divide gene table in function of the chromosome
  list_genome_table = lapply(unique(genome$chr),function(chr) {
    tt = genome[genome$chr == chr,c("TSS","gene_ENS")]
    return(tt)
  })
  names(list_genome_table) = unique(genome$chr)

  # Return gene associate enhancer
  list_sub_enhancer_table = split(1:length(enhancer_table),
    rep_len(1:nCore, length(enhancer_table)))

  list_enhancer_table = mclapply(list_sub_enhancer_table, function(pos,table) {
    tt = table[pos,]
    # results = lapply(1:length(tt),count_number_gene_associate,tt,genome_gene_file)
    results = lapply(1:length(tt), comparison_position_gene_enhancer, tt, list_genome_table)

    tt$gene_association = unlist(lapply(results,function(x) {return(x[1])}))
    tt$distance = unlist(lapply(results,function(x){return(x[2])}))
    tt$gene_list = unlist(lapply(results,function(x){return(x[3])}))
    return(tt)
  }, table = enhancer_table, mc.cores = nCore)

  table_final = unlist(as(list_enhancer_table,"GRangesList"))
  return(table_final)
}

# comparison_position_gene_enhancer = function(genome,start,start_500kb,end,end_500kb){
#
#   TSS = as.numeric(genome["start"])
#
#   ## Position of enhancer
#   start = as.numeric(start)
#   end = as.numeric(end)
#   start_500kb = as.numeric(start_500kb)
#   end_500kb = as.numeric(end_500kb)
#
#   if (TSS < start && TSS > start_500kb ) { distance = abs(start - TSS) }
#   else if (TSS > end && TSS < end_500kb) { distance = abs(end - TSS) }
#   else { distance = NA }
#   return(distance)
# }
#
# find_gene_associate_enhancer = function(row,enhancer,genome) {
#
#   chr_value_enhancer = seqnames(enhancer[row,])@values
#   chr_genome = genome[genome$chr == chr_value_enhancer,]
#
#   ## Return position from gene in function enhancer
#   results = unlist(apply(chr_genome,1,comparison_position_gene_enhancer,start(enhancer[row,]),
#     enhancer[row,]$start_500kb,end(enhancer[row,]),enhancer[row,]$end_500kb))
#   chr_genome$distance = results
#
#   ## Return gene with minimal distance from enhancer start or end
#   chr_genome = chr_genome[!is.na(chr_genome$distance),]
#   distance_minimal = min(chr_genome$distance,na.rm = T)
#   gene_associate_enhancer = chr_genome[chr_genome$distance == distance_minimal,"gene_ENS"]
#   if(length(gene_associate_enhancer) == 0) {
#     return(list("NA","NA"))
#   } else if(length(gene_associate_enhancer) > 1){
#     return(list(paste0(gene_associate_enhancer,collapse = ";"),distance_minimal))
#   } else {
#     return(list(gene_associate_enhancer,distance_minimal))
#   }
# }

## Nouvelle methode
# comparison_position_gene_enhancer = function(gene,genome_table,start,start_500kb,end,end_500kb){
#
#   table = genome_table[gene,]
#   ## Return TSS position of gene
#   if(table$strand == "+") {
#     TSS = table$start
#   } else {
#     TSS = table$end
#   }
#
#   ## Position of enhancer
#   start = as.numeric(start)
#   end = as.numeric(end)
#   start_500kb = as.numeric(start_500kb)
#   end_500kb = as.numeric(end_500kb)
#
#   if (TSS < start && TSS > start_500kb ) {
#     distance = abs(start - TSS)
#     name = table$gene_ENS
#     return(paste0(c(distance,name), collapse = ";"))
#   } else if (TSS > end && TSS < end_500kb) {
#     distance = abs(end - TSS)
#     name = table$gene_ENS
#     return(paste0(c(distance,name), collapse = ";"))
#   }
# }

comparison_position_gene_enhancer = function(enhancer,enhancer_table,list_genome_table){

  ## Chromosome associate position
  chr_value_enhancer = seqnames(enhancer_table[enhancer,])@values

  genome = list_genome_table[[chr_value_enhancer]]

  ## Position of enhancer
  start = start(enhancer_table[enhancer,])
  end = end(enhancer_table[enhancer,])
  start_500kb =  enhancer_table[enhancer,]$start_500kb
  end_500kb = enhancer_table[enhancer,]$end_500kb

  tt = genome[genome$TSS < start,]
  sub_genome = tt[tt$TSS > start_500kb,]
  sub_genome$distance = abs(sub_genome$TSS - start)
  tt = genome[genome$TSS > end,]
  sub_genome2 = tt[tt$TSS < end_500kb,]
  sub_genome2$distance = abs(sub_genome2$TSS - end)
  genome_res = rbind(sub_genome,sub_genome2)

  nb_gene = nrow(genome_res)
  distance = paste0(genome_res$distance, collapse = ";")
  gene = paste0(genome_res$gene_ENS, collapse = ";")

  return(list(nb_gene,distance,gene))
}
