## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ChromENVEE)

## -----------------------------------------------------------------------------
state_number = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14",
  "U15","U16","U17","U18")
state_name = c("TSS","TSSFlnk","TSSFlnkD","Tx","TxWk","EnhG","EnhG","EnhA","EnhWk",
  "ZNF/Rpts","Het","TssBiv","EnhBiv","ReprPC","ReprPCWk","Quies","Quies","Quies")
color_value = c("#B71C1C","#E65100","#E65100","#43A047","#1B5E20","#99FF66","#99FF66", "#F5B041","#FFEB3B","#48C9B0","#B39DDB","#880E4F","#666633","#424949","#7B7D7D","#D0D3D4",
  "#D0D3D4","#D0D3D4")

## -----------------------------------------------------------------------------
data(genome_gene_file)
head(genome_gene_file)


## -----------------------------------------------------------------------------
data(data_chromatin_state)
head(data_chromatin_state)

## -----------------------------------------------------------------------------
summary_chromatin_state = plotChromatinState(data_chromatin_state,state_name, state_number,
  merge = T, plot = T, color = color_value, filename = "figures/chromatin_state_distribution.png")

## -----------------------------------------------------------------------------
head(summary_chromatin_state)

## -----------------------------------------------------------------------------
data(list_table_enhancer)
list_table_enhancer[[1]]

## -----------------------------------------------------------------------------
table_enhancer_gene = enhancerAnnotation(list_table_enhancer[[1]],genome = genome_gene_file,interval = 500000, nCore = 1)
table_enhancer_gene

## -----------------------------------------------------------------------------
table_enhancer_gene$sample_name = "RS"
table_enhancer_gene

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotGeneAssociation(table_enhancer_gene, all = F)

## -----------------------------------------------------------------------------
gene_expression = read.table("~/Documents/Mouse/RNAseq/DE_analysis_3_factors_Kit_SC_SCII_RS_spike/RS_CTL_gene_expression.tsv",sep = "\t",h = T)
head(gene_expression)
table_enhancer_gene_expression = enhancerExpression(table_enhancer_gene,gene_expression_table = gene_expression)
table_enhancer_gene_expression

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotDistanceExpression(table_enhancer_gene_expression,color = color_value,
    state_name = state_name, state_number = state_number)

## ---- fig.width = 7,fig.asp = 0.3---------------------------------------------
plotGeneDistance(table_enhancer_gene_expression)

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotEnhancerExpression(table_enhancer_gene_expression, scale = "log10",
  color = color_value, state_name = state_name, state_number = state_number)

## -----------------------------------------------------------------------------
state_order_reduce = c("TSS","TSSFlnk","TSSFlnk","Tx","Tx","EnhG","EnhG","EnhA","EnhW","ZNF.Rpts",
  "Het","TssBiv","EnhBiv","ReprPC","ReprPC","Quies","Quies","Quies")

## -----------------------------------------------------------------------------
data(table_gene_bed)
head(table_gene_bed)


## -----------------------------------------------------------------------------
#table_overlapping = geneEnvironment(table_gene_bed,data_chromatin_state, unique(state_order_reduce))
#rownames(table_overlapping) = table_overlapping$gene_ENS

## ----echo = FALSE-------------------------------------------------------------
#head(table_overlapping)

## -----------------------------------------------------------------------------
#result_umap = predominentState(table_overlapping, state = unique(state_order_reduce),
#  header = unique(state_order_reduce) ,neighbors = 32, metric = "euclidean", dist = 0.5)
#head(result_umap)

## ----echo = FALSE-------------------------------------------------------------
sessionInfo()

