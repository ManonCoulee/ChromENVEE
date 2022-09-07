## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Loading package
library(ChromENVEE)

## -----------------------------------------------------------------------------
state_number = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18")
state_name = c("TSSA","TSSFlnk","TSSFlnkD","Tx","TxWk","EnhG","EnhG","EnhA","EnhWk","ZNF/Rpts","Het","TssBiv","EnhBiv","ReprPC","ReprPCWk","Quies","Quies","Quies")
color_value = c("#B71C1C","#E65100","#E65100","#43A047","#1B5E20","#99FF66","#99FF66","#F5B041","#FFEB3B","#48C9B0","#B39DDB","#880E4F","#666633","#424949","#7B7D7D","#D0D3D4","#D0D3D4","#D0D3D4")

## -----------------------------------------------------------------------------
data(genome_file)

## ----echo = FALSE-------------------------------------------------------------
head(genome_file)

## -----------------------------------------------------------------------------
data(chromatin_state)

## ----echo = FALSE-------------------------------------------------------------
head(chromatin_state)

## -----------------------------------------------------------------------------
summary_chromatin_state = plotChromatinState(chromatin_state, state_name = state_name, state_number = state_number, merge = TRUE, plot = FALSE, color = color_value, filename = "")

## -----------------------------------------------------------------------------
head(summary_chromatin_state)

## -----------------------------------------------------------------------------
data(list_table_enhancer)

## ----echo = FALSE-------------------------------------------------------------
list_table_enhancer[[1]]

## -----------------------------------------------------------------------------
table_enhancer_gene = enhancerAnnotation(list_table_enhancer[[1]],genome = genome_file,interval = 500000, nCore = 1)

## ----echo = FALSE-------------------------------------------------------------
table_enhancer_gene
table_enhancer_gene$sample_name = "RS"

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotGeneAssociation(table_enhancer_gene, all = F)

## -----------------------------------------------------------------------------
data(geneExpression)

## ----echo = FALSE-------------------------------------------------------------
head(geneExpression)

## -----------------------------------------------------------------------------
table_enhancer_gene_expression = enhancerExpression(table_enhancer_gene, gene_expression_table = geneExpression)

## ----echo = FALSE-------------------------------------------------------------
table_enhancer_gene_expression

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotDistanceExpression(table_enhancer_gene_expression, color = color_value, state_name = state_name, state_number = state_number)

## ---- fig.width = 7,fig.asp = 0.3---------------------------------------------
plotGeneDistance(table_enhancer_gene_expression)

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotEnhancerExpression(table_enhancer_gene_expression, scale = "log10", color = color_value, state_name = state_name, state_number = state_number)

## -----------------------------------------------------------------------------
# list_table_enhancer_gene = lapply(list_table_enhancer,enhancerAnnotation, genome = genome_file, interval = 500000, nCore = 1)

# list_table_enhancer_gene_expression = lapply(list_table_enhancer_gene, enhancerExpression, gene_expression_table = geneExpression)

## -----------------------------------------------------------------------------
# plotGeneAssociation(list_table_enhancer_gene_sample, all = T)

## -----------------------------------------------------------------------------
# plotDistanceExpression(list_table_enhancer_gene_expression, color = color_value, state_name = state_name, state_number = state_number)

## -----------------------------------------------------------------------------
# plotGeneDistance(list_table_enhancer_gene_expression)

## -----------------------------------------------------------------------------
# plotEnhancerExpression(list_table_enhancer_gene_expression, scale = "log10", color = color_value, state_name = state_name, state_number = state_number)

## -----------------------------------------------------------------------------
state_order_reduce = c("TSSA","TSSFlnk","TSSFlnk","Tx","Tx","EnhG","EnhG","EnhA","EnhWk","ZNF.Rpts","Het","TssBiv","EnhBiv","ReprPC","ReprPC","Quies","Quies","Quies")

## -----------------------------------------------------------------------------
data(geneExpression)
data(chromatin_state)

## -----------------------------------------------------------------------------
table_overlapping = geneEnvironment(geneExpression, chromatin_state, state_order = unique(state_order_reduce), interval = 3000)
rownames(table_overlapping) = table_overlapping$gene_ENS

## ----echo = FALSE-------------------------------------------------------------
head(table_overlapping)

## -----------------------------------------------------------------------------
result_umap = predominentState(table_overlapping, state = unique(state_order_reduce),
 header = unique(state_order_reduce) ,neighbors = 32, metric = "euclidean", dist = 0.5)

## ----echo = FALSE-------------------------------------------------------------
head(result_umap)

## ----echo = FALSE-------------------------------------------------------------
sessionInfo()

