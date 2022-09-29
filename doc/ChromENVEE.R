## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Loading package
library(ChromENVEE)

## -----------------------------------------------------------------------------
stateNumber = c("U1","U2","U3","U4","U5","U6","U7","U8","U9","U10","U11","U12","U13","U14","U15","U16","U17","U18")
stateName = c("TSSA","TSSFlnk","TSSFlnkD","Tx","TxWk","EnhG","EnhG","EnhA","EnhWk","ZNF/Rpts","Het","TssBiv","EnhBiv","ReprPC","ReprPCWk","Quies","Quies","Quies")
colorValue = c("#B71C1C","#E65100","#E65100","#43A047","#1B5E20","#99FF66","#99FF66","#F5B041","#FFEB3B","#48C9B0","#B39DDB","#880E4F","#666633","#424949","#7B7D7D","#D0D3D4","#D0D3D4","#D0D3D4")

## -----------------------------------------------------------------------------
data(genomeFile)

## ----echo = FALSE-------------------------------------------------------------
head(genomeFile)

## -----------------------------------------------------------------------------
data(chromatinState)

## ----echo = FALSE-------------------------------------------------------------
head(chromatinState)

## -----------------------------------------------------------------------------
summary_chromatin_state = plotChromatinState(chromatinState, stateName = stateName, stateNumber = stateNumber, merge = TRUE, plot = FALSE, color = colorValue, filename = "")

## -----------------------------------------------------------------------------
head(summary_chromatin_state)

## -----------------------------------------------------------------------------
data(listTableEnhancer)

## ----echo = FALSE-------------------------------------------------------------
listTableEnhancer[[1]]

## -----------------------------------------------------------------------------
table_enhancer_gene = enhancerAnnotation(listTableEnhancer[[1]],genome = genomeFile, interval = 500000, nCore = 1)

## ----echo = FALSE-------------------------------------------------------------
table_enhancer_gene
table_enhancer_gene$sample_name = "RS"

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotGeneAssociation(table_enhancer_gene, all = FALSE)

## -----------------------------------------------------------------------------
data(geneExpression)

## ----echo = FALSE-------------------------------------------------------------
head(geneExpression)

## -----------------------------------------------------------------------------
table_enhancer_gene_expression = enhancerExpression(table_enhancer_gene, geneExpressionTable = geneExpression)

## ----echo = FALSE-------------------------------------------------------------
table_enhancer_gene_expression

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotDistanceExpression(table_enhancer_gene_expression, color = colorValue, stateName = stateName, stateNumber = stateNumber)

## ---- fig.width = 7,fig.asp = 0.3---------------------------------------------
plotGeneDistance(table_enhancer_gene_expression)

## ---- fig.width = 7,fig.asp = 0.6---------------------------------------------
plotEnhancerExpression(table_enhancer_gene_expression, scale = "log10", color = colorValue, stateName = stateName, stateNumber = stateNumber, ylab = "gene expression log10(CPM)")

## -----------------------------------------------------------------------------
# list_table_enhancer_gene = lapply(listTableEnhancer, enhancerAnnotation, genome = genomeFile, interval = 500000, nCore = 1)

# list_table_enhancer_gene_expression = lapply(list_table_enhancer_gene, enhancerExpression, gene_expression_table = geneExpression)

## -----------------------------------------------------------------------------
# plotGeneAssociation(list_table_enhancer_gene_sample, all = T)

## -----------------------------------------------------------------------------
# plotDistanceExpression(list_table_enhancer_gene_expression, color = colorValue, stateName = stateName, stateNumber = stateNumber)

## -----------------------------------------------------------------------------
# plotGeneDistance(list_table_enhancer_gene_expression)

## -----------------------------------------------------------------------------
# plotEnhancerExpression(list_table_enhancer_gene_expression, scale = "log10", color = colorValue, stateName = stateName, stateNumber = stateNumber, ylab = "gene expression log10(CPM)")

## -----------------------------------------------------------------------------
stateOrderReduce = c("TSSA","TSSFlnk","TSSFlnk","Tx","Tx","EnhG","EnhG","EnhA","EnhWk","ZNF.Rpts","Het","TssBiv","EnhBiv","ReprPC","ReprPC","Quies","Quies","Quies")

## -----------------------------------------------------------------------------
data(geneExpression)
data(chromatinState)

## -----------------------------------------------------------------------------
table_overlapping = geneEnvironment(geneExpression, chromatinState, stateOrder = unique(stateOrderReduce), interval = 3000)
rownames(table_overlapping) = table_overlapping$gene_ENS

## ----echo = FALSE-------------------------------------------------------------
head(table_overlapping)

## -----------------------------------------------------------------------------
result_umap = predominentState(table_overlapping, state = unique(stateOrderReduce),
 header = unique(stateOrderReduce), neighbors = 32, metric = "euclidean", dist = 0.5)

## ----echo = FALSE-------------------------------------------------------------
head(result_umap)

## ----echo = FALSE-------------------------------------------------------------
sessionInfo()

