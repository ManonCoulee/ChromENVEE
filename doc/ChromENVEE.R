## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Loading package
library(ChromENVEE)

## -----------------------------------------------------------------------------
data(colorTable)

## -----------------------------------------------------------------------------
data(genomeFile)

## ----echo = FALSE-------------------------------------------------------------
head(genomeFile)

## -----------------------------------------------------------------------------
data(chromatinState)

## ----echo = FALSE-------------------------------------------------------------
head(chromatinState)

## -----------------------------------------------------------------------------
summary_chromatin_state = plotChromatinState(chromatinState, merge = TRUE, plot = FALSE,
colorTable = colorTable, filename = "")

## ----echo = FALSE-------------------------------------------------------------
head(summary_chromatin_state)

## -----------------------------------------------------------------------------
data(listTableEnhancer)

## ----echo = FALSE-------------------------------------------------------------
listTableEnhancer[[1]]

## -----------------------------------------------------------------------------
table_enhancer_gene = enhancerAnnotation(listTableEnhancer[[1]], genome = genomeFile,
interval = 500000, nCore = 1)

## ----echo = FALSE-------------------------------------------------------------
head(table_enhancer_gene)

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotGeneAssociation(table_enhancer_gene, all = FALSE)

## -----------------------------------------------------------------------------
data(geneExpression)

## ----echo = FALSE-------------------------------------------------------------
head(geneExpression)

## -----------------------------------------------------------------------------
table_enhancer_gene_expression = enhancerExpression(table_enhancer_gene,
geneExpressionTable = geneExpression)

## ----echo = FALSE-------------------------------------------------------------
head(table_enhancer_gene_expression)

## ---- fig.width = 10,fig.asp = 0.3--------------------------------------------
plotGeneDistance(table_enhancer_gene_expression, limit = 500000, xlab = "",
ylab = "distance enhancer-gene (bp)")

## ---- fig.width = 10,fig.asp = 0.4--------------------------------------------
plotEnhancerExpression(table_enhancer_gene_expression, scale = "log10",
colorTable = colorTable, ylab = "gene expression log10(CPM)")

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotDistanceExpression(table_enhancer_gene_expression, colorTable = colorTable,
limit = 500000)

## -----------------------------------------------------------------------------
list_table_enhancer_gene = lapply(listTableEnhancer, enhancerAnnotation,
genome = genomeFile, interval = 500000, nCore = 1)
listTableEnhancerGeneExpression = lapply(list_table_enhancer_gene, enhancerExpression,
geneExpressionTable = geneExpression)

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotGeneAssociation(listTableEnhancerGeneExpression, all = TRUE)

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotGeneDistance(listTableEnhancerGeneExpression, limit = 500000,
xlab = "", ylab = "distance enhancer-gene (bp)")

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotEnhancerExpression(listTableEnhancerGeneExpression, scale = "log10",
colorTable = colorTable, ylab = "gene expression log10(CPM)")

## ---- fig.width = 10,fig.asp = 0.5--------------------------------------------
plotDistanceExpression(listTableEnhancerGeneExpression, colorTable = colorTable,
limit = 500000)

## -----------------------------------------------------------------------------
data(geneExpression)
data(chromatinState)

## -----------------------------------------------------------------------------
table_overlapping = geneEnvironment(geneExpression, chromatinState,
stateOrder = unique(colorTable$stateName), interval = 3000)
rownames(table_overlapping) = table_overlapping$gene_ENS

## ----echo = FALSE-------------------------------------------------------------
head(table_overlapping)

## -----------------------------------------------------------------------------
result_umap = predominantState(table_overlapping, state = unique(colorTable$stateName),
header = unique(colorTable$stateName), neighbors = 32, metric = "euclidean", dist = 0.5)

## ----echo = FALSE-------------------------------------------------------------
head(result_umap)

## ----echo = FALSE-------------------------------------------------------------
sessionInfo()

