

setwd('/Users/kkalantar/Documents/Research/NMLHC/geo_data')


#
# Generate the ANCESTRY dataset
#

ANCESTRY_metadata1 <- GEOquery::getGEO(filename='GSE81046-GPL11154_series_matrix.txt')
ANCESTRY_metadata2 <- GEOquery::getGEO(filename='GSE81046-GPL16791_series_matrix.txt')

ANCESTRY_p1 <- pData(ANCESTRY_metadata1)
ANCESTRY_p2 <- pData(ANCESTRY_metadata2)
extra_col <- colnames(ANCESTRY_p1)[!(colnames(ANCESTRY_p1) %in% colnames(ANCESTRY_p2))]
ANCESTRY_p1 <- ANCESTRY_p1[,colnames(ANCESTRY_p1) != extra_col]

ANCESTRY_meta <- rbind(ANCESTRY_p1, ANCESTRY_p2)

ANCESTRY_genecounts <- read.table("GSE81046_genes.tsv", sep="\t", row.names = 1, header = TRUE)

ANCESTRY_genecounts <- ANCESTRY_genecounts[, rownames(ANCESTRY_meta)]
ANCESTRY_DATASET <- ExpressionSet(assayData = as.matrix(ANCESTRY_genecounts), phenoData = AnnotatedDataFrame(ANCESTRY_meta) )

colors <- c("magenta","turquoise","orange")
p <- prcomp(t(exprs(ANCESTRY_DATASET)))
plot(p$x[,1],p$x[,2], col = colors[ANCESTRY_DATASET$characteristics_ch1.2] ,cex=1.2, lwd = 2, main = "ANCESTRY", xlab = "PC1",ylab = "PC2")
saveRDS( ANCESTRY_DATASET, "ANCESTRY_DATASET.rds" )



#
# Generate the ANTIVIRAL dataset
# 

ANTIVIRAL_metadata <- GEOquery::getGEO(filename = "GSE92904_series_matrix.txt")
ANTIVIRAL_genecounts <- read.table("GSE92904_genes.tsv", sep = "\t", row.names = 1, header = TRUE)

ANTIVIRAL_DATASET <- ExpressionSet(assayData = as.matrix(ANTIVIRAL_genecounts), phenoData = AnnotatedDataFrame(pData(ANTIVIRAL_metadata)) )

colors <- c("magenta","turquoise","orange")
p <- prcomp(t(exprs(ANTIVIRAL_DATASET)))
plot(p$x[,1],p$x[,3], col = colors[ANTIVIRAL_DATASET$characteristics_ch1.3] ,cex=1.2, lwd = 2, main = "ANTIVIRAL", xlab = "PC1",ylab = "PC3")

saveRDS( ANTIVIRAL_DATASET, "ANTIVIRAL_DATASET.rds" )



#
# Generate the TB_BLOOD dataset
#


TBBLOOD_metadata <- GEOquery::getGEO(filename = "GSE94438_series_matrix.txt")
TBBLOOD_genecounts <- read.table("GSE94438_genes.tsv", sep=",", row.names = 1, header = TRUE)
TBBLOOD_genecounts <- TBBLOOD_genecounts[,!(colnames(TBBLOOD_genecounts) %in% "symbol")]

TBBLOOD_metadata$characteristics_ch1.9 <- unlist(lapply(TBBLOOD_metadata$characteristics_ch1.1, function(x){gsub("code: ", "X", x)})) 
TBBLOOD_pdata <- pData(TBBLOOD_metadata)
library(dplyr)
TBBLOOD_pdata_unique <- TBBLOOD_pdata %>% distinct(characteristics_ch1.9, .keep_all = TRUE)

TBBLOOD_pdata_unique$characteristics_ch1.9 == colnames(TBBLOOD_genecounts)

TBBLOOD_genecounts <- TBBLOOD_genecounts[,TBBLOOD_pdata_unique$characteristics_ch1.9] 
colnames(TBBLOOD_genecounts) <- TBBLOOD_pdata_unique$geo_accession
rownames(TBBLOOD_pdata_unique) <- TBBLOOD_pdata_unique$geo_accession

TBBLOOD_DATASET <- ExpressionSet(assayData = as.matrix(TBBLOOD_genecounts), phenoData = AnnotatedDataFrame(TBBLOOD_pdata_unique))

colors <- c("magenta","turquoise","orange")
p <- prcomp(t(exprs(TBBLOOD_DATASET)))
plot(p$x[,1],p$x[,2], col = colors[TBBLOOD_DATASET$characteristics_ch1.6] ,cex=1.2, lwd = 2, main = "TBBLOOD", xlab = "PC1",ylab = "PC2")

saveRDS( TBBLOOD_DATASET, "TBBLOOD_DATASET.rds" )
