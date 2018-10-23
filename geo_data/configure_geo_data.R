
source("/Users/kkalantar/Documents/Research/NMLHC/microarray_format_utils.R")
library(hgu133a2.db)           # for conversion from probe IDs to gene IDs
library(hgu133plus2.db)        # ""
library(illuminaHumanv4.db)    # ""
library(sva)                   # for COMBAT batch correction


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





  

#sub <- subset_known_dataset(ANCESTRY_DATASET, "Listeria", "Non-infected", "characteristics_ch1.2")





#
# Tsalik original data
#

library(sva)
# variable: OG$characteristics_ch1.1; viral, bacterial, non-infectious
tsalik_OG <- GEOquery::getGEO(filename = "/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE63990_series_matrix.txt")
batch_corrected_OG <- ComBat(exprs(tslik_OG),tslik_OG$`batch:ch1`)
genenames_batch_corrected_OG <- detect_genecnt_platform(batch_corrected_OG)

TSALIK_ORIGINAL <- ExpressionSet(assayData = as.matrix(genenames_batch_corrected_OG), phenoData = AnnotatedDataFrame(pData(tsalik_OG)) )
saveRDS(TSALIK_ORIGINAL, "/Users/kkalantar/Documents/Research/NMLHC/geo_data/TSALIK_ORIGINAL.rds")





# Each infectious agent represents a unique combination of pathogen-associated molecular patterns that interact with specific 
# pattern-recognition receptors expressed on immune cells. Therefore, we surmised that the blood immune cells of individuals with 
# different infections might bear discriminative transcriptional signatures. Gene expression profiles were obtained for 131 
# peripheral blood samples from pediatric patients with acute infections caused by influenza A virus, Gram-negative 
# (Escherichia coli) or Gram-positive (Staphylococcus aureus and Streptococcus pneumoniae) bacteria. 

GSE6269_GPL96 <- GEOquery::getGEO(filename="/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE6269-GPL96_series_matrix.txt") 
GSE6269_GPL96_gn <- detect_genecnt_platform(exprs(GSE6269_GPL96))
Tsalik_GSE6269_GPL96 <- ExpressionSet(assayData = as.matrix(GSE6269_GPL96_gn), phenoData = AnnotatedDataFrame(pData(GSE6269_GPL96)))
saveRDS(Tsalik_GSE6269_GPL96,"/Users/kkalantar/Documents/Research/NMLHC/geo_data/Tsalik_GSE6269_GPL96.rds")

# I'm not planning to use the validation sets from this GEO dataset
#g2 <- GEOquery::getGEO(filename="/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE6269-GPL570_series_matrix.txt")    
#g3 <- GEOquery::getGEO(filename="/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE6269-GPL2507_series_matrix.txt”)    
#c1 <- Biobase::combine(g, g2)
#c2 <- Biobase::combine(c1, g3)


# variable: characteristics_ch1.4; contrast Adenovirus, HHV6, Enterovirus, Rhinovirus, E.coli, NSSA, None
GSE40396 <- GEOquery::getGEO(filename="/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE40396_series_matrix.txt")
GSE40396_gn <- detect_genecnt_platform(exprs(GSE40396))
Tsalik_GSE40396 <- ExpressionSet(assayData = as.matrix(GSE40396_gn), phenoData = AnnotatedDataFrame(pData(GSE40396)))
saveRDS(Tsalik_GSE40396, "/Users/kkalantar/Documents/Research/NMLHC/geo_data/Tsalik_GSE40396.rds")


# variable: characteristics_ch1.1; contrast bacterial, none, Influenza, RSV
GSE42026 <- GEOquery::getGEO(filename = "/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE42026_series_matrix.txt")
GSE42026_gn <- detect_genecnt_platform(exprs(GSE42026))
Tsalik_GSE42026 <- ExpressionSet(assayData = as.matrix(GSE42026_gn), phenoData = AnnotatedDataFrame(pData(GSE42026)))
saveRDS(Tsalik_GSE42026, "/Users/kkalantar/Documents/Research/NMLHC/geo_data/Tsalik_GSE42026.rds")


# variable: source_name_ch1; bacterial, Severe, Vaccine
GSE20346 <-  GEOquery::getGEO(filename = "/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE20346_series_matrix.txt")
GSE20346_gn <- detect_genecnt_platform(exprs(GSE20346))
Tsalik_GSE20346 <- ExpressionSet(assayData = as.matrix(GSE20346_gn), phenoData = AnnotatedDataFrame(pData(GSE20346)))
saveRDS(Tsalik_GSE20346, "/Users/kkalantar/Documents/Research/NMLHC/geo_data/Tsalik_GSE20346.rds")


# variable: characteristics_ch1.2; lung cancer, Control, TB, Active Sarcoid, Non-active sarcoidosis
GSE42834 <- GEOquery::getGEO(filename = "/Users/kkalantar/Documents/Research/NMLHC/geo_data/GSE42834_series_matrix.txt")
GSE42834_gn <- detect_genecnt_platform(exprs(GSE42834))
Tsalik_GSE42834 <- ExpressionSet(assayData = as.matrix(GSE42834_gn), phenoData = AnnotatedDataFrame(pData(GSE42834)))
saveRDS(Tsalik_GSE42834, "/Users/kkalantar/Documents/Research/NMLHC/geo_data/Tsalik_GSE42834.rds")


### NOTE - I might need a way to say "contrast A vs. all others" when I generate the datasets in the script (id for GSE42834 it would be TB v. all others)




