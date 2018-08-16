setwd("/Users/kkalantar/Documents/Research/Biomarker/Exp1_HostBenchmark")
source("HostExpnFunctions.R")

geo_data <- GEOquery::getGEO(filename='/Users/kkalantar/Documents/Research/Biomarker/Exp1_HostBenchmark/GSE60244_series_matrix.txt')
dataset <- subset_known_dataset(geo_data, "BACTERIA", "VIRUS")
project_pca(dataset,c("red","blue"))
result <- run_scramble_classification(dataset, 5, 1000, .2)
plot_res_heatmaps(result)

g <- GEOquery::getGEO(filename='/Users/kkalantar/Documents/Research/Biomarker/Exp1_HostBenchmark/GSE33341-GPL1261_series_matrix.txt')
dataset <- subset_known_dataset(g, "S. aureus", "E. coli")
project_pca(dataset,c("red","blue"))
result <- run_scramble_classification(dataset, 5, 1000, .2)
plot_res_heatmaps(result)

data <- filtered_eset  # from mBALPkg
TRAINING_NAMES <- c("TA.212","TA.225","TA.298","TA.304","TA.314","TA.315","TA.335","TA.337","TA.343","TA.350", 
                    "TA.349","TA.273","TA.331","TA.221","TA.220","TA.215","TA.270","TA.241","TA.211","TA.218")  # same as in mBAL study
DEgenes <- read.table("Documents/Research/Biomarker/DEgenes.csv")  # taken from mBAL study
x <- generate_simulated_dataset(data, TRAINING_NAMES, DEgenes)
project_pca(x,c("red","blue"))
result <- run_scramble_classification(x, 5, 1000, .2)
plot_res_heatmaps(result)




