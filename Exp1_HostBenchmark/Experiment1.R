setwd("/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark")
source("/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/HostExpnFunctions.R")

geo_data <- GEOquery::getGEO(filename='/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE60244_series_matrix.txt')
dataset <- subset_known_dataset(geo_data, "BACTERIA", "VIRUS")
project_pca(dataset,c("red","blue"))
result <- run_scramble_classification(dataset, 5, 1000, .2)
plot_res_heatmaps(result)

g <- GEOquery::getGEO(filename='/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE33341-GPL1261_series_matrix.txt')
dataset <- subset_known_dataset(g, "S. aureus", "E. coli")
project_pca(dataset,c("red","blue"))
result <- run_scramble_classification(dataset, 5, 1000, .2)
plot_res_heatmaps(result)

data <- filtered_eset  # from mBALPkg
TRAINING_NAMES <- c("TA.212","TA.225","TA.298","TA.304","TA.314","TA.315","TA.335","TA.337","TA.343","TA.350", 
                    "TA.349","TA.273","TA.331","TA.221","TA.220","TA.215","TA.270","TA.241","TA.211","TA.218")  # same as in mBAL study
DEgenes <- read.table("Documents/Research/NMLHC/Exp1_HostBenchmark/data/DEgenes.csv")  # taken from mBAL study
x <- generate_simulated_dataset(data, TRAINING_NAMES, DEgenes)
td <- data[,TRAINING_NAMES]
td2 <- t(t(exprs(td))/colSums(exprs(td)))
exprs(td) <- td2

project_pca(x,c("red","blue"))
project_pca2(td, x,c("red","blue"))
project_pca3(td, x,c("red","blue"))
result <- run_scramble_classification(x, 5, 1000, .2)
plot_res_heatmaps(result)


library(scales)

project_pca <- function(input_dataset, c){
  p <- prcomp(t(exprs(input_dataset$train_data)))
  #print(alpha(c[as.integer(input_dataset$train_labels == 1) + 1], .5))
  plot(p$x[,1],p$x[,2], 
       col=c[as.integer(input_dataset$train_labels == 1) + 1], 
       xlab="PC1",ylab="PC2", cex=1.6, pch=16)
  
  s <- scale(t(exprs(input_dataset$test_data)), p$center, p$scale) %*% p$rotation 
  points(s[,1],s[,2], col=c[as.integer(input_dataset$test_labels == 1)  + 1], cex=1.6)
}



project_pca2 <- function(td, input_dataset, c){
  p <- prcomp(t(exprs(td))) #t(exprs(input_dataset$train_data)))
  plot(p$x[,1],p$x[,2],
       col=c[as.integer(td$effective_group == 1) + 1],
       xlab="PC1",ylab="PC2", cex=1.6, pch=16)

  s2 <- scale(t(exprs(input_dataset$train_data)), p$center, p$scale) %*% p$rotation
  points(s2[,1],s2[,2], col=c[as.integer(input_dataset$train_labels == 1) + 1], cex = 1.6, pch = 4)
  print(s2[,1:2])


  s <- scale(t(exprs(input_dataset$test_data)), p$center, p$scale) %*% p$rotation
  points(s[,1],s[,2], col=c[as.integer(input_dataset$test_labels == 1)  + 1], cex=1.6)

}



project_pca3 <- function(td, input_dataset, c){
  
  p <- prcomp(t(exprs(input_dataset$train_data)))
  plot(p$x[,1],p$x[,2],
       col=c[as.integer(input_dataset$train_labels == 1) + 1],
       xlab="PC1",ylab="PC2", cex=1.6, pch=16)
  
  s <- scale(t(exprs(input_dataset$test_data)), p$center, p$scale) %*% p$rotation
  points(s[,1],s[,2], col=c[as.integer(input_dataset$test_labels == 1)  + 1], cex=1.6)
  
  s2 <- scale(t(exprs(td)), p$center, p$scale) %*% p$rotation
  points(s2[,1],s2[,2], col=c[as.integer(td$effective_group == 1) + 1], cex = 1.6, pch = 4)
  print(s2[,1:2])
  
  
}
project_pca3(td, x,c("red","blue"))


# MAKE A FUNCTION THAT DOES THIS
data <- filtered_eset  # from mBALPkg
TRAINING_NAMES <- c("TA.212","TA.225","TA.298","TA.304","TA.314","TA.315","TA.335","TA.337","TA.343","TA.350", 
                    "TA.349","TA.273","TA.331","TA.221","TA.220","TA.215","TA.270","TA.241","TA.211","TA.218")  # same as in mBAL study
DEgenes <- read.table("Documents/Research/NMLHC/Exp1_HostBenchmark/data/DEgenes.csv")  # taken from mBAL study
x <- generate_simulated_dataset(data, TRAINING_NAMES, DEgenes)