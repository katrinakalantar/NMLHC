#source("http://bioconductor.org/biocLite.R")
library(R.matlab)      # enables cross-talk with matlab server
library(PopED)         # contains feval function
library(matlab)        # contains repmat function equivalent to matlab version
library(pROC)          # contains roc function
library(glmnet)        # contains glmnet function for regularized regression in standard R library
library(gplots)        # contains heatmap.2 function
library(jsonlite)      # enables reading .json input parameters
library(pheatmap)      # enables plotting heatmap from pheatmap package
library(grid)          # plot multiple heatmaps on one page
library(gridExtra)     # plot multiple heatmaps on one page
library(ggplot2)       # contains the "theme()" function used for controlling margins of multiple pheatmaps
library(RColorBrewer)  # contains color palettes used in line plots
library(purrr)         # contains the possibly() function to enable continued execution after mablab error
library(permute)       # contains the shuffle() function used in nmlhc_matlab2R_functions.R
library(reshape2)      # contains the melt() function used in ggplot 
library(mBALPkg)       # contains the data required to simulate data from mini-BAL cohort
library(SimSeq)        # package required for simulation from an existing dataset
library(SuperLearner)
library(caret)
library(ROCR)
library(cvAUC)
librrary(randomForest)
library(NoiseFiltersR)
library(biomaRt)       # required for converting genenames
library(GSA)           # required for reading the .gmt file for collapsing gene names
library(statmod)
library(parallel)


source("/Users/kkalantar/Documents/Research/NMLHC/nmlhc_matlab2R_functions.R")
source("/Users/kkalantar/Documents/Research/NMLHC/microarray_format_utils.R")
source("/Users/kkalantar/Documents/Research/NMLHC/pylogger.R")  #sourced from: https://gist.github.com/jonathancallahan/3ed51265d3c6d56818458de95567d3ae

pca_cols <- c("turquoise","green","magenta")
pca_cols2 <- c("turquoise3","green3","magenta3")

#
# set the experiment directory - this will contain the input parameter file and is where all the output will be written
#

EXPERIMENT_DIR <- "/Users/kkalantar/Documents/Research/NMLHC/EXPERIMENTS/Experiment_1/"
parameters <- read_json(paste(EXPERIMENT_DIR, "parameters.json", sep=""), simplifyVector = TRUE)
logger.setup(infoLog = paste(EXPERIMENT_DIR, "info.txt"))
init_log(parameters)

#
# set global options, quit script if required parameters were not provided
# 

stopifnot(check_params(parameters, c("ITER", "CLS", "DATASET_PARAM", 
                                     "use_PCs","feature_select",
                                     "common_reg","common_sn","common_maxIter",
                                     "DS_SIZE_RANGE","DIM_RANGE","EXP_RANGE","EXP_RANGE_J")))

ITER = parameters$ITER                              
CLS = parameters$CLS                               # number of classes
DATASET_PARAM = parameters$DATASET_PARAM           # data generation method 
use_PCs = parameters$use_PCs
feature_select = parameters$feature_select 
common_reg = parameters$common_reg                 # regularization type
common_sn = parameters$common_sn      
common_maxIter = parameters$common_maxIter         # max iterations for the algorithm

# load all GEO datasets up front (bc this is slow) 
list_of_geo_datasets <- list()     
for(dat in rownames(parameters$datasets)){
  d <- parameters$datasets[dat,]
  if(d$type == "geo"){
    if(grepl("series", d$series_filename)){
      list_of_geo_datasets[[d$name]] <- GEOquery::getGEO(filename=d$series_filename)      
    }else if(grepl("rds", d$series_filename)){
      list_of_geo_datasets[[d$name]] <- readRDS(d$series_filename)
    }
  }
}

#list_of_geo_datasets[[1]] <- list_of_geo_datasets[[1]][,list_of_geo_datasets[[1]]$`infection:ch1` == "Non-infected"]


DS_SIZE_RANGE   = parameters$DS_SIZE_RANGE 
DIM_RANGE       = parameters$DIM_RANGE
EXP_RANGE       = parameters$EXP_RANGE
EXP_RANGE_J     = parameters$EXP_RANGE_J

# preallocating error storage
err_lr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
err_rlr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
err_gammalr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);

# preallocating AUC storage
auc_lr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
auc_rlr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
auc_gammalr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);

CMs_MF = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
CMs_MFauc = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);
CMs_CF = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, parameters$datasets$name, EXP_RANGE_J, EXP_RANGE);

features_selected <- list() 


#
# Run analysis with given parameters, track results
#

for(dimr in seq(1:length(DIM_RANGE))){
  for(dsize in seq(1:length(DS_SIZE_RANGE))){
    for(k in seq(1:ITER)){           # number of iterations to run (to obtain mean performance)
      
      # BRAINSTORNING - create list of datasets here (include gendataset bc all relevant params have been identified, 
      #                 draw train/test sets which will remain the same for all iterations of flipping...seems OK to me)
      DIM = DIM_RANGE[dimr]  # n_features
      DS_SIZE = DS_SIZE_RANGE[dsize]
      
      for(d in seq(1:1)){

        for(j in seq(1:length(EXP_RANGE_J))){
          for(i in seq(1:length(EXP_RANGE))){
            
            flip_j = EXP_RANGE_J[j]
            flip_i = EXP_RANGE[i]
            
            dataset <- feval('subset_geo_cv', parameters$datasets[1,]$name, list_of_geo_datasets, parameters$datasets[1,]$source_variable, flip_i, flip_j)

            logger.info(msg = paste(c("MAIN - ITER - ", "DATA = ", d, ", K = ", k, ", flip_j = ", flip_j, ", flip_i = ", flip_i), collapse="" ))

            
            cv_results <- list()
            
            for(cv in seq(1:length(dataset))){
              # MOVE ALL THE ACTION INTO HERE!!! 
              
              Xt = dataset[[cv]]$x
              Xt = Xt[,colSums(Xt) != 0]  # select only the genes with > 0 reads
              yt = dataset[[cv]]$y
              
              Xs = dataset[[cv]]$xx
              Xs = Xs[,colSums(Xs) != 0]  # select only the genes with > 0 reads
              ys = dataset[[cv]]$tt
              
              print(dim(Xt))
              print(dim(Xs))
              
              # a priori feature selection - http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
              #features_to_keep <- apriori_feature_selection(Xt)
              #feature_names <- intersect(intersect(colnames(Xt), colnames(Xs)), features_to_keep)
              
              # use only features with > 0 reads in training set and test set
              feature_names <- intersect(colnames(Xt), colnames(Xs))
              Xs = Xs[,feature_names]      # select only the genes that were present in the training set 
              Xt = Xt[,feature_names]
              
              logger.info(msg = sprintf("MAIN - DATA - TRAIN - N = %i samples; G1 = %i (%f), G2 = %i (%f)", length(yt), table(yt)[1], table(yt)[1]/length(yt), table(yt)[2], table(yt)[2]/length(yt)))
              logger.info(msg = sprintf("MAIN - DATA - TEST - N = %i samples; G1 = %i (%f), G2 = %i (%f)", length(ys), table(ys)[1], table(ys)[1]/length(ys), table(ys)[2], table(ys)[2]/length(ys)))
              
              Xt = scale(Xt)  # CHECK ON THIS
              Xs = scale(Xs)  # CHECK ON THIS
              
              yz = dataset[[cv]]$yz
              fdz = dataset[[cv]]$ff
              
              pca_res <- prcomp(Xt)
              Xt_pcatrans <- pca_res$x
              
              # filter with PC
              print("FIRST FILTRATION:")
              filter <- make_ensemble_parallel(Xt_pcatrans, y = unlist(lapply(yz, function(x){if(x==1){return("one")}else if(x==2){return("two")}})), 
                                               c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"

              
              # in this plot, the larger points are the ones that get kept.
              plot(pca_res$x[,1],pca_res$x[,2], col = pca_cols[yt] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz < 0) + 1],
                   cex = c(1.0, 2.2)[as.integer(filter$MF) + 1], main = "Iter1")
              
              a <- flag_flipped_samples(filter$MF, fdz)
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - %s", toString(a$table)))
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - accuracy = %s", toString(a$overall["Accuracy"])))
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - AUC = %s", toString(pROC::roc((fdz==1),rowSums(filter$full_res))$auc)))
              
              idx_keep <- which(filter$MF)
              Xt_confident <- Xt[idx_keep,]
              yz_confident <- yz[idx_keep]
              yt_confident <- yt[idx_keep]
              fdz_confident <- fdz[idx_keep]
              
              ####### LOOPING AGAIN!!! #######
              print("SECOND FILTRATION:")
              pca_res2 <- prcomp(Xt_confident)
              Xt_pcatrans2 <- pca_res2$x
              filter2 <- make_ensemble_parallel(Xt_pcatrans2, y = unlist(lapply(yz_confident, function(x){if(x==1){return("one")}else if(x==2){return("two")}})),
                                               c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"
              
              # in this plot, the larger points are the ones that get kept.
              plot(pca_res2$x[,1],pca_res2$x[,2], col = pca_cols2[yt_confident] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz_confident < 0) + 1],
                   cex = c(1.0, 2.2)[as.integer(filter2$MF) + 1], main = "Iter2")
              #text(pca_res2$x[,1] + 2,pca_res2$x[,2] + 1, rowSums(filter2$full_res))
              a2 <- flag_flipped_samples(filter2$MF, fdz_confident)
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - %s", toString(a2$table)))
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - accuracy = %s", toString(a2$overall["Accuracy"])))
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - AUC = %s", toString(pROC::roc((fdz_confident==1),rowSums(filter2$full_res))$auc)))

              idx_keep <- which(filter2$MF)
              Xt_confident <- Xt_confident[idx_keep,]
              yz_confident <- yz_confident[idx_keep]
              yt_confident <- yt_confident[idx_keep]
              fdz_confident <- fdz_confident[idx_keep]

              ####### ####### ####### #######
              
              # clean model
              CV <- cv.glmnet(Xt_confident, y = (yz_confident == 2), alpha = .2)
              model_clean <- glmnet(Xt_confident, y = (yz_confident == 2), family = "binomial", lambda = CV$lambda.min, alpha = 1)
              res_train <- predict(model_clean, Xt_confident, type = "response")
              roc_train <- pROC::roc((yz_confident == 2), res_train[,1])
              logger.info(msg = paste("CV", cv, "Train AUC: ", roc_train$auc))
              res_test <- predict(model_clean, Xs, type = "response")
              roc_test <- pROC::roc((ys==2), res_test[,1])
              auc_rlr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
              logger.info(msg = paste("CV", cv, "Test AUC: ", roc_test$auc))
              
              cv_results[[cv]] <- res_test[,1]
              
            }
            
            final_ROC <- pROC::roc(unlist(lapply(dataset, function(x){return(x$tt)}))==2, unlist(cv_results))
            logger.info(msg = paste("Final ROC:", final_ROC$auc))
            
            newdat <- rbind(Xt,Xs)
            newpca <- prcomp(newdat)
            newpca_trans <- newpca$x
            #set.seed(20)
            clusters <- kmeans(newpca_trans, 2)
            par(mfrow = c(1,4))
            plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[clusters$cluster], main = "k-means cluster prediction")                       # cluster result scatterplot
            plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[c(dataset[[10]]$y, dataset[[10]]$tt)], main = "truth")   # truth scatterplot
            plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[c(yz, dataset[[10]]$tz)], main = "mislabelled",
                 pch = c(16,1)[as.integer(as.factor(c(yz, dataset[[10]]$tz)) == as.factor(c(dataset[[10]]$y, dataset[[10]]$tt))) + 1])                # mislabelled scatterplot
            
            colfunc <- colorRampPalette(c("red", "blue"))
            plot(newpca_trans[,1], newpca_trans[,2], col=c("red","blue")[round(unlist(cv_results)[names(newpca_trans[,1])]) + 1],
                 main = "KK_alg prediction", pch = c(16,1)[as.integer(as.factor(round(unlist(cv_results)[names(newpca_trans[,1])]) + 1) == as.factor(c(dataset[[10]]$y, dataset[[10]]$tt))) + 1])       # algorithm-predicted scatterplot
            barplot(unlist(cv_results), col = c("red","blue")[(unlist(lapply(dataset, function(x){return(x$tt)}))==2) + 1])             # barplot
            barplot(unlist(clusters$cluster)[names(unlist(cv_results))], col = c("red","blue")[(unlist(lapply(dataset, function(x){return(x$tt)}))==2) + 1])            
            
            
            cm_final_kmeans <- confusionMatrix(as.factor(clusters$cluster), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            logger.info(msg = paste("k-means achieves accuracy:", cm_final_kmeans))
            cm_final_kkalg <- confusionMatrix(as.factor(round(unlist(cv_results)[names(newpca_trans[,1])]) + 1), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            logger.info(msg = paste("kk-algorithm achieves accuracy:", cm_final_kkalg))
            
            cm_final_mislabelling <- confusionMatrix(as.factor(c(yz, dataset[[10]]$tz)), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            logger.info(msg = paste("final_mislabelling achieves accuracy:", cm_final_mislabelling))
            
            
            
            
            
            # cv_results2 <- list()
            # 
            # for(cv in seq(1:length(dataset))){
            #   # MOVE ALL THE ACTION INTO HERE!!!
            #   
            #   logger.info(msg="IN SECOND ROUND OF CLASSIFICATION")
            #   
            #   Xt = dataset[[cv]]$x
            #   Xt = Xt[,colSums(Xt) != 0]  # select only the genes with > 0 reads
            #   yt = dataset[[cv]]$y
            #   
            #   Xs = dataset[[cv]]$xx
            #   Xs = Xs[,colSums(Xs) != 0]  # select only the genes with > 0 reads
            #   ys = dataset[[cv]]$tt
            #   
            #   print(dim(Xt))
            #   print(dim(Xs))
            #   
            #   # a priori feature selection - http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
            #   #features_to_keep <- apriori_feature_selection(Xt)
            #   #feature_names <- intersect(intersect(colnames(Xt), colnames(Xs)), features_to_keep)
            #   
            #   # use only features with > 0 reads in training set and test set
            #   feature_names <- intersect(colnames(Xt), colnames(Xs))
            #   Xs = Xs[,feature_names]      # select only the genes that were present in the training set 
            #   Xt = Xt[,feature_names]
            #   
            #   Xt = scale(Xt)  # CHECK ON THIS
            #   Xs = scale(Xs)  # CHECK ON THIS
            #   
            #   yz = as.integer(unlist(cv_results)[rownames(Xt)]>.5)+1
            #   fdz = (as.integer(yt == (as.integer(unlist(cv_results)[rownames(Xt)]>.5)+1)) + 1)   # ?????? 
            #   
            #   pca_res <- prcomp(Xt)
            #   Xt_pcatrans <- pca_res$x
            #   
            #   # filter with PC
            #   print("FIRST FILTRATION:")
            #   #filter <- make_ensemble_parallel(Xt_pcatrans, y = unlist(lapply(yz, function(x){if(x==1){return("one")}else if(x==2){return("two")}})), 
            #   #                                 c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"
            #   filter <- make_ensemble(Xt_pcatrans, y = unlist(lapply(yz, function(x){if(x==1){return("one")}else if(x==2){return("two")}})), 
            #                           c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"
            #   
            #   # in this plot, the larger points are the ones that get kept.
            #   plot(pca_res$x[,1],pca_res$x[,2], col = pca_cols[yt] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz)],
            #        cex = c(1.0, 2.2)[as.integer(filter$MF) + 1], main = "Iter1")
            #   
            #   #a <- flag_flipped_samples(filter$MF, fdz)
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - %s", toString(a$table)))
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - accuracy = %s", toString(a$overall["Accuracy"])))
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - AUC = %s", toString(pROC::roc((fdz==1),rowSums(filter$full_res))$auc)))
            #   
            #   idx_keep <- which(filter$MF)
            #   Xt_confident <- Xt[idx_keep,]
            #   yz_confident <- yz[idx_keep]
            #   yt_confident <- yt[idx_keep]
            #   fdz_confident <- fdz[idx_keep]
            #   
            #   ####### LOOPING AGAIN!!! #######
            #   print("SECOND FILTRATION:")
            #   pca_res2 <- prcomp(Xt_confident)
            #   Xt_pcatrans2 <- pca_res2$x
            #   #filter2 <- make_ensemble_parallel(Xt_pcatrans2, y = unlist(lapply(yz_confident, function(x){if(x==1){return("one")}else if(x==2){return("two")}})),
            #   #                                 c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"
            #   filter2 <- make_ensemble_parallel(Xt_pcatrans2, y = unlist(lapply(yz_confident, function(x){if(x==1){return("one")}else if(x==2){return("two")}})),
            #                                     c("rf","svmLinear3","regLogistic","knn","nnet"),multiple = FALSE)  #,"nnet","nnet",, "mlp","rFerns" "RFlda",hdda  , ,"nnet" "knn",,"nnet", ,"nnet","regLogistic"
            #   # in this plot, the larger points are the ones that get kept.
            #   plot(pca_res2$x[,1],pca_res2$x[,2], col = pca_cols2[yt_confident] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz_confident)],
            #        cex = c(1.0, 2.2)[as.integer(filter2$MF) + 1], main = "Iter2")
            #   #text(pca_res2$x[,1] + 2,pca_res2$x[,2] + 1, rowSums(filter2$full_res))
            #   #a2 <- flag_flipped_samples(filter2$MF, fdz_confident)
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - %s", toString(a2$table)))
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - accuracy = %s", toString(a2$overall["Accuracy"])))
            #   #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF2 - AUC = %s", toString(pROC::roc((fdz_confident==1),rowSums(filter2$full_res))$auc)))
            #   
            #   idx_keep <- which(filter2$MF)
            #   Xt_confident <- Xt_confident[idx_keep,]
            #   yz_confident <- yz_confident[idx_keep]
            #   yt_confident <- yt_confident[idx_keep]
            #   #fdz_confident <- fdz_confident[idx_keep]
            #   
            #   ####### ####### ####### #######
            #   
            #   # clean model
            #   CV <- cv.glmnet(Xt_confident, y = (yz_confident == 2), alpha = .2)
            #   model_clean <- glmnet(Xt_confident, y = (yz_confident == 2), family = "binomial", lambda = CV$lambda.min, alpha = 1)
            #   res_train <- predict(model_clean, Xt_confident, type = "response")
            #   roc_train <- pROC::roc((yz_confident == 2), res_train[,1])
            #   logger.info(msg = paste("CV", cv, "Train AUC: ", roc_train$auc))
            #   res_test <- predict(model_clean, Xs, type = "response")
            #   roc_test <- pROC::roc((ys==2), res_test[,1])
            #   auc_rlr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
            #   logger.info(msg = paste("CV", cv, "Test AUC: ", roc_test$auc))
            #   
            #   cv_results2[[cv]] <- res_test[,1]
            #   
            # }
            # 
            # #final_ROC2 <- pROC::roc(unlist(lapply(dataset, function(x){return(x$tt)}))==2, unlist(cv_results2))
            # final_ROC2 <- pROC::roc(unlist(lapply(dataset, function(x){return(x$tt)}))==2, unlist(cv_results2)[unlist(lapply(dataset, function(x){return(rownames(x$xx))}))])
            # logger.info(msg = paste("Final ROC2:", final_ROC2$auc))
            # 
            # newdat <- rbind(Xt,Xs)
            # newpca <- prcomp(newdat)
            # newpca_trans <- newpca$x
            # #set.seed(20)
            # clusters <- kmeans(newpca_trans, 2)
            # par(mfrow = c(1,4))
            # plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[clusters$cluster], main = "k-means cluster prediction")                       # cluster result scatterplot
            # plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[c(dataset[[10]]$y, dataset[[10]]$tt)], main = "truth")   # truth scatterplot
            # plot(newpca_trans[,1], newpca_trans[,2], col= c("red","blue")[c(yz, dataset[[10]]$tz)], main = "mislabelled",
            #      pch = c(16,1)[as.integer(as.factor(c(yz, dataset[[10]]$tz)) == as.factor(c(dataset[[10]]$y, dataset[[10]]$tt))) + 1])                # mislabelled scatterplot
            # 
            # colfunc <- colorRampPalette(c("red", "blue"))
            # plot(newpca_trans[,1], newpca_trans[,2], col=c("red","blue")[round(unlist(cv_results2)[names(newpca_trans[,1])]) + 1],
            #      main = "KK_alg prediction", pch = c(16,1)[as.integer(as.factor(round(unlist(cv_results2)[names(newpca_trans[,1])]) + 1) == as.factor(c(dataset[[10]]$y, dataset[[10]]$tt))) + 1])       # algorithm-predicted scatterplot
            # barplot(unlist(cv_results2), col = c("red","blue")[(unlist(lapply(dataset, function(x){return(x$tt)}))==2) + 1])             # barplot
            # barplot(unlist(clusters$cluster)[names(unlist(cv_results2))], col = c("red","blue")[(unlist(lapply(dataset, function(x){return(x$tt)}))==2) + 1])            
            # 
            # 
            # cm_final_kmeans2 <- confusionMatrix(as.factor(clusters$cluster), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            # logger.info(msg = paste("k-means2 achieves accuracy:", cm_final_kmeans2))
            # cm_final_kkalg2 <- confusionMatrix(as.factor(round(unlist(cv_results2)[names(newpca_trans[,1])]) + 1), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            # logger.info(msg = paste("kk-algorithm2 achieves accuracy:", cm_final_kkalg2))
            # cm_final_mislabelling2 <- confusionMatrix(as.factor(c(yz, dataset[[10]]$tz)), as.factor(c(dataset[[10]]$y, dataset[[10]]$tt)), positive = "2")
            # logger.info(msg = paste("final_mislabelling achieves accuracy:", cm_final_kkalg2))
            # 
            
            
            
            
            
            
            
          }
        }
      }
    }
  }
}


save_data(auc_lr, get_variable_name(auc_lr), EXPERIMENT_DIR)
save_data(auc_rlr, get_variable_name(auc_rlr), EXPERIMENT_DIR)
save_data(err_lr, get_variable_name(err_lr), EXPERIMENT_DIR)
save_data(err_rlr, get_variable_name(err_rlr), EXPERIMENT_DIR)


#
# Plotting Functions
#
# 
my_palette <- colorRampPalette(c("purple3", "white"))(n = 1001)
my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)    # continuous color scale for heatmaps
discrete_pal <- brewer.pal(6, "Spectral")                          # color palette for discrete variables - line plots for different conditions




# P-value rLR v. gammaLR
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr[1,1,,,,][,x], auc_gammalr[1,1,,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr[1,1,,,,][,x], auc_lr[1,1,,,,][,x], alternative = "greater")$p.value)}))


unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,1,,,][,x], auc_gammalr_test[1,1,1,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,1,,,][,x], auc_lr_test[1,1,1,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,2,,,][,x], auc_gammalr_test[1,1,2,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,2,,,][,x], auc_lr_test[1,1,2,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,3,,,][,x], auc_gammalr_test[1,1,3,,,][,x], alternative = "greater")$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,3,,,][,x], auc_lr_test[1,1,3,,,][,x], alternative = "greater")$p.value)}))


unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,1,,,][,x], auc_gammalr_test[1,1,1,,,][,x])$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,1,,,][,x], auc_lr_test[1,1,1,,,][,x])$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,2,,,][,x], auc_gammalr_test[1,1,2,,,][,x])$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,2,,,][,x], auc_lr_test[1,1,2,,,][,x])$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,3,,,][,x], auc_gammalr_test[1,1,3,,,][,x])$p.value)}))
unlist(lapply(seq(1:4),function(x){return(wilcox.test(auc_rlr_test[1,1,3,,,][,x], auc_lr_test[1,1,3,,,][,x])$p.value)}))



plot_confusion_matrix(CMs_MF, colorRampPalette(c("white","dodgerblue"))(n = 401), max=400)
#plot_confusion_matrix(CMs_CF, colorRampPalette(c("white","purple"))(n = 301), max = 300)


plot_params <- list("DS_SIZE_RANGE" = 2000, "DIM_RANGE" = 200, "DATASET" = "DATASETS_viral_bacterial", "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)



colMeans(auc_gammalr_test[1,1,1,,,] - auc_rlr_test[1,1,1,,,])
colMeans(auc_lr_test[1,1,1,,,] - auc_rlr_test[1,1,1,,,])


# TEST SET
colors <- brewer.pal(5, "Dark2")
matplot(t(auc_rlr[1,1,1,,,] - auc_gammalr[1,1,1,,,]), col = colors, cex = 1.8, pch=0, lwd=2,
        main = "rLR - gammaLR", xlab = "J_FLIP", ylab = "Difference", ylim = c(-.5, .5))
abline(h = 0, lty = 2)
matplot(t(auc_rlr[1,1,1,,,] - auc_lr[1,1,1,,,]), col = colors, cex = 1.8, pch=0, lwd=2,
        main = "rLR - LR", xlab = "J_FLIP", ylab = "Difference", ylim = c(-.5, .5))
abline(h = 0, lty = 2)

# VALIDATION SET(S)
matplot(t(auc_rlr_test[1,1,3,,,] - auc_gammalr_test[1,1,3,,,]), col = colors, cex = 1.8, pch=0, lwd=2,
        main = "rLR - gammaLR", xlab = "J_FLIP", ylab = "Difference", ylim = c(-.5, .5))
matplot(t(auc_rlr_test[1,1,2,,,] - auc_gammalr_test[1,1,2,,,]), col = colors, cex = 1.8, pch=1, lwd=2, add = TRUE)
matplot(t(auc_rlr_test[1,1,1,,,] - auc_gammalr_test[1,1,1,,,]), col = colors, cex = 1.8, pch=2, lwd=2, add = TRUE)
abline(h = 0, lty = 2)

matplot(t(auc_rlr_test[1,1,3,,,] - auc_lr_test[1,1,3,,,]), col = colors, cex = 1.8, pch=0, lwd=2,
        main = "rLR - LR", xlab = "J_FLIP", ylab = "Difference", ylim = c(-.5, .5))
matplot(t(auc_rlr_test[1,1,2,,,] - auc_lr_test[1,1,2,,,]), col = colors, cex = 1.8, pch=1, lwd=2, add = TRUE)
matplot(t(auc_rlr_test[1,1,1,,,] - auc_lr_test[1,1,1,,,]), col = colors, cex = 1.8, pch=2, lwd=2, add = TRUE)
abline(h = 0, lty = 2)



