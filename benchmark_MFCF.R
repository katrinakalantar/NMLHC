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



source("/Users/kkalantar/Documents/Research/NMLHC/nmlhc_matlab2R_functions.R")
source("/Users/kkalantar/Documents/Research/NMLHC/microarray_format_utils.R")
source("/Users/kkalantar/Documents/Research/NMLHC/pylogger.R")  #sourced from: https://gist.github.com/jonathancallahan/3ed51265d3c6d56818458de95567d3ae

pca_cols <- c("turquoise","green","magenta")

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
for(dat in rownames(parameters$testdata)){
  d <- parameters$testdata[dat,]
  if(d$type == "geo"){
    if(grepl("series", d$series_filename)){
      list_of_geo_datasets[[d$name]] <- GEOquery::getGEO(filename=d$series_filename)      
    }else if(grepl("rds", d$series_filename)){
      list_of_geo_datasets[[d$name]] <- readRDS(d$series_filename)
    }
  }
}

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
      
      list_of_datasets <- list()
      for(x in rownames(parameters$datasets)){
        d <- parameters$datasets[x,]
        if(d$type == "gen"){
          list_of_datasets[[d$name]] <- feval("generate_data", as.numeric(CLS), as.numeric(DIM), as.numeric(DS_SIZE), 1000, as.numeric(d$class_sep), "gen")
          logger.info(msg=paste(c("DATA - GEN - ", "CLS: ", CLS, ", DIM: ", DIM, ", DS_SIZE: ", DS_SIZE, ", CLASS_SEP: ", d$class_sep), collapse=""))
        }else if(d$type == "geo"){
          list_of_datasets[[d$name]] <- feval('subset_geo', d$name, list_of_geo_datasets, d$source_variable)
          logger.info(msg=paste("DATA - GEO - ", d$name, sep=""))
        }else if(d$type == "sim"){
          list_of_datasets[[d$name]] <- feval('simulate_data', d$base_dataset, d$n_samples)
          logger.info(msg=paste("DATA - SIM - ", d$name, sep=""))
        }else{
          list_of_datasets[[d$name]] <- NULL
        }
      }
      
      
      # 
      # Set up test datasets 
      #
      list_of_test_sets <- list()
      for(x in rownames(parameters$testdata)){
        d <- parameters$testdata[x,]
        if(d$type == "geo"){
          list_of_test_sets[[d$name]] <- feval('subset_geo', d$name, list_of_geo_datasets, d$source_variable)
          logger.info(msg=paste("DATA - GEO - ", d$name, sep=""))
        }
      }
      
      
      for(d in seq(1:length(list_of_datasets))){
        
        dataset <- list_of_datasets[[names(list_of_datasets)[d]]]
        
        # if(parameters$collapse_pathways){
        #   g <- GSA.read.gmt(parameters$geneset_file)    # this is the .gmt file used for any pathway analysis
        #   dataset$x <- convert_genenames(dataset$x)     # convert the gene names prior to collapsing pathways, even if that variable isn't set in parameters.json
        #   dataset$xx <- convert_genenames(dataset$xx)
        #   keep_pathways <- unique(unlist(lapply(parameters$pathway_keywords, function(x){return(g$geneset.names[grepl(x,g$geneset.names)])})))
        #   dataset$x <- collapse_pathways2(dataset$x, g, keep_pathways)
        #   dataset$xx <- collapse_pathways2(dataset$xx, g, keep_pathways)
        # }else if(parameters$convert_genenames){
        #   dataset$x <- convert_genenames(dataset$x)
        #   dataset$xx <- convert_genenames(dataset$xx)
        # }
        
        features_selected[[names(list_of_datasets)[d]]] <- list("unflipped" = create_new_feature_set(EXP_RANGE_J, EXP_RANGE, colnames(dataset$x)),
                                                                "flipped_wilcox" = create_new_feature_set(EXP_RANGE_J, EXP_RANGE, colnames(dataset$x)),
                                                                "flipped_ttest" = create_new_feature_set(EXP_RANGE_J, EXP_RANGE, colnames(dataset$x)))
        
        for(j in seq(1:length(EXP_RANGE_J))){
          for(i in seq(1:length(EXP_RANGE))){
            
            flip_j = EXP_RANGE_J[j]
            flip_i = EXP_RANGE[i]
            
            logger.info(msg = paste(c("MAIN - ITER - ", "DATA = ", d, ", K = ", k, ", flip_j = ", flip_j, ", flip_i = ", flip_i), collapse="" ))
            
            Xt = dataset$x
            Xt = Xt[,colSums(Xt) != 0]  # select only the genes with > 0 reads
            yt = dataset$y
            
            Xs = dataset$xx
            Xs = Xs[,colSums(Xs) != 0]  # select only the genes with > 0 reads
            ys = dataset$tt
            
            # use only features with > 0 reads in training set and test set
            feature_names <- intersect(colnames(Xt), colnames(Xs))
            Xs = Xs[,feature_names]      # select only the genes that were present in the training set 
            Xt = Xt[,feature_names]
            
            logger.info(msg = sprintf("MAIN - DATA - TRAIN - N = %i samples; G1 = %i (%f), G2 = %i (%f)", length(yt), table(yt)[1], table(yt)[1]/length(yt), table(yt)[2], table(yt)[2]/length(yt)))
            logger.info(msg = sprintf("MAIN - DATA - TEST - N = %i samples; G1 = %i (%f), G2 = %i (%f)", length(ys), table(ys)[1], table(ys)[1]/length(ys), table(ys)[2], table(ys)[2]/length(ys)))
            
            Xt = scale(Xt)  # CHECK ON THIS
            Xs = scale(Xs)  # CHECK ON THIS
            
            a = inject_label_noiseR(yt, flip_i, flip_j)
            yz = a[[1]]
            fdz = a[[2]]
            
            pca_res <- prcomp(Xt)
            plot(pca_res$x[,1],pca_res$x[,2], col = pca_cols[yt] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz < 0) + 1],)
            
            Xt_pcatrans <- pca_res$x
            
            # filter with PC
            filter <- make_ensemble(Xt_pcatrans, y = unlist(lapply(yz, function(x){if(x==1){return("one")}else if(x==2){return("two")}})), 
                                    c("rf", "knn", "svmLinear3"),multiple = FALSE)  #, , ,"nnet","regLogistic" 
            
            # in this plot, the larger points are the ones that get kept.
            plot(pca_res$x[,1],pca_res$x[,2], col = pca_cols[yt] ,lwd = 2, xlab = "PC1",ylab = "PC2", pch = c(16, 1)[as.integer(fdz < 0) + 1],
                 cex = c(1.0, 2.0)[as.integer(filter$MF) + 1])
            
            if(length(table(fdz)) > 1){
              print("F")
              a <- flag_flipped_samples(filter$MF, fdz)
              logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - MF - %s", toString(a$table)))
              CMs_MF[dimr, dsize, d, k, j, i] = paste(as.vector(a$table), collapse="_")
              #a <- flag_flipped_samples(filter$CF, fdz)
              #logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - CF - %s", toString(a$table)))
              #CMs_CF[dimr, dsize, d, k, j, i] = paste(as.vector(a$table), collapse="_")
            }else{
              print("G")
              CMs_MF[dimr, dsize, d, k, j, i] = paste(c(1,1,1,1), collapse="_")
              #CMs_CF[dimr, dsize, d, k, j, i] = paste(c(1,1,1,1), collapse="_")
              logger.info(msg = "flip_i and flip_j = 0, so no QC metric available:")
              logger.info(msg = paste("MF: ", toString(table(filter$MF))))
              logger.info(msg = paste("CF: ", toString(table(filter$CF))))
            }
            
            idx_keep <- which(filter$MF)
            Xt_confident <- Xt[idx_keep,]
            yz_confident <- yz[idx_keep]
            
            idx_sub <- shuffle(seq(1:dim(Xt)[1]))[1:round(dim(Xt)[1]/2)] 
            Xt_sub <- Xt[idx_sub,]
            yt_sub <- yt[idx_sub,]
            yz_sub <- yz[idx_sub,]
            
            # clean model
            CV <- cv.glmnet(Xt_confident, y = (yz_confident == 2), alpha = .5)
            model_clean <- glmnet(Xt_confident, y = (yz_confident == 2), family = "binomial", lambda = CV$lambda.min, alpha = .5)
            res_train <- predict(model_clean, Xt_confident, type = "response")
            roc_train <- pROC::roc((yz_confident == 2), res_train[,1])
            res_test <- predict(model_clean, Xs, type = "response")
            roc_test <- pROC::roc((ys==2), res_test[,1])
            auc_rlr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
            
            # # subset unflipped model
            # CV <- cv.glmnet(Xt_sub_unflipped, y = (yz_sub_unflipped == 2), alpha = .5)
            # model_sub <- glmnet(Xt_sub_unflipped, y = (yz_sub_unflipped == 2), family = "binomial", lambda = CV$lambda.min, alpha = .5)
            # res_train <- predict(model_sub, Xt_sub_unflipped, type = "response")
            # roc_train <- pROC::roc(yz_sub_unflipped, res_train[,1])
            # res_test <- predict(model_sub, Xs, type = "response")
            # roc_test <- pROC::roc(ys, res_test[,1])
            # auc_gammalr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
            
            # subset unflipped model
            CV <- cv.glmnet(Xt_sub, y = (yt_sub == 2), alpha = .5)
            model_sub <- glmnet(Xt_sub, y = (yt_sub == 2), family = "binomial", lambda = CV$lambda.min, alpha = .5)
            res_train <- predict(model_sub, Xt_sub, type = "response")
            roc_train <- pROC::roc(yt_sub, res_train[,1])
            res_test <- predict(model_sub, Xs, type = "response")
            roc_test <- pROC::roc((ys==2), res_test[,1])
            auc_gammalr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
            
            # noisy model
            CV <- cv.glmnet(Xt, y = (yz == 2), alpha = .5)
            model_full <- glmnet(Xt, y = (yz == 2), family = "binomial", lambda = CV$lambda.min, alpha = .5)
            res_train <- predict(model_full, Xt, type = "response")
            roc_train <- pROC::roc(yz, res_train[,1])
            res_test <- predict(model_full, Xs, type = "response")
            roc_test <- pROC::roc((ys==2), res_test[,1])
            auc_lr[dimr, dsize, d, k, j, i] = roc_test$auc  # save result
            
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

# close(matlab)

# # THIS IS A HEATMAP PLOT
# plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "simulate_mbal", "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
# plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)
# 
# plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = NULL, "EXP_RANGE" = NULL, "EXP_RANGE_J" = .2, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
# plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)
# 
# plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data1", "EXP_RANGE" = NULL, "EXP_RANGE_J" = .2, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
# plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)
# 
# # THIS IS A LINE PLOT
# plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data1", "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file2", sep=""), "performance_metric" = "AUC")
# plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)
# plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data2", "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file3", sep=""), "performance_metric" = "AUC")
# plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)


# 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , , x, 1, 1])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, , 1])  # DS_size as a function of error rate in J
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, 1 , ])  # DS_size as a function of error rate in I, with 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 1, x, 1 , ])  # all DIM x all error rate I, for x iterations, with DS_SIZE = DS_SIZE[1] and ERROR_J = ERROR_J[1]
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 2, x, 1 , ])  # all DIM as a function of error rate in I, with DS_SIZE = DS_SIZE[2] and ERROR_J = ERROR_J[1]
# 



# length(intersect(rownames(head(rfs_flipped$res, n=1000)),rownames(head(rfs_unflipped$res, n=1000))))/1000
# length(intersect(rownames(head(rfs_unflipped$res, n=1000)),union( rownames(head(rfs_flipped_ttest$res, n=1000)),rownames(head(rfs_flipped_roc$res, n=1000)))))/1000
# length(union(rownames(head(rfs_flipped_ttest$res, n=1000)),rownames(head(rfs_flipped_roc$res, n=1000))))

# overlap <- plot_feature_similarity(features_selected, c(100, 500, 1000, 2500, 5000, 10000))#, "fileroot" = paste(EXPERIMENT_DIR, "features_byflip", sep=""))



plot_confusion_matrix(CMs_MF, colorRampPalette(c("white","dodgerblue"))(n = 401), max=400)
#plot_confusion_matrix(CMs_CF, colorRampPalette(c("white","purple"))(n = 301), max = 300)


plot_params <- list("DS_SIZE_RANGE" = 2000, "DIM_RANGE" = 200, "DATASET" = "DATASETS_viral_bacterial", "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)



