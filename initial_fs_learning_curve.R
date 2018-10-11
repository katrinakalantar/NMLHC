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

source("/Users/kkalantar/Documents/Research/NMLHC/nmlhc_matlab2R_functions.R")
source("/Users/kkalantar/Documents/Research/NMLHC/pylogger.R")  #sourced from: https://gist.github.com/jonathancallahan/3ed51265d3c6d56818458de95567d3ae


#
# set the experiment directory - this will contain the input parameter file and is where all the output will be written
#

EXPERIMENT_DIR <- "/Users/kkalantar/Documents/Research/NMLHC/EXPERIMENTS/Experiment_1/"

parameters <- read_json(paste(EXPERIMENT_DIR, "parameters.json", sep=""), simplifyVector = TRUE)
logger.setup(infoLog = paste(EXPERIMENT_DIR, "info.txt"))
init_log(parameters)


#
# Open the Matlab server
#

options(matlab="/Applications/MATLAB_R2018a.app/bin/matlab")  # set the location to be used by all matlab calls
matlab <- Matlab()
if(isOpen(matlab)){ close(matlab) }      # if a previous server is open, close it now.
Matlab$startServer()
isOpen <- open(matlab)
load_matlab_libraries()   # load the .m files required for running the analysis / estimation procedures

logger.info(paste("MAIN - matlab server status ", toString(isOpen)))


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

DS_SIZE_RANGE   = parameters$DS_SIZE_RANGE 
DIM_RANGE       = parameters$DIM_RANGE
EXP_RANGE       = parameters$EXP_RANGE
EXP_RANGE_J     = parameters$EXP_RANGE_J

#result$auc$auc



#
# Run analysis with given parameters, track results
#

      
## THIS ASSUMES THAT LEARNING CURVE WILL ONLY BE GENERATED ON more currated datasets
list_of_datasets <- list()
for(x in rownames(parameters$datasets)){
  d <- parameters$datasets[x,]
  if(d$type == "gen"){
    DIM = 20000
    DS_SIZE = 100
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

for(d in seq(1:length(list_of_datasets))){
  
  dataset <- list_of_datasets[[names(list_of_datasets)[d]]]
  
  learning_curve_iters <- seq(1,dim(dataset$x)[1],5)
  learning_curve_iters <- learning_curve_iters[2:length(learning_curve_iters)]
  
  # preallocating error storage
  err_lr = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  err_rlr = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  err_gammalr = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  
  # preallocating AUC storage
  auc_full = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  auc_sub = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  auc_sub_filt = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  auc_full_filt = create_new_result_set_learning_curve(parameters$datasets$name, learning_curve_iters);
  
  Xt = dataset$x
  yt = dataset$y
  Xs = dataset$xx
  ys = dataset$tt
  
  for(k in seq(1:ITER)){           # number of iterations to run for the learning curve
    
    # LEARNING CURVE FOR FEATURE SELECTION
    for(s in seq(1:length(learning_curve_iters))){
      
      logger.info(msg = paste(c("FUNCTION - GENERATING LC - DATASET: ", d, ", ITER: ", k, 
                                ", LC: ", learning_curve_iters[s]), collapse=""))
      
      subs <- learning_curve_iters[s]
      
      # experiment 2
      flip_i <- .2
      flip_j <- .2
      a = inject_label_noise(yt, flip_i, flip_j)
      yz = a[[1]]
      fdz = a[[2]] # default is -1, if it's been flipped it changes to +1
      
      Xtsub <- NULL
      ytsub <- c(1)
      yzsub <- c(1)
      fdzsub <- c(1)
      while(length(unique(ytsub)) < 2){             # resample until you have at least one of each class
        sub_index <- sample(seq(dim(Xt)[1]), subs)
        Xtsub <- Xt[sub_index,]
        ytsub <- yt[sub_index]
        
        # experiment 2
        yzsub <- yz[sub_index]    
        fdzsub <- fdz[sub_index]
      }
      
      # experiment 2
      # identify correctly labeled samples
      # assume you can identify correctly labelled training samples with 100% accuracy
      Xtsub_correctlabels <- Xtsub[fdzsub < 0,]
      yzsub_correctlabels <- yzsub[fdzsub < 0]
      logger.info(msg = paste(c("# of correct samples: ", sum(fdzsub < 0), ", # of total samples: ", learning_curve_iters[s]), collapse = ""))
      
      # experiment 2
      # standardize the subset of the subset data where labels are correct
      a = standardise(Xtsub_correctlabels, Xs)
      Xtsub_correctlabels = a[[1]]
      Xssub_correctlabels = a[[2]]
      
      # standardize the subset data
      a = standardise(Xtsub, Xs)
      Xtsub = a[[1]]
      Xssub = a[[2]]
      
      # standardize the full datasets
      a = standardise(Xt, Xs)
      Xt = a[[1]]
      Xs = a[[2]]
      
      colnames(Xt) <- colnames(dataset$x)   # must have labels for the features_selection() function to work properly
      colnames(Xs) <- colnames(dataset$xx)
      colnames(Xssub) <- colnames(dataset$xx)
      colnames(Xtsub) <- colnames(dataset$x)
      rownames(Xt) <- rownames(dataset$x)
      rownames(Xtsub) <- rownames(dataset$x)[sub_index]
      rownames(Xs) <- rownames(dataset$xx)
      
      # experiment 2
      colnames(Xtsub_correctlabels) <- colnames(dataset$x)
      colnames(Xssub_correctlabels) <- colnames(dataset$xx)
      rownames(Xssub) <- rownames(dataset$xx)
      rownames(Xtsub_correctlabels) <- rownames(Xtsub)[fdzsub < 0]
      rownames(Xssub_correctlabels) <- rownames(dataset$xx)
      
      # experiment 1
      # compute features based on the LC subset (as opposed to the full dataset)
      #rfs_unflipped <- run_FS(Xtsub, ytsub, "wilcoxon")
      #keep_vars <- head(rownames(rfs_unflipped), n = 5000)
      
      # experiment 2
      # compute features based on the LC subset where labels are correct (as opposed to the total subset)
      rfs_unflipped <- run_FS(Xtsub_correctlabels, yzsub_correctlabels, "wilcoxon")
      keep_vars <- head(rownames(rfs_unflipped), n = 1000)
      
      # experiment 1
      #Xtsub_filt <- Xtsub[,keep_vars]
      #Xtfull_filt <- Xt[,keep_vars]
      #Xssub_filt <- Xssub[,keep_vars]
      #Xs_filt <- Xs[,keep_vars]
      
      # experiment 2
      Xtsub_correctlabels_filt <- Xtsub_correctlabels[,keep_vars]
      Xtsub_filt <- Xtsub[,keep_vars]
      Xssub_correctlabels_filt <- Xssub_correctlabels[,keep_vars]
      Xssub_filt <- Xssub[,keep_vars]
      
      # create a new winit for training the noised model
      winit = randn(dim(Xt)[2] + 1, 1)
      winit_filt = randn(dim(Xtsub_filt)[2] + 1, 1)
      
      #
      # EXPERIMENT 2
      #
      
      # LR training on the full dataset without feature selection up front
      options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      options$estG = FALSE
      options$regFunc = common_reg
      options$sn = common_sn
      ginit = cbind(c(1,0),c(0,1))
      result = run_rlr("rlr", winit, ginit, Xtsub_correctlabels, yzsub_correctlabels, options, Xssub_correctlabels, ys)
      auc_sub[d, k, s] <- result$auc$auc
      
      # LR training on just the subset of data after filtering features based on subset of CORRECT SAMPLES
      options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      options$estG = FALSE
      options$regFunc = common_reg
      options$sn = common_sn
      ginit = cbind(c(1,0),c(0,1))
      result = run_rlr("rlr", winit_filt, ginit, Xtsub_correctlabels_filt, yzsub_correctlabels, options, Xssub_correctlabels_filt, ys)
      auc_sub_filt[d, k, s] <- result$auc$auc
      
      # LR training on the full dataset after filtering features based on subset
      options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      options$estG = TRUE
      options$regFunc = common_reg
      options$sn = common_sn
      ginit = cbind(c(.8,.2),c(.2,.8))
      result = run_rlr("rlr", winit_filt, ginit, Xtsub_filt, yzsub, options, Xssub_filt, ys)
      auc_full_filt[d, k, s] <- result$auc$auc
      
      
      
      # #
      # # EXPERIMENT 1
      # #
      # 
      # # # LR training on the full dataset without feature selection up front
      # # options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      # # options$estG = FALSE
      # # options$regFunc = common_reg
      # # options$sn = common_sn
      # # ginit = cbind(c(1,0),c(0,1))
      # # result = run_rlr("rlr", winit, ginit, Xt, yt, options, Xs, ys)
      # # auc_full[d, k, s] <- result$auc$auc
      # 
      # # LR training on just the subset of data without feature selection up front
      # options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      # options$estG = FALSE
      # options$regFunc = common_reg
      # options$sn = common_sn
      # ginit = cbind(c(1,0),c(0,1))
      # result = run_rlr("rlr", winit, ginit, Xtsub, ytsub, options, Xssub, ys)
      # auc_sub[d, k, s] <- result$auc$auc
      # 
      # # LR training on just the subset of data after filtering features based on subset
      # options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      # options$estG = FALSE
      # options$regFunc = common_reg
      # options$sn = common_sn
      # ginit = cbind(c(1,0),c(0,1))
      # result = run_rlr("rlr", winit_filt, ginit, Xtsub_filt, ytsub, options, Xssub_filt, ys)
      # auc_sub_filt[d, k, s] <- result$auc$auc
      # 
      # # LR training on the full dataset after filtering features based on subset
      # options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      # options$estG = FALSE
      # options$regFunc = common_reg
      # options$sn = common_sn
      # ginit = cbind(c(1,0),c(0,1))
      # result = run_rlr("rlr", winit_filt, ginit, Xtfull_filt, yt, options, Xs_filt, ys)
      # auc_full_filt[d, k, s] <- result$auc$auc
      
      
    }
  }
}

close(matlab)

#plot_multiple_learning_curves(list("AUC FULL FILT" = auc_full_filt, "AUC SUB FILT" = auc_sub_filt), c("magenta","orange"), learning_curve_iters)
#plot_multiple_learning_curves(list("AUC FULL FILT" = auc_full_filt, "AUC SUB FILT" = auc_sub_filt, "AUC SUB" = auc_sub), c("magenta","orange","gold"), learning_curve_iters)
#plot_multiple_learning_curves(list("AUC FULL FILT" = auc_full_filt, "AUC SUB FILT" = auc_sub_filt, "AUC SUB" = auc_sub, "AUC FULL" = auc_full), c("magenta","orange","gold", "blue"), learning_curve_iters, plot_all_points = FALSE)
#plot_multiple_learning_curves(list("AUC FULL FILT" = auc_full_filt, "AUC SUB" = auc_sub), c("magenta","gold"), learning_curve_iters)


save_data(auc_full, get_variable_name(auc_full), EXPERIMENT_DIR)
save_data(auc_full_filt, get_variable_name(auc_full_filt), EXPERIMENT_DIR)
save_data(auc_sub_filt, get_variable_name(auc_sub_filt), EXPERIMENT_DIR)
save_data(auc_sub, get_variable_name(auc_sub), EXPERIMENT_DIR)

# Experiment 2
pdf(paste(EXPERIMENT_DIR, "experiment2_v1.pdf"), height = 5, width = 8)
plot_multiple_learning_curves(list("AUC FULL FILT (rLR)" = auc_full_filt, "AUC SUB FILT (LR)" = auc_sub_filt), c("green","turquoise"), learning_curve_iters, success = 1:7, title = "Learning Curve:\n BACTERIA_VIRUS, n_features_filtered = 1000")
plot_multiple_learning_curves(list("AUC SUB (LR)" = auc_sub, "AUC SUB FILT (LR)" = auc_sub_filt), c("blue","turquoise"), learning_curve_iters, success = 1:7, title = "Learning Curve:\n BACTERIA_VIRUS, n_features_filtered = 1000")
plot_multiple_learning_curves(list("AUC SUB (LR)" = auc_sub, "AUC FULL FILT (rLR)" = auc_full_filt), c("blue","green"), learning_curve_iters, success = 1:7, title = "Learning Curve:\n BACTERIA_VIRUS, n_features_filtered = 1000")
plot_multiple_learning_curves(list("AUC FULL FILT (rLR)" = auc_full_filt, "AUC SUB FILT (LR)" = auc_sub_filt, "AUC SUB (LR)" = auc_sub), c("green","turquoise","blue"), learning_curve_iters, success = 1:7, title = "Learning Curve:\n BACTERIA_VIRUS, n_features_filtered = 1000")
dev.off()
