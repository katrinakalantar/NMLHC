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
    list_of_geo_datasets[[d$name]] <- GEOquery::getGEO(filename=d$series_filename)
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


features_selected <- list() 



run_FS <- function(x, y, method){
  # method options are "ttest", "entropy", "bhattacharyya", "roc", "wilcoxon"
  setVariable(matlab, x = x)
  setVariable(matlab, y = y)
  setVariable(matlab, method=method)
  evaluate(matlab, "[idx, z] = rankfeatures(x', y, 'Criterion', method);")
  
  feature_values <- getVariable(matlab, "z")[[1]]
  rownames(feature_values) <- colnames(x)
  
  
  print(head(feature_values))
  ranks <- getVariable(matlab, "idx")[[1]]  
  
  result <- as.matrix(feature_values[colnames(x)[ranks],])
  
  return(result)
}







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
      #DIM = DIM_RANGE[1]
      #DS_SIZE = DS_SIZE_RANGE[1]
      list_of_datasets <- list()
      for(x in rownames(parameters$datasets)){
        d <- parameters$datasets[x,]
        if(d$type == "gen"){
          list_of_datasets[[d$name]] <- feval("generate_data", as.numeric(CLS), as.numeric(DIM), as.numeric(DS_SIZE), 1000, as.numeric(d$class_sep), "gen")
          logger.info(msg=paste(c("DATA - GEN - ", "CLS: ", CLS, ", DIM: ", DIM, ", DS_SIZE: ", DS_SIZE, ", CLASS_SEP: ", d$class_sep), collapse=""))
        }else if(d$type == "geo"){
          list_of_datasets[[d$name]] <- feval('subset_geo', d$name, list_of_geo_datasets)
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
        
        features_selected[[names(list_of_datasets)[d]]] <- list("unflipped" = create_new_feature_set(EXP_RANGE_J, EXP_RANGE, colnames(dataset$x)),
                                           "flipped" = create_new_feature_set(EXP_RANGE_J, EXP_RANGE, colnames(dataset$x)))
        
        
        for(j in seq(1:length(EXP_RANGE_J))){
          for(i in seq(1:length(EXP_RANGE))){
            
            
            # DO STUFF
            #DIM = DIM_RANGE[dimr];  
            #DS_SIZE = DS_SIZE_RANGE[dsize];
            flip_j = EXP_RANGE_J[j]
            flip_i = EXP_RANGE[i]
            
            logger.info(msg = paste(c("MAIN - ITER - ", "DIM = ", DIM, ", DS_SIZE = ", DS_SIZE, ", DATASET = ", d,
                                      ", flip_j = ", flip_j, ", flip_i = ", flip_i, ", K = ", k), collapse="" ))
            
            # obtain the data to be used in this analyses
            # dataset <- feval(DATASET_PARAM, as.numeric(CLS), as.numeric(DIM), as.numeric(DS_SIZE), 1000, 1, "gen")
            
            Xt = dataset$x
            yt = dataset$y
            Xs = dataset$xx
            ys = dataset$tt
            
            a = standardiseR(Xt, Xs)
            Xt = a[[1]]
            Xs = a[[2]]
            
            a = inject_label_noise(yt, flip_i, flip_j)
            yz = a[[1]]
            fdz = a[[2]]
            
            
            # DO feature selection on the original dataset, save features for evaluation
            rfs_unflipped <- run_FS(Xt, yt, "wilcoxon")
            features_selected[[names(list_of_datasets)[d]]][["unflipped"]][k,j,i, ] <- rownames(rfs_unflipped)
            # WHY WONT THIS WORK??? vvv
            #names(features_selected[[names(list_of_datasets)[d]]][["unflipped"]][k,j,i, ]) <- rownames(rfs_unflipped)
            
            # Do feature selection on the modified/flipped dataset, save features for evaluation
            rfs_flipped <- run_FS(Xt, yz, "wilcoxon")
            features_selected[[names(list_of_datasets)[d]]][["flipped"]][k,j,i, ] <- rownames(rfs_flipped)
            # WHY WONT THIS WORK??? vvv
            #names(features_selected[[names(list_of_datasets)[d]]][["unflipped"]][k,j,i, ]) <- rownames(rfs_flipped)
            
            
            
            
            
            
            # create a new winit for training the noised model
            winit = randn(dim(Xt)[2] + 1, 1)
            
            # rLR (using true labels) to give baseline of number of features
            options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
            options$estG = FALSE
            options$regFunc = common_reg
            options$sn = common_sn
            ginit = cbind(c(1,0),c(0,1))
            result = run_rlr("rlr", winit, ginit, Xt, yt, options, Xs, ys)
            
            # rLR (using flipped labels) with no estimation of G
            options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
            options$estG = FALSE
            options$regFunc = common_reg
            options$sn = common_sn
            result = run_rlr("rlr", winit, cbind(c(1,0),c(0,1)), Xt, yz, options, Xs, ys)
            err_lr[dimr, dsize, d, k, j, i] = result$error
            auc_lr[dimr, dsize, d, k, j, i] = result$auc$auc
            
            # rLR (using flipped labels) with estimation of G 
            options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
            options$estG = TRUE;
            options$regFunc = common_reg;
            options$sn = common_sn;
            rr = .2;
            # result <- run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys) # this was the original fn
            result <- NULL
            
            try(
              #logger.info(msg = paste("failed function repeating function: ", attempt, sep = ""))
              result <- run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
            )
            
            #attempt <- 0
            #while(is.null(result) && attempt <= 3){
            #  attempt <- attempt + 1
            #  try(
            #    #logger.info(msg = paste("failed function repeating function: ", attempt, sep = ""))
            #    result <- run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
            #  )
            #  logger.info(msg = paste("failed function repeating function: ", attempt, sep = ""))
            #}
            # result = tryCatch({
            #   run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
            # }, error = function(cond){
            #   logger.info(msg = "ERROR")
            #   logger.info(msg = cond)
            #   logger.info(msg = list(error = NULL, auc = list(auc=NULL)))
            #   print(list(error = NULL, auc = list(auc=NULL)))
            #   return(list(error = NULL, auc = list(auc=NULL)))
            # })
            
            #result <- possibly_run_rlr()
            
            logger.info(msg = paste("CHECK RESULT of tryCatch", result$error))
            
            #result = run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
            
            if(!is.null(result)){
              err_rlr[dimr, dsize, d, k, j, i] = result$error
              auc_rlr[dimr, dsize, d, k, j, i] = result$auc$auc
            }
            
            
            # # gammaLR (using flipped labels) with estimation of G 
            # options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
            # options$estG = TRUE;
            # rr = .2;
            # result = run_rlr("gammalr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
            # err_gammalr[[k]][i,j] = result$error
            # auc_gammalr[[k]][i,j] = result$auc
            
            
            
            # run the EM algorithm (R implementation)
            
            
            
            
            # # run standard LASSO regression
            # # with true labels
            # model <- glmnet(Xt, as.numeric(yt), family="binomial", alpha = 1) 
            # model$beta[,ncol(model$beta)-3]  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, just selected third-to-last lambda row randomly for now
            # # with flipped labels
            # model <- glmnet(Xt, as.numeric(yz), family="binomial", alpha = 1) 
            # model$beta[,ncol(model$beta)-3]  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, just selected third-to-last lambda row randomly for now
            
            
            
            
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
#my_palette3 <- colorRampPalette(c("blue", "white","red"))(n=100)
my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)    # continuous color scale for heatmaps
discrete_pal <- brewer.pal(6, "Spectral")                          # color palette for discrete variables - line plots for different conditions

close(matlab)

# THIS IS A HEATMAP PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "simulate_mbal", "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)

plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = NULL, "EXP_RANGE" = NULL, "EXP_RANGE_J" = .2, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)

plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data1", "EXP_RANGE" = NULL, "EXP_RANGE_J" = .2, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)

# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data1", "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file2", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "DATASET" = "generate_data2", "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file3", sep=""), "performance_metric" = "AUC")
plot_relevant_data_withdataset(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, parameters$datasets$name, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)


# 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , , x, 1, 1])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, , 1])  # DS_size as a function of error rate in J
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, 1 , ])  # DS_size as a function of error rate in I, with 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 1, x, 1 , ])  # all DIM x all error rate I, for x iterations, with DS_SIZE = DS_SIZE[1] and ERROR_J = ERROR_J[1]
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 2, x, 1 , ])  # all DIM as a function of error rate in I, with DS_SIZE = DS_SIZE[2] and ERROR_J = ERROR_J[1]
# 



plot_feature_similarity <- function(features_selected, i, j){
  for(d in features_selected){
    for(N in c(100, 500, 1000, 2500, 5000, 10000)){
      jaccard <- c(); spearman <- c(); 
      for(k in seq(1:dim(d$flipped)[1])){
        jaccard <- c(jaccard, length( intersect( names(head(sort(d$flipped[k,j,i,]),n=N)), names(head(sort(d$unflipped[k,j,i,]),n=N))) ) / N)
        spearman <- c(spearman, cor(d$flipped[k,j,i,], d$unflipped[k,j,i,], method="spearman"))
      }
      print(mean(jaccard))
      print(std(jaccard))
    }
    print(mean(spearman))
  }
}




plot_feature_similarity(features_selected, 4, 2)

