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

source("/Users/kkalantar/Documents/Research/NMLHC/nmlhc_matlab2R_functions.R")
source("/Users/kkalantar/Documents/Research/NMLHC/pylogger.R")  #sourced from: https://gist.github.com/jonathancallahan/3ed51265d3c6d56818458de95567d3ae


#
# set the experiment directory - this will contain the input parameter file and is where all the output will be written
#

EXPERIMENT_DIR <- "/Users/kkalantar/Documents/Research/NMLHC/EXPERIMENTS/Experiment_1/"

parameters <- read_json(paste(EXPERIMENT_DIR, "parameters2.json", sep=""), simplifyVector = TRUE)
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


ITER = parameters$ITER #2                             
CLS = parameters$CLS #2                              # number of classes
DATASET_PARAM = parameters$DATASET_PARAM # 'generate_data'      # data generation method 

use_PCs = parameters$use_PCs #FALSE
feature_select = parameters$feature_select #TRUE

common_reg = parameters$common_reg #'lasso'                 # regularization type
common_sn = parameters$common_sn # 1e-8;                    # 
common_maxIter = parameters$common_maxIter # 1000                # max iterations for the algorithm

# parameters iterated
#EXP_RANGE       = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # y-axis of heatmap
#EXP_RANGE_J     = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # x-axis of heatmap

DS_SIZE_RANGE   = parameters$DS_SIZE_RANGE #c(200)
DIM_RANGE       = parameters$DIM_RANGE #c(2000)
EXP_RANGE       = parameters$EXP_RANGE #c(10, 50, 100, 500, 1000, 2500, 5000) # this range was used for sample sizes
EXP_RANGE_J     = parameters$EXP_RANGE_J #c(.3, .2); # x-axis of heatmap

#err_lr = list(); err_lr_nonoise = list(); err_rlr = list(); err_gammalr =  list ();     # preallocating error storage
#auc_lr = list(); auc_lr_nonoise = list(); auc_rlr = list(); auc_gammalr = list();       # preallocating AUC storage

err_lr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);
err_rlr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);
err_gammalr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);

auc_lr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);
auc_rlr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);
auc_gammalr = create_new_result_set(DIM_RANGE, DS_SIZE_RANGE, EXP_RANGE_J, EXP_RANGE);

#
# Run analysis with given parameters, track results
#

for(dimr in seq(1:length(DIM_RANGE))){
  for(dsize in seq(1:length(DS_SIZE_RANGE))){
    for(k in seq(1:ITER)){           # number of iterations to run (to obtain mean performance)
      for(j in seq(1:length(EXP_RANGE_J))){
        for(i in seq(1:length(EXP_RANGE))){
          
          DIM = DIM_RANGE[dimr];  # n_features
          DS_SIZE = DS_SIZE_RANGE[dsize];
          flip_j = EXP_RANGE_J[j]
          flip_i = EXP_RANGE[i]
          
          logger.info(msg = paste(c("MAIN - ITER - ", "DIM = ", DIM, ", DS_SIZE = ", DS_SIZE, 
                                    ", flip_j = ", flip_j, ", flip_i = ", flip_i, ", K = ", k), collapse="" ))
          
          # obtain the data to be used in this analysess
          dataset <- feval(DATASET_PARAM, as.numeric(CLS), as.numeric(DIM), as.numeric(DS_SIZE), 1000, 1, "gen")
          
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
          
          # create a new winit for training the noised model
          winit = randn(dim(Xt)[2] + 1, 1)
          
          # rLR (using true labels) to give baseline of number of features
          options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
          options$estG = parameters$estG #FALSE
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
          err_lr[dimr, dsize, k, j, i] = result$error
          auc_lr[dimr, dsize, k, j, i] = result$auc$auc
          
          # rLR (using flipped labels) with estimation of G 
          options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
          options$estG = TRUE;
          options$regFunc = common_reg;
          options$sn = common_sn;
          rr = .2;
          result = run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
          err_rlr[dimr, dsize, k, j, i] = result$error
          auc_rlr[dimr, dsize, k, j, i] = result$auc$auc
          
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


save_data(auc_lr, get_variable_name(auc_lr), EXPERIMENT_DIR)
save_data(auc_rlr, get_variable_name(auc_rlr), EXPERIMENT_DIR)
save_data(err_lr, get_variable_name(err_lr), EXPERIMENT_DIR)
save_data(err_rlr, get_variable_name(err_rlr), EXPERIMENT_DIR)


#
# Plotting Functions
#

my_palette <- colorRampPalette(c("purple3", "white"))(n = 1001)
#my_palette3 <- colorRampPalette(c("blue", "white","red"))(n=100)
my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)    # continuous color scale for heatmaps
discrete_pal <- brewer.pal(6, "Spectral")                          # color palette for discrete variables - line plots for different conditions

# THIS IS A HEATMAP PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file1", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette2)

# THIS IS A HEATMAP PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file4", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = err_lr, rlr = err_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette)

# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file2", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)

# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3, "fileroot" = paste(EXPERIMENT_DIR, "file5", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = err_lr, rlr = err_rlr, gammalr = err_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, colorpal = discrete_pal, color_by_pval = TRUE)


# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = .1, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file3", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J,  colorpal = discrete_pal, color_by_pval = TRUE)
# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = .1, "EXP_RANGE_J" = NULL, "fileroot" = paste(EXPERIMENT_DIR, "file6", sep=""), "performance_metric" = "AUC")
plot_relevant_data(list(lr = err_lr, rlr = err_rlr, gammalr = err_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J,  colorpal = discrete_pal, color_by_pval = TRUE)

close(matlab)




# 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , , x, 1, 1])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, , 1])  # DS_size as a function of error rate in J
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ 1, , x, 1 , ])  # DS_size as a function of error rate in I, with 
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 1, x, 1 , ])  # all DIM x all error rate I, for x iterations, with DS_SIZE = DS_SIZE[1] and ERROR_J = ERROR_J[1]
# a = lapply(seq(dim(result_array)[3]), function(x) result_array[ , 2, x, 1 , ])  # all DIM as a function of error rate in I, with DS_SIZE = DS_SIZE[2] and ERROR_J = ERROR_J[1]
# 


## TESTING OUT WILCOX TEST FUNCTION
# raw_data <- list()
# result_arrays <- list(err_lr,err_rlr, err_gammalr)
# for(i in seq(1:length(result_arrays))){
#   a = lapply(seq(dim(result_arrays[[i]])[3]), function(x) result_arrays[[i]][index_1, index_2, x, index_4, index_5])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
#   raw_data[[i]] <- do.call(rbind, a)
# }
# ew <- evaluate_wilcox(raw_data)
# # if you have a matrix w that is the wilcox p-value results...
# which(w<0.05,arr.ind = T) # gives you the index values that are significantly different







