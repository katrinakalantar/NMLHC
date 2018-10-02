library(R.matlab)   # enables cross-talk with matlab server
library(PopED)      # contains feval function
library(matlab)     # contains repmat function equivalent to matlab version
library(pROC)       # contains roc function
library(glmnet)     # contains glmnet function for regularized regression in standard R library
library(gplots)     # contains heatmap.2 function

source("/Users/kkalantar/Documents/Research/NMLHC/nmlhc_matlab2R_functions.R")
source("/Users/kkalantar/Documents/Research/NMLHC/pylogger.R")  #sourced from: https://gist.github.com/jonathancallahan/3ed51265d3c6d56818458de95567d3ae

#
# Open the Matlab server
#

options(matlab="/Applications/MATLAB_R2018a.app/bin/matlab")  # set the location to be used by all matlab calls
matlab <- Matlab()
if(isOpen(matlab)){ close(matlab) }      # if a previous server is open, close it now.
Matlab$startServer()
isOpen <- open(matlab)
load_matlab_libraries()   # load the .m files required for running the analysis / estimation procedures

#
# set global options
# 

#### could we have a "params" object that just gets loaded and then we have params$iter, params$cls, params$ds_size, params$dataset_param, params$n_features, etc.

#INPUT_DIR = ;
#OUTPUT_DIR = ;
#logger.setup(infoLog = paste(OUTPUT_DIR, "info.txt"))
# logger.info(msg = "")


ITER = 2                             # number of iterations to run (to obtain mean performance)
CLS = 2                              # number of classes
DATASET_PARAM = 'generate_data'      # data generation method 

use_PCs = FALSE
feature_select = TRUE

common_reg = 'lasso'                 # regularization type
common_sn = 1e-8;                    # 
common_maxIter = 1000                # max iterations for the algorithm

# parameters iterated
#EXP_RANGE       = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # y-axis of heatmap
#EXP_RANGE_J     = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # x-axis of heatmap

DS_SIZE_RANGE   = c(200)
DIM_RANGE       = c(2000)
EXP_RANGE       = c(10, 50, 100, 500, 1000, 2500, 5000) # this range was used for sample sizes
EXP_RANGE_J     = c(.3, .2); # x-axis of heatmap

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
    for(k in seq(1:ITER)){
      for(j in seq(1:length(EXP_RANGE_J))){
        for(i in seq(1:length(EXP_RANGE))){
          
          # DO STUFF IN HERE
          
          DIM = DIM_RANGE[dimr];  # n_features
          DS_SIZE = DS_SIZE_RANGE[dsize];
          flip_j = EXP_RANGE_J[j]
          flip_i = EXP_RANGE[i]
          
          # obtain the data to be used in this analysess
          dataset <- feval(DATASET_PARAM, CLS, DIM, DS_SIZE, 1000, 1, "gen")
          
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


#
# Plotting Functions
#

my_palette <- colorRampPalette(c("purple3", "white"))(n = 100)
my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)
my_palette3 <- colorRampPalette(c("blue", "white","red"))(n=100)

# THIS IS A HEATMAP PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" = NULL)
plot_relevant_data(list(lr = auc_lr, rlr = auc_rlr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, fileroot, performance_metric = "AUC", colorpal = my_palette2)

# THIS IS A LINE PLOT
plot_params <- list("DS_SIZE_RANGE" = 20, "DIM_RANGE" = 100, "EXP_RANGE" = NULL, "EXP_RANGE_J" =.3)
plot_relevant_data(list(lr = auc_lr, rlr = auc_rlr, gammalr = auc_gammalr), plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, fileroot, performance_metric = "AUC", colorpal = my_palette2, color_by_pval = TRUE)

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





