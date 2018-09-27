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

ITER = 5                             # number of iterations to run (to obtain mean performance)
CLS = 2                              # number of classes
# DIM = 1000                         # dimensionality of the dataset generated
DS_SIZE = 200                        # dataset size 
DATASET_PARAM = 'generate_data'      # data generation method 
n_features = 2000                    # number of features

use_PCs = FALSE
feature_select = TRUE

common_reg = 'lasso'                 # regularization type
common_sn = 1e-8;                    # 
common_maxIter = 1000                # max iterations for the algorithm

# parameters iterated
#EXP_RANGE       = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # y-axis of heatmap
#EXP_RANGE_J     = c(0, .01, .025, .05, .1, .15, .2, .3, .4, .5);  # x-axis of heatmap
#EXP_RANGE       = c(10 20 30 40 50 100 200 300 400 500 1000); #this range was used for sample sizes
EXP_RANGE       = c(10, 50, 100, 500, 1000, 2500, 5000) # this range was used for sample sizes
EXP_RANGE_J     = c(.3, .2); # x-axis of heatmap

err_lr = list(); err_lr_nonoise = list(); err_rlr = list(); err_gammalr =  list ();                # preallocating error storage
auc_lr = list(); auc_lr_nonoise = list(); auc_rlr = list(); auc_gammalr = list();  # preallocating AUC storage


#
# Run analysis with given parameters, track results
#

for(k in seq(1:ITER)){
  
  # add a new matrix for tracking results in this iteration
  err_lr[[k]]             = new_results_matrix(EXP_RANGE, EXP_RANGE_J)
  err_rlr[[k]]            = new_results_matrix(EXP_RANGE, EXP_RANGE_J)
  auc_lr[[k]]             = new_results_matrix(EXP_RANGE, EXP_RANGE_J)
  auc_rlr[[k]]            = new_results_matrix(EXP_RANGE, EXP_RANGE_J)
  
  for(i in seq(1:length(EXP_RANGE))){
    
    flip_i = 0;
    DIM = EXP_RANGE[i]
    n_features = DIM
    DS_SIZE = 100
    
    for(j in seq(1:length(EXP_RANGE_J))){
      
      flip_j = EXP_RANGE_J[j]
      
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
      err_lr[[k]][i,j] = result$error
      auc_lr[[k]][i,j] = result$auc$auc
      
      # rLR (using flipped labels) with estimation of G 
      options <- list(regFunc = common_reg, sn = common_sn, maxIter = common_maxIter)
      options$estG = TRUE;
      options$regFunc = common_reg;
      options$sn = common_sn;
      rr = .2;
      result = run_rlr("rlr", winit, cbind(c(1-rr, rr), c(rr, 1-rr)), Xt, yz, options, Xs, ys)
      err_rlr[[k]][i,j] = result$error
      auc_rlr[[k]][i,j] = result$auc$auc
      
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


close(matlab)



#
# Plotting Functions
#

my_palette <- colorRampPalette(c("purple3", "white"))(n = 100)
my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)
my_palette3 <- colorRampPalette(c("blue", "white","red"))(n=100)

mean_auc_lr = Reduce("+", auc_lr) / length(auc_lr)
heatmap.2(mean_auc_lr, Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette,sepcolor="black",
          colsep=0:ncol(mean_auc_lr)+1,rowsep=0:nrow(mean_auc_lr)+1, sepwidth=c(0.01,0.01), main = "AUC (Test)",
          cellnote=round(mean_auc_lr,2), notecol='black')

mean_err_lr = Reduce("+", err_lr) / length(err_lr)
heatmap.2(mean_err_lr, Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette2,sepcolor="black",
          colsep=0:ncol(mean_err_lr)+1,rowsep=0:nrow(mean_err_lr)+1, sepwidth=c(0.01,0.01), main = "Error (Test)",
          cellnote=round(mean_err_lr,2), notecol='black')



