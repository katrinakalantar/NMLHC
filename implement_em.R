
#
# Implementation of EM Algorithm, per documentation here: https://github.com/ihh/noisy-logit/blob/master/paper.pdf
# Katrina Kalantar
# September 24, 2018
#

library('e1071')                                        # contains sigmoid function
library("mlbench")                                      # contains BreastCancer dataset (and other ML datasets)
library("glmnet")                                       # contains regularized regression function
library("datasets")                                     # contains the iris dataset
library("DirichletReg")

#
# Helper functions
#

modify_input_labels <- function(Ytf, z){
  col_one <- cumsum(round(length(Ytf[Ytf==1])*z[,1]))
  col_two <- cumsum(round(length(Ytf[Ytf==2])*z[,2]))
  idx_one <- which(Ytf==1)
  idx_two <- which(Ytf==2)
  Ytf[c(idx_one[1:col_one[1]-1], idx_two[1:col_two[1]-1])] <- 1
  Ytf[c(idx_one[col_one[1]:col_one[2]-1], idx_two[col_two[1]:col_two[2]-1])] <- 2
  Ytf[c(idx_two[col_two[2]:col_two[3]-1], idx_one[col_one[2]:col_one[3]-1])] <- 3
  Ytf[c(idx_two[col_two[3]:col_two[4]-1], idx_one[col_one[3]:col_one[4]-1])] <- 4
  return(Ytf)
}
recompute_z <- function(Ytf, Ytf2){
  return(cbind(table(Ytf2[Ytf==1])/length(Ytf2[Ytf==1]), table(Ytf2[Ytf==2])/length(Ytf2[Ytf==2])))
}



#
# Set-up Data - Iris dataset - pre-specifying z-matrix and flipping labels in the probabilities given
#

data(iris)

z <- cbind(c(.8,.1,.09,.01), c(.2,.7,.03,.07))

iris_mod <- iris[sample(nrow(iris)),]
iris_mod <- iris_mod[iris_mod$Species %in% c("setosa","virginica"),] # filter out versicolor (more similar to virginica)

train <- sample(seq(1,dim(iris_mod)[1]),.7*dim(iris_mod)[1])  # create dummy training / test dataset
Xt <- t(data.matrix(iris_mod[train,1:4]))
Yt <- as.character(iris_mod[train,5])
Ytf <- as.integer(as.factor(Yt))
Ytf2 <- modify_input_labels(Ytf, z)
z_fin <- recompute_z(Ytf,Ytf2)

print("Z estimated by EM algorithm should eventually be similar to:")
print(z_fin)

Xs <- t(data.matrix(iris_mod[!(rownames(iris_mod) %in% train),1:4]))
Ys <- as.character(iris_mod[!(rownames(iris_mod) %in% train),5])
Ysf <- as.integer(as.factor(Ys))
Ysf2 <- modify_input_labels(Ysf, z)
z_test <- recompute_z(Ysf, Ysf2)


#
# Run EM
#

# algorithm parameters
n_classes_observed = 4                                  # number of observed classes (clinician labels)
n_true_classes = 2                                      # number of true classes
max_iter = 10                                           # maximum number of iterations of EM before quitting with "did not converge"

### EM Algorithm Implementation

# create list variables to store values on every iteration so that I can track errors
probs <- list(); prob_list <- list();
sigma <- list(); B <- list(); zmat <- list(); u <- list()
P_u <- list(); P_z <- list(); E_u <- list(); E_z <- list(); E_theta <- list()
continue = TRUE

# Set theta_1 to some sensible values, where theta_1 = (Z, u)
m <- matrix(runif(n_true_classes * n_classes_observed), nrow=n_classes_observed)    # initialize Z probability parameter matrix
zmat[[1]] <- t(t(m)/colSums(m))                                                     # normalize the columns to sum to 1
u[[1]] <- as.matrix(runif(nrow(Xs), 0, 1))                                          # initialize u weight vector


#
# Functions definitiones that are used multiple times
#

compute_Pu <- function(u){
  return(prod(exp(-abs(u))))
}
# THERE may be SOMETHING WRONG WITH THIS FUNCTION - dirichlet v. just post-normalization
compute_Pz <- function(z){
  return(((1-sum(z[,1])) * prod(z[,1] ^ c(1,0,0,0))) * ((1-sum(z[,2])) * prod(z[,2] ^ c(0,1,0,0))))
}
compute_Eu <- function(sigma, B, u){
  P_u <- compute_Pu(u)
  E_u <- log(P_u) + rowSums((1-B)*log(1-sigma) + B*log(sigma))
  return(E_u)
}

# NOTE: is it OK to have added one (1) to all of these log values??
compute_Ez <- function(sigma, B, z, Ytf){  
  P_z <- compute_Pz(z)
  #E_z <- log(P_z + 1) + rowSums((1-B)*log(z[,1][Ytf] + 1) + B*log(z[,2][Ytf] + 1)) 
  E_z <- log(P_z) + rowSums((1-B)*log(z[,1][Ytf]) + B*log(z[,2][Ytf])) 
  return(E_z)
}


#
# EM Algorithm Main Loop
#

for(i in seq(1:max_iter)){
  print(paste("------- iteration:", i, " -------"))
  
  #
  # E-step
  #
  
  # calculate B[[i]] using current theta (u and z)
  sigma[[i]] <- sigmoid(t(u[[i]]) %*% Xt)     
  B[[i]] = 1 / (1 + (((1-sigma[[i]]) * zmat[[i]][Ytf,1]) / (sigma[[i]] * zmat[[i]][Ytf,2])))   # lines in paper: P(bn=1|theta, xn, cn)
  
  # compute other parameters that we want to track
  P_u[[i]] = compute_Pu(u[[i]])
  P_z[[i]] = compute_Pz(zmat[[i]])
  E_u[[i]] = compute_Eu(sigma[[i]], B[[i]], u[[i]])
  E_z[[i]] = compute_Ez(sigma[[i]], B[[i]], zmat[[i]], Ytf) 
  E_theta[[i]] = E_u[[i]] + E_z[[i]]
  prob_list[[i]] = P_u[[i]] * P_z[[i]] * prod(as.numeric((sigma[[i]] * zmat[[i]][,2][Ytf]) + (1-sigma[[i]] * zmat[[i]][,1][Ytf])))
  
  # verify increase in likelihoods per iteration; also, check if time to break loop
  if(i > 1){
    print(paste("E_u increases each interation:", (E_u[[i]] > E_u[[i-1]]))) # E_u should increase after "u <- u_new"
    print(paste("E_z increases each interation:", (E_z[[i]] > E_z[[i-1]]))) # E_z should increase after "zmat <- zmat_new"
    #if(prob_list[[i]] / prob_list[[i-1]] < 0.000000001) break               # check for exit criteria
  }
  
  #
  # M-step
  #
  
  ### Set z(i+1) ← argmaxz(Ez) by counting & normalizing
  print(zmat[[i]])
  gamma <- matrix(data = NA, nrow = nrow(zmat[[i]]), ncol = ncol(zmat[[i]]))
  for(obs in seq(1:n_classes_observed)){
    gamma[obs, 1] = as.integer(obs==1) + sum(1 - B[[i]][Ytf == obs])
    gamma[obs, 2] = as.integer(obs==2) + sum(B[[i]][Ytf == obs])
  }
  zmat[[i + 1]] <- t(t(gamma)/colSums(gamma))      # normalize the columns to sum to 1

  ### Set u(i+1) ← argmaxu(Eu) using glmnet() with weights
  pseudo_dataset <- cbind(Xt, Xt)
  pseudo_labels <- cbind(rep(0, ncol(Xt)), rep(1, ncol(Xt)))
  weights_f <- cbind((1-B[[i]]), B[[i]])
  model <- glmnet(t(pseudo_dataset), as.numeric(pseudo_labels), family="binomial", weights = as.numeric(weights_f), alpha = 1) 
  u[[i + 1]] <- model$beta[,ncol(model$beta)-3]  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, just selected third-to-last lambda row randomly for now
  
}

# plots should show increasing E_u and E_z per iteration - THIS IS NOT THE CASE!! SOMETHING IS GOING WRONG
par(mfrow=c(1,2))
plot(unlist(E_u), main = "E_u")
plot(unlist(E_z), main = "E_z")
