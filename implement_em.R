
#
# Implementation of EM Algorithm, per documentation here: https://github.com/ihh/noisy-logit/blob/master/paper.pdf
# Katrina Kalantar
# September 12, 2018
#

library('e1071')                                        # contains sigmoid function
library("mlbench")                                      # contains BreastCancer dataset (and other ML datasets)
library("glmnet")   


#
# Set-up Data 
#

data(BreastCancer) 

# create dummy training / test dataset
train <- sample(seq(1,699),400)

Xt <- t(data.matrix(BreastCancer[train,seq(2,ncol(BreastCancer)-1)]))
Xt[is.na(Xt)] <- 0
Yt <- as.character(BreastCancer[train,ncol(BreastCancer)])
Yt[seq(1,round(length(Yt)/4))] = c("maybe","unlikely")[as.integer(Yt[seq(1,round(length(Yt)/4))] == "benign")+1]  # convert benign to unlikely, malignant to maybe
# NOTE: order the levels in accordance with the "true" labels
# (ie 1 == negative, 2 == positive, 3... are additional observed classes that do not correspond to a true latent label)
Ytf <- as.integer(factor(Yt, levels = c("benign","malignant","maybe","unlikely")))

Xs<- t(data.matrix(BreastCancer[!(rownames(BreastCancer) %in% train),seq(2,ncol(BreastCancer)-1)]))
Xs[is.na(Xs)] <- 0
Ys <- as.character(BreastCancer[!(rownames(BreastCancer) %in% train),ncol(BreastCancer)])
Ys[seq(1,round(length(Ys)/4))] = c("maybe","unlikely")[as.integer(Ys[seq(1,round(length(Ys)/4))] == "benign")+1]  # convert benign to unlikely, malignant to maybe
Ysf <- as.integer(factor(Ys), levels = c("benign","malignant","maybe","unlikely"))


#
# Run EM
#

# algorithm parameters
n_classes_observed = 4                                  # number of observed classes (clinician labels)
n_true_classes = 2                                      # number of true classes
max_iter = 1000                                         # maximum number of iterations of EM before quitting with "did not converge"

### EM Algorithm Implementation

# Set theta_1 to some sensible values, where theta_1 = (Z, u)
m <- matrix(rnorm(n_true_classes * n_classes_observed), nrow=n_classes_observed)    # initialize Z probability parameter matrix
zmat <- exp(m)/rowSums(exp(m))
u <- as.matrix(runif(nrow(Xs), 0, 1))                   # initialize u weight vector

zmat_orig <- zmat                                       # save the original parameters theta_1 for reference
u_orig <- u

probs <- list()
continue = TRUE

### for i in {1,2,3...}
for(i in seq(1:max_iter)){
  
  sigma <- sigmoid(t(u) %*% Xt)
  B = 1 / (1 + (((1-sigma) * zmat[Ytf,1]) / (sigma * zmat[Ytf,2])))
  
  ### Set z(i+1) ← argmaxz(Ez) by counting & normalizing
  gamma <- matrix(data = NA, nrow = nrow(zmat), ncol = ncol(zmat))
  for(obs in seq(1:n_classes_observed)){
    gamma[obs, 1] = as.integer(obs==1) + sum(1 - B[Ytf == obs])
    gamma[obs, 2] = as.integer(obs==2) + sum(B[Ytf == obs])
  }
  zmat_new <- t(t(gamma)/colSums(gamma))      # or; cbind(gamma[,1] / sum(gamma[,1]), gamma[,2] / sum(gamma[,2]))
  zmat_new <- zmat_new/rowSums(zmat_new)      # normalize rows?? NOT SURE IF THIS IS CORRECT
  
  ### Set u(i+1) ← argmaxu(Eu) using glmnet() with weights
  pseudo_dataset <- cbind(Xt, Xt)
  pseudo_labels <- cbind(rep(0, ncol(Xt)), rep(1, ncol(Xt)))
  weights_f <- cbind((1-B), B)
  model <- glmnet(t(pseudo_dataset), as.numeric(pseudo_labels), family="binomial", weights = as.numeric(weights_f)) 
  u_new <- model$beta[,ncol(model$beta)-3]  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, just selected third-to-last lambda row randomly for now
  
  P_b1 = sigmoid(t(u_new) %*% Xt)
  P_b0 = 1 - P_b1
  P_u = prod(exp(-abs(u_new)))
  P_z = (1-sum(zmat_new[,1])) * prod(zmat_new[,1] ^ c(0,0,0,1)) * (1-sum(zmat_new[,2])) * prod(zmat_new[,2] ^ c(1,0,0,0))
  new_prob = P_u * P_z * prod(as.numeric((P_b1 * zmat_new[,2][Ytf]) + (P_b0 * zmat_new[,1][Ytf])))

  P_b1_old = sigmoid(t(u) %*% Xt)
  P_b0_old = 1 - P_b1_old
  P_u_old = prod(exp(-abs(u)))
  P_z_old = (1-sum(zmat[,1])) * prod(zmat[,1] ^ c(0,0,0,1)) * (1-sum(zmat[,2])) * prod(zmat[,2] ^ c(1,0,0,0))
  previous_prob = P_u_old * P_z_old * prod(as.numeric((P_b1_old * zmat[,2][Ytf]) + (P_b0_old * zmat[,1][Ytf])))
  
  probs[[i]] <- new_prob/previous_prob
  zmat <- zmat_new
  u <- u_new
  
  # print statements for debugging 
  print(i)
  print((new_prob/previous_prob))
  print(new_prob)
  print(zmat)
  
  # the error seems to go up/down wildly in the initial rounds (1-10), so require at least 10 rounds for now
  if( i>2 & abs(new_prob)/abs(previous_prob) < .0001) break
  
}

plot(seq(2:length(probs)), probs[2:length(probs)], ylim = c(-10,10), xlab = "EM Iteration", ylab = "prob")
print(zmat)