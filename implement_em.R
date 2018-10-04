
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

Ytf <- Ytf2

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
max_iter = 50                                           # maximum number of iterations of EM before quitting with "did not converge"

### EM Algorithm Implementation

# create list variables to store values on every iteration so that I can track errors
probs <- list(); prob_list <- list();
sigma <- list(); B <- list(); zmat <- list(); u <- list(); gammalist <- list();
P_u <- list(); P_z <- list(); E_u <- list(); E_z <- list(); E_theta <- list()
continue = TRUE


# THESE TRACK whether new z and u values increase the expected value at each update.
z_better <- list()  # variable to track whether the update to z caused an increase in E_z
z_better2 <- list()
u_better <- list()  # variable to track whether the update to u caused an increase in E_u 


# Set theta_1 to some sensible values, where theta_1 = (Z, u)
m <- matrix(runif(n_true_classes * n_classes_observed), nrow=n_classes_observed)    # initialize Z probability parameter matrix
zmat[[1]] <- t(t(m)/colSums(m))                                                     # normalize the columns to sum to 1
u[[1]] <- as.matrix(runif(nrow(Xs), 0, 1)) *c(-1,1,-1,1)                                         # initialize u weight vector


#
# Functions definitions that are used multiple times
#

compute_Pu <- function(u){
  return(prod(exp(-abs(u))))
}
# THERE may be SOMETHING WRONG WITH THIS FUNCTION - dirichlet v. just post-normalization
compute_Pz <- function(z){
  #return(z[1,1] * z[2,2])  # equivalent to: return((prod(z[,1] ^ c(1,0,0,0))) * ( prod(z[,2] ^ c(0,1,0,0))))
  return((prod(z[,1] ^ c(1,.001,.001,.001))) * ( prod(z[,2] ^ c(.001,1,.001,.001))))
}
compute_Eu <- function(sigma, B, u){
  P_u <- compute_Pu(u)
  E_u <- log(P_u) + rowSums( ((1-B)*log(1-sigma)) + (B*log(sigma)) )
  return(E_u)
}
compute_Ez <- function(sigma, B, z, Ytf){  
  P_z <- compute_Pz(z)
  E_z <- log(P_z) + rowSums( ((1-B)*log(z[Ytf,1])) + (B*log(z[Ytf,2])) )
  return(E_z)
}
compute_B <- function(sigma, z, Ytf){
  return( 1 / (   1 + (  ( (1-sigma) * z[Ytf,1] ) / ( sigma * z[Ytf,2] )  )   ) )
}

# there are two formulae for Ez, attempting to evaluate the second one listed in the paper.pdf
compute_Ez2 <- function(sigma, B, z, Ytf){
  gamma <- matrix(data = NA, nrow = nrow(z), ncol = ncol(z))
  for(obs in seq(1:n_classes_observed)){
    gamma[obs, 1] = max(as.integer(obs==1), .001) + sum(1 - B[Ytf == obs])
    gamma[obs, 2] = max(as.integer(obs==2), .001) + sum(B[Ytf == obs])
  }
  Ez = sum(colSums(gamma * z))
  return(Ez)
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
  B[[i]] = compute_B(sigma[[i]], zmat[[i]], Ytf)    # lines in paper: P(bn=1|theta, xn, cn)
  
  # compute other parameters that we want to track
  P_u[[i]] = compute_Pu(u[[i]])
  P_z[[i]] = compute_Pz(zmat[[i]])
  E_u[[i]] = compute_Eu(sigma[[i]], B[[i]], u[[i]])
  E_z[[i]] = compute_Ez(sigma[[i]], B[[i]], zmat[[i]], Ytf) 
  E_theta[[i]] = E_u[[i]] + E_z[[i]]
  
  prob_list[[i]] = P_u[[i]] * P_z[[i]] * prod(as.numeric((sigma[[i]] * zmat[[i]][Ytf,2]) + (1-sigma[[i]] * zmat[[i]][Ytf,1])))
  #prob_list[[i]] = P_u[[i]] * P_z[[i]] * prod( (B[[i]]*zmat[[i]][Ytf,2]) + ((1-B[[i]])*zmat[[i]][Ytf,1])  )
  
  
  
  # verify increase in likelihoods per iteration; also, check if time to break loop
  if(i > 1){
    print(paste("E_u increases each interation:", (E_u[[i]] > E_u[[i-1]])))  # E_u should increase after "u <- u_new"
    print(paste("E_z increases each interation:", (E_z[[i]] > E_z[[i-1]])))  # E_z should increase after "zmat <- zmat_new"
    #if(prob_list[[i]] / prob_list[[i-1]] < 0.000000001) break               # check for exit criteria
  }
  
  
  #
  # M-step
  #
  
  #
  # Set z(i+1) ← argmax_z(Ez) by counting & normalizing
  # 
  print(compute_Ez(sigma[[i]], B[[i]], zmat[[i]], Ytf))
  gamma <- matrix(data = NA, nrow = nrow(zmat[[i]]), ncol = ncol(zmat[[i]]))
  for(obs in seq(1:n_classes_observed)){
    gamma[obs, 1] = max(as.integer(obs==1), .001) + sum(1 - B[[i]][Ytf == obs])
    gamma[obs, 2] = max(as.integer(obs==2), .001) + sum(B[[i]][Ytf == obs])
  }
  gammalist[[i]] <- gamma
  zmat[[i + 1]] <- t(t(gamma)/colSums(gamma))      # normalize the columns to sum to 1
  print(compute_Ez(sigma[[i]], B[[i]], zmat[[i+1]], Ytf))
  z_better[[i]] <- compute_Ez(sigma[[i]], B[[i]], zmat[[i+1]], Ytf) - compute_Ez(sigma[[i]], B[[i]], zmat[[i]], Ytf)
  z_better2[[i]] <- compute_Ez2(sigma[[i]], B[[i]], zmat[[i+1]], Ytf) - compute_Ez2(sigma[[i]], B[[i]], zmat[[i]], Ytf)
  
  
  #
  # Set u(i+1) ← argmax_u(Eu) using glmnet() with weights
  #
  pseudo_dataset <- cbind(Xt, Xt)
  pseudo_labels <- c(rep(1, ncol(Xt)), rep(2, ncol(Xt)))
  weights_f <- c((1-B[[i]]), B[[i]])
  
  # original method, naively picking lambda
  #model <- glmnet(t(pseudo_dataset), as.factor(pseudo_labels), family="binomial", weights = as.numeric(weights_f), alpha = 1, lambda = .0004) 
  #u[[i + 1]] <- as.matrix(model$beta)  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, set to 0.0004 for now (randomly)
  
  # second attempted method using the cv.glmnet function to auto-select lambda
  cvfit = cv.glmnet(t(pseudo_dataset), as.factor(pseudo_labels), family = "binomial", 
                    weights = as.numeric(weights_f), alpha = 1, intercept = FALSE)
  u[[i + 1]] <- as.matrix(coef(cvfit, s = "lambda.min"))[2:(dim(Xt)[1] + 1),]                               # remove the "intercept" term from the result because we aren't fitting this value anyway
  #u_better[[i]] <- compute_Eu(sigma[[i]], B[[i]], u[[i + 1]]) - compute_Eu(sigma[[i]], B[[i]], u[[i]])     # computing u with original sigma/Beta values
  s <- sigmoid(t(u[[i + 1]]) %*% Xt)                                                                        # recompute sigma for use in new E_u 
  b <- compute_B(s, zmat[[i+1]], Ytf)                                                                       # recompute b for use in new E_u
  u_better[[i]] <- compute_Eu(s, b, u[[i + 1]]) - compute_Eu(sigma[[i]], B[[i]], u[[i]])
  
   
}


# if z and u are returning better results at each iteration, the values of z_better and u_better 
# (the difference between Expected Values at each iteration) should always be positive.
# compute the percentage of iterations where the update incrases Expected Value
print("z_better")
unlist(z_better)
sum(unlist(z_better) >= 0)/length(unlist(z_better))
print("u_better")
unlist(u_better)
sum(unlist(u_better) >= 0)/length(unlist(u_better))


# plots should show E_u and E_z per iteration - this DOES NOT necessarily increase each iteration - which is OK
par(mfrow=c(1,2))
plot(unlist(E_u), main = "E_u")
plot(unlist(E_z), main = "E_z")




# DEBUGGING - 
# variables that I can look at...
# u, zmat, B, P_z, P_u
# i think i'm calculating P_z and P_u correctly
# computation of B seems to be correct
# WHY ARE THE MATRIX COLUMNS appearing to SWITCH PLACES in the result??? Is it just ending on a bad esimtation of z?

# QUESTION - the way I have it impelented, the gamma_jk version of E_z isn't the same as the original E_z (what is the gamma there??)


# for(i in seq(2:10)){
#   print("zmat:")
#   print(zmat[[i]])
#   print("Ez")
#   print(E_z[[i]])
#   print("gammas")
#   print(gammalist[[i]])
# }
