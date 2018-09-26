
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
# Set-up Data - Iris dataset
#

data(iris)
iris_mod <- iris[sample(nrow(iris)),]
keep <- c("setosa","virginica")
iris_mod <- iris_mod[iris_mod$Species %in% keep,]

train <- sample(seq(1,dim(iris_mod)[1]),.7*dim(iris_mod)[1])  # create dummy training / test dataset

Xt <- t(data.matrix(iris_mod[train,1:4]))
Yt <- as.character(iris_mod[train,5])
Yt[seq(1,round(length(Yt)/4))] = c("unlikely","maybe")[as.integer(Yt[seq(1,round(length(Yt)/4))] == "setosa")+1]  # convert setosa to maybe, verginica to unlikely
# NOTE: order the levels in accordance with the "true" labels
# (ie 1 == negative, 2 == positive, 3... are additional observed classes that do not correspond to a true latent label)
Ytf <- as.integer(factor(Yt, levels = c("setosa","virginica","maybe","unlikely")))

training_dist <- table(Ytf)
names(training_dist) <- c("setosa","virginica","maybe","unlikely")
training_dist/length(Ytf)
true_dist <- table(iris_mod$Species)/length(iris_mod$Species)

Xs <- t(data.matrix(iris_mod[!(rownames(iris_mod) %in% train),1:4]))
Ys <- as.character(iris_mod[!(rownames(iris_mod) %in% train),5])
Ys[seq(1,round(length(Ys)/4))] = c("unlikely","maybe")[as.integer(Ys[seq(1,round(length(Ys)/4))] == "setosa")+1]  # convert setosa to maybe, verginica to unlikely
Ysf <- as.integer(factor(Ys), levels = c("setosa","virginica","maybe","unlikely"))


#
# Run EM
#

# algorithm parameters
n_classes_observed = 4                                  # number of observed classes (clinician labels)
n_true_classes = 2                                      # number of true classes
max_iter = 10                                        # maximum number of iterations of EM before quitting with "did not converge"

### EM Algorithm Implementation

# Set theta_1 to some sensible values, where theta_1 = (Z, u)
m <- matrix(runif(n_true_classes * n_classes_observed), nrow=n_classes_observed)    # initialize Z probability parameter matrix
#zmat_orig <- m/rowSums(m)
zmat_orig <- t(t(m)/colSums(m))      # normalize the columns to sum to 1
u_orig <- as.matrix(runif(nrow(Xs), 0, 1))                                          # initialize u weight vector

# create list variables to store values on every iteration so that I can track errors
probs <- list(); prob_list <- list();
sigma <- list(); B <- list(); zmat <- list(); u <- list()
P_u <- list(); P_z <- list(); E_u <- list(); E_z <- list(); E_theta <- list()
continue = TRUE

#
# Functions definitiones that are used multiple times
#

compute_Pu <- function(u){
  return(prod(exp(-abs(u))))
}
# THERE IS SOMETHING WRONG WITH THIS FUNCTION
compute_Pz <- function(z){
  return(((1-sum(z[,1])) * prod(z[,1] ^ c(1,0,0,0))) * ((1-sum(z[,2])) * prod(z[,2] ^ c(0,1,0,0))))
}
compute_Eu <- function(sigma, B, u){
  P_u <- compute_Pu(u)
  E_u <- log(P_u) + rowSums((1-B)*log(1-sigma) + B*log(sigma))
  return(E_u)
}
compute_Ez <- function(sigma, B, z, Ytf){  # NOTE: is it OK to have added one (1) to all of these log values??
  P_z <- compute_Pz(z)
  E_z <- log(P_z + 1) + rowSums((1-B)*log(z[,1][Ytf] + 1) + B*log(z[,2][Ytf] + 1)) 
  #E_z <- log(P_z) + rowSums((1-B)*log(z[,1][Ytf]) + B*log(z[,2][Ytf])) 
  return(E_z)
}

# compute E_u and E_z for the initial values.
u[[1]] = u_orig
zmat[[1]] = zmat_orig
sigma[[1]] <- sigmoid(t(u[[1]]) %*% Xt)       
B[[1]] = 1 / (1 + (((1-sigma[[1]]) * zmat[[1]][Ytf,1]) / (sigma[[1]] * zmat[[1]][Ytf,2])))   # lines in paper: P(bn=1|theta, xn, cn)
P_u[[1]] = compute_Pu(u[[1]])
P_z[[1]] = compute_Pz(zmat[[1]])
E_u[[1]] = compute_Eu(sigma[[1]], B[[1]], u[[1]])
E_z[[1]] = compute_Ez(sigma[[1]], B[[1]], zmat[[1]], Ytf) 
E_theta[[1]] = E_u[[1]] + E_z[[1]]
prob_list[[1]] = P_u[[1]] * P_z[[1]] * prod(as.numeric((sigma[[1]] * zmat[[1]][,2][Ytf]) + (1-sigma[[1]] * zmat[[1]][,1][Ytf])))

#print(paste("DEBUG: dim(sigma) = ", dim(sigma[[1]])[1], "x", dim(sigma[[1]])[2]))
#print(paste("DEBUG: dim(B) = ", dim(B[[1]])[1], "x", dim(B[[1]])[2]))

for(i in seq(1:max_iter)){
  print(paste("------- iteration:", i, " -------"))
  
  ### Set z(i+1) ← argmaxz(Ez) by counting & normalizing
  print(zmat[[i]])
  print(rowSums(zmat[[i]]))
  gamma <- matrix(data = NA, nrow = nrow(zmat[[i]]), ncol = ncol(zmat[[i]]))
  for(obs in seq(1:n_classes_observed)){
    gamma[obs, 1] = as.integer(obs==1) + sum(1 - B[[i]][Ytf == obs])
    gamma[obs, 2] = as.integer(obs==2) + sum(B[[i]][Ytf == obs])
  }
  zmat_new <- t(t(gamma)/colSums(gamma))      # normalize the columns to sum to 1
  #zmat_new <- zmat_new/rowSums(zmat_new)      # normalize rows?? NOT SURE IF THIS IS CORRECT
  #zmat_new <- gamma/rowSums(gamma)      # normalize the columns to sum to 1
  print(zmat_new)
  zmat[[i + 1]] = zmat_new
  
  P_z[[i + 1]] = compute_Pz(zmat[[i + 1]])
  E_z[[i + 1]] = compute_Ez(sigma[[i]], B[[i]], zmat[[i + 1]], Ytf)  # compute Ez after updating zmat
  
  ### Set u(i+1) ← argmaxu(Eu) using glmnet() with weights
  pseudo_dataset <- cbind(Xt, Xt)
  pseudo_labels <- cbind(rep(0, ncol(Xt)), rep(1, ncol(Xt)))
  weights_f <- cbind((1-B[[i]]), B[[i]])
  model <- glmnet(t(pseudo_dataset), as.numeric(pseudo_labels), family="binomial", weights = as.numeric(weights_f), alpha = 1) 
  u[[i + 1]] <- model$beta[,ncol(model$beta)-3]  # NOTE: SHOULD GET LAMBDA PARAM BY CROSS-VALIDATION, just selected third-to-last lambda row randomly for now
  
  sigma[[i + 1]] <- sigmoid(t(u[[i + 1]]) %*% Xt)
  B[[i + 1]] <- 1 / (1 + (((1-sigma[[i + 1]]) * zmat[[i + 1]][Ytf,1]) / (sigma[[i + 1]] * zmat[[i + 1]][Ytf,2]))) 
  P_u[[i + 1]] = compute_Pu(u[[i+1]])
  E_u[[i + 1]] = compute_Eu(sigma[[i+1]], B[[i + 1]], u[[i+1]])  # compute Eu after updating u
  
  print(paste("E_u increases each interation:", (E_u[[i+1]] > E_u[[i]]))) # E_u should increase after "u <- u_new"
  print(paste("E_z increases each interation:", (E_z[[i+1]] > E_z[[i]]))) # E_z should increase after "zmat <- zmat_new"
  
  prob_list[[i + 1]] = P_u[[i+1]] * P_z[[i+1]] * prod(as.numeric((sigma[[i+1]] * zmat[[i+1]][,2][Ytf]) + (1-sigma[[i+1]] * zmat[[i+1]][,1][Ytf])))
  
  E_theta[[i + 1]] = E_u[[i + 1]] + E_z[[i + 1]]
  
  probs[[i + 1]] <- prob_list[[i+1]]/prob_list[[i]]
  
  # the error seems to go up/down wildly in the initial rounds (1-10), so require at least 10 rounds for now
  #if( i>2 & abs(new_prob)/abs(previous_prob) < .0001) break
  #if( abs(new_prob)/abs(previous_prob) < .0001) break
  #if(prob_list[[i + 1]] / prob_list[[i]] < 0.000000001) break
  
}

plot(seq(1:length(probs)), probs[1:length(probs)], ylim = c(-10,10), xlab = "EM Iteration", ylab = "prob")
print(zmat)