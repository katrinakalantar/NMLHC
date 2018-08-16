library(GEOquery)
library(gplots)
library(RColorBrewer)
library(mBALPkg)
library(RCurl)
library(Rfast)
library(scales)
library(SimSeq)
library(glmnet)


split_train_test <- function(input){
  s <- shuffle(seq(1:dim(input)[2]))
  training_set <- input[,s[1:round(dim(input)[2]*.7)]]
  test_set <- input[,s[(round(dim(input)[2]*.7)+1):dim(input)[2]]]
  return(list("training_set" = training_set, "test_set" = test_set))
}

subset_known_dataset <- function(geo_data, pos_regex, neg_regex){
  
  pos_split <- split_train_test(geo_data[, grep(pos_regex, geo_data$source_name_ch1)])
  neg_split <- split_train_test( geo_data[, grep(neg_regex, geo_data$source_name_ch1)])
  
  full_train <- combine(pos_split$training_set, neg_split$training_set)
  full_test <- combine(pos_split$test_set, neg_split$test_set)
  
  true_labels_train <- rep(0, length(full_train$source_name_ch1))
  true_labels_train[grep(pos_regex, full_train$source_name_ch1)] <- 1
  
  true_labels_test <- rep(0, length(full_test$source_name_ch1))
  true_labels_test[grep(pos_regex, full_test$source_name_ch1)] <- 1
  
  return(list("train_data" = full_train, "test_data" = full_test, "train_labels" = true_labels_train, "test_labels" = true_labels_test))
  
}

project_pca <- function(input_dataset, c){
  p <- prcomp(t(exprs(input_dataset$train_data)))
  plot(p$x[,1],p$x[,2], 
       col=alpha(c[as.integer(input_dataset$train_labels == 1) + 1], .5), 
       xlab="PC1",ylab="PC2", cex=1.6, pch=16)
  
  s <- scale(t(exprs(input_dataset$test_data)), p$center, p$scale) %*% p$rotation 
  points(s[,1],s[,2], col=alpha(c[as.integer(input_dataset$test_labels == 1) + 1], .8), cex=1.6)
}


run_scramble_classification <- function(dataset, num_iterations, most_var, glm_alpha){
  
  master_train_res <- list()
  master_test_res <- list()
  master_count_res <- list()
  
  for(count in seq(1:num_iterations)){
    
    # NOTE: THIS LOGIC REQUIRES THAT ALL POSITIVES are first, followed by ALL NEGATIVE cases
    all_labels <- list()
    for(i in seq(0,10)){
      all_labels[[i+1]] <- list()
      for(j in seq(0,10)){
        this_swap <- dataset$train_labels
        
        # swap i A's to B's  
        originally_positive <- this_swap[this_swap == 1]
        originally_positive[sample(seq(1:length(originally_positive)),  round((((i*10)/2)/100)*length(originally_positive)))] <- 0
        # swap j B's to A's
        originally_negative <- this_swap[this_swap == 0]
        originally_negative[sample(seq(1:length(originally_negative)), round((((j*10)/2)/100)*length(originally_negative)))] <- 1
        
        # recombine data
        this_swap <- c(originally_positive, originally_negative)
        print('final')
        print(this_swap)
        all_labels[[i+1]][[j+1]] <- this_swap
      }
    }
    
    mm_train <- matrix(0, length(seq(0,10)), length(seq(0,10)))
    counts_train <- matrix(0, length(seq(0,10)), length(seq(0,10)))
    mm_test <- matrix(0, length(seq(0,10)), length(seq(0,10)))
    
    for(i in seq(0,10)){
      for(j in seq(0,10)){
        
        print(paste("i:",i))
        print(paste("j:",j))
        
        rv <- rowVars(exprs(dataset$train_data))
        names(rv) <- rownames(exprs(dataset$train_data))
        #merged <- merged[names(tail(sort(rv),1000)),]
        
        # build model for each label set in list
        model <- glmnet( t(exprs(dataset$train_data[names(tail(sort(rv),most_var)),])), as.integer(as.factor(all_labels[[i+1]][[j+1]])), family = c("binomial"), alpha = glm_alpha, lambda = .6 )
        counts_train[i+1, j+1] <- sum(abs(model$beta) > 0)
        
        # evaluate performance on training set
        p_train <- predict(model, t(exprs(dataset$train_data[names(tail(sort(rv),most_var)),])))
        x <- pROC::roc(dataset$train_labels,  p_train[,1], ci = TRUE)
        print(x$auc[[1]])
        mm_train[i+1,j+1] <- x$auc[[1]]
        
        # evaluate performance on test set
        p_test <- predict(model, t(exprs(dataset$test_data[names(tail(sort(rv),most_var)),])))
        x <- pROC::roc(dataset$test_labels, p_test, ci = TRUE)
        print(x$auc[[1]])
        mm_test[i+1,j+1] <- x$auc[[1]]
        
      }
    }
    
    colnames(mm_train) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    rownames(mm_train) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    colnames(mm_test) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    rownames(mm_test) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    colnames(counts_train) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    rownames(counts_train) <- lapply(seq(0,10), function(x){return(paste((x*10)/2, "%", sep=""))})
    
    mm_train <- round(mm_train,2)
    mm_test <- round(mm_test,2)
    
    master_train_res[[count]] <- mm_train
    master_test_res[[count]] <- mm_test
    master_count_res[[count]] <- counts_train
    
  }
  
  train_res_mean <- apply(simplify2array(master_train_res), 1:2, mean)
  train_res_sd <- apply(simplify2array(master_train_res), 1:2, sd)
  test_res_mean <- apply(simplify2array(master_test_res), 1:2, mean)
  test_res_sd <- apply(simplify2array(master_test_res), 1:2, sd)
  counts_train_res_mean <- apply(simplify2array(master_count_res), 1:2, mean)
  counts_train_res_sd <- apply(simplify2array(master_count_res), 1:2, sd)
  
  return(list("training_result" = train_res_mean, "test_result" = test_res_mean, "counts_model" = counts_train_res_mean))
}


plot_res_heatmaps <- function(results){
  
  my_palette <- colorRampPalette(c("purple3", "white"))(n = 100)
  my_palette2 <- colorRampPalette(c("white", "green4"))(n = 1001)
  my_palette3 <- colorRampPalette(c("blue", "white","red"))(n=100)
  
  heatmap.2(results$training_result, Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette,sepcolor="black",
            colsep=0:ncol(results$training_result)+1,rowsep=0:nrow(results$training_result)+1, sepwidth=c(0.01,0.01), main = "Train AUC",
            cellnote=results$training_result, notecol='black')
  heatmap.2(results$test_result, Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette,sepcolor="black",
            colsep=0:ncol(results$test_result)+1,rowsep=0:nrow(results$test_result)+1, sepwidth=c(0.01,0.01), main = "Test AUC",
            cellnote=results$test_result, notecol='black')
  heatmap.2(round(results$training_result-results$test_result,2), Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette3,sepcolor="black",
            colsep=0:ncol(results$training_result)+1,rowsep=0:nrow(results$training_result)+1, sepwidth=c(0.01,0.01), main = "Diff AUC (Train - Test)",
            cellnote=round(results$training_result-results$test_result,2), notecol='black')
  heatmap.2(results$counts_model, Rowv=FALSE, trace="none", Colv=FALSE, col=my_palette2,sepcolor="black",
            colsep=0:ncol(results$counts_model)+1,rowsep=0:nrow(results$counts_model)+1, sepwidth=c(0.01,0.01), main="Counts",
            cellnote=results$counts_model, notecol='black')
  
}


generate_simulated_data <- function(negative_reference_data, positive_reference_data, n = 1, DEgenes){
  new_negatives <- c()
  new_positives <- c()
  
  min_dim <- min(dim(negative_reference_data)[2], dim(positive_reference_data)[2])/2
  
  for(i in seq(1:n)){
    sd1 <- SimData(cbind(exprs(negative_reference_data), exprs(positive_reference_data)), 
                   treatment = c(rep(0,as.integer(dim(negative_reference_data)[2])),rep(1,as.integer(dim(positive_reference_data)[2]))), 
                   genes.select = rep(TRUE,length(rownames(exprs(positive_reference_data)))),
                   genes.diff = (rownames(exprs(positive_reference_data)) %in% sample(DEgenes$V1, 800)),
                   n.diff = 800, k.ind=min_dim, sort.method="unpaired", switch.trt=1)
    sd2 <- SimData(cbind(exprs(negative_reference_data), exprs(positive_reference_data)), 
                   treatment = c(rep(0,as.integer(dim(negative_reference_data)[2])),rep(1,as.integer(dim(positive_reference_data)[2]))),
                   genes.select = rep(TRUE,length(rownames(exprs(positive_reference_data)))),
                   genes.diff = (rownames(exprs(positive_reference_data)) %in% sample(DEgenes$V1, 800)),
                   n.diff = 800, k.ind=min_dim, sort.method="unpaired", switch.trt=0)
    if(i == 1){
      new_negatives <- sd1$counts[,1:min_dim] 
      new_positives <- sd2$counts[,(min_dim + 1):(min_dim*2)]
    }else{
      new_negatives <- cbind(new_negatives, sd1$counts[,1:min_dim] )
      new_positives <- cbind(new_positives,sd2$counts[,(min_dim + 1):(min_dim*2)])    
    }
  }
  
  return_dataset <- cbind(new_positives, new_negatives)
  colnames(return_dataset) <- c(lapply(seq(1:ncol(new_positives)), function(x){return(paste("Sim",x,"_Group","1",sep=""))}), lapply(seq(1:ncol(new_negatives)), function(x){return(paste("Sim",x,"_Group","0",sep=""))}))
  return(return_dataset)
}

make_sim_eset <- function(input_data){
  pd <- cbind(as.integer(grepl("Group1", colnames(input_data))), rep("sim", length(colnames(input_data))))
  rownames(pd) <- colnames(input_data)
  colnames(pd) <- c("classification","simulated_status")
  return(ExpressionSet(input_data, phenoData = AnnotatedDataFrame(as.data.frame(pd))))
}

generate_simulated_dataset <- function(input_data, training_names, DEgenes){
  training_set <- input_data[,training_names]
  test_set <- input_data[,!(colnames(input_data) %in% training_names)]
  
  training_dataset <- generate_simulated_data(training_set[,training_set$effective_group == 4],  training_set[,training_set$effective_group == 1], 2, DEgenes)
  test_dataset <- generate_simulated_data(test_set[,test_set$effective_group == 4], test_set[,test_set$effective_group == 1], 3, DEgenes)
  
  training_dataset_normalized <- make_sim_eset(t(t(training_dataset)/colSums(training_dataset)))
  test_dataset_normalized <- make_sim_eset(t(t(test_dataset)/colSums(test_dataset)))
  
  return(list("train_data" = training_dataset_normalized, "test_data" = test_dataset_normalized, "train_labels" = training_dataset_normalized$classification, "test_labels" = test_dataset_normalized$classification))
}