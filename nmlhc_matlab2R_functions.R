

convert_genenames <- function(data, gmt_file){
  
  genenames <- colnames(data)
  
  # # create the genemap variable, so we don't rely on biomaRt
  # 
  # #ran this section of code once to generate the HSapiens_gene_ensembl.csv file
  # listMarts(host="www.ensembl.org")
  # mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
  # genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
  #                     filters = "ensembl_gene_id",
  #                     values = colnames(x_mod),
  #                     mart)
  # write.csv(genemap,"/Users/kkalantar/Documents/Research/NMLHC/reference/HSapiens_gene_ensembl.csv",quote=FALSE,row.names=FALSE)
  
  genemap <- read.csv("/Users/kkalantar/Documents/Research/NMLHC/reference/HSapiens_gene_ensembl.csv", header = TRUE)
  
  # create an index in which you are rearranging ensmbl gene ids from genemap into same order as they are in genenames
  idx <- match(genenames, genemap$ensembl_gene_id)  
  hgnc_symbol <- genemap$hgnc_symbol[idx]
  
  # convert original genecounts matrix, filtered for PC-genes, to have HGNC rownames (necessary to run CIBERSORT)
  pc_hgnc_genecounts = data
  colnames(pc_hgnc_genecounts) = hgnc_symbol
  pc_hgnc_genecounts <- pc_hgnc_genecounts[,!(colnames(pc_hgnc_genecounts) %in% c(""))]
  
  return(pc_hgnc_genecounts)
}

collapse_pathways <- function(data, gmt_file){
  pathway_sums <- lapply(gmt_file$genesets, function(x){return(rowSums(data[,colnames(data) %in% unlist(x)]))})
  names(pathway_sums) <- gmt_file$geneset.names
  final_matrix <- do.call(cbind, pathway_sums)
  return(final_matrix)
}

collapse_pathways2 <- function(data, gmt_file, keep_pathways){
  print("In collapse_pathways2")
  #names(gmt_file$genesets) <- gmt_file$geneset.names
  #pathway_genes <- unique(unlist(lapply(gmt_file$genesets, function(x){return(unlist(x))})))
  pathways_genes <- as.vector(unique(unlist(lapply(seq(1:length(gmt_file$genesets)), function(x){if(gmt_file$geneset.names[x] %in% keep_pathways){return(gmt_file$genesets[x])}}))))
  print(length(pathways_genes))
  return(data[,colnames(data) %in% pathways_genes])
}


plot_confusion_matrix <- function(CM, colorpal = colorRampPalette(c("purple3", "white"))(n = 201), max=300 ){
  a <- CM[1,1,1,1,,]
  gl <- 0
  plot_list <- list()
  for(i in colnames(a)){
    for(j in rownames(a)){
      gl <- gl + 1
      ints <- as.numeric(strsplit(a[j,i],"_")[[1]])
      #m <- matrix(ints/sum(ints), ncol = 2)
      m <- matrix(ints, ncol=2)
      #m <- t(t(m)/colSums(m))
      colnames(m) <- c("REF_0", "REF_1")
      rownames(m) <- c("PRED_0", "PRED_1")
      hmp <- pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,col=colorpal,
                      border_color = "black", fontsize_number = 10, cellheight=25, cellwidth=25, breaks = seq(0,max,by=1))
      plot_list[[gl]] <- hmp[[4]]
    }
  }
  
  #par(mar=c(5,15,5,5))
  g <- do.call(grid.arrange, plot_list)
}


#' Generate predictions based on an ensemble of classification algorithms; 
#' 
#' @param x training data; matrix of features/covariates used for prediction.
#' @param y training labels; must be factor of non-integer values i.e. "one"/"two" instead of 1/2
#' @param library A vector of strings indicating the algorithms to be included in the ensemble. Available algorithms can be queried here: https://topepo.github.io/caret/available-models.html
#' @param multiple Boolean value indicating whether to run multiple iterations of cross-validation to generate ensemble predictions; default = TRUE, with 5 iterations 
#' @return A list containing: "MF", samples to be kept by majority filter; "CF", samples to be kept by consensus filter, "full_res", the full list of discordant predictions
#' @examples
#' make_ensemble(x, y, c("regLogistic", "rf", "knn", "svmLinear3", "nnet"), multiple = TRUE)
make_ensemble <- function(x, y, library, multiple = TRUE){
  list_of_results <- list()
  for(algorithmL in library){
    ctrl <- trainControl(method = "cv", savePred=T) #, classProb=T
    
    if(multiple){
      print("generating ensemble with repeats")
      ctrl <- trainControl(method = "repeatedcv", savePred=T, number = 10, repeats = 5) #classProb=T, 
    }
    
    mod <- train(x, y, method = algorithmL, trControl = ctrl)
    all_predictions <- mod$pred
    predictions_to_use <- all_predictions[rowSums(do.call(cbind, lapply(names(mod$bestTune), function(x){print(x); return(all_predictions[,x] == mod$bestTune[[x]])}))) == length(mod$bestTune),]
    predictions_to_use <- predictions_to_use[order(predictions_to_use$rowIndex),]   # MAJOR FIX 10/15
    
    if(multiple){
      # collapse predictions_to_use
      n <- do.call(rbind, lapply(unique(predictions_to_use$rowIndex), function(x){b = predictions_to_use[predictions_to_use$rowIndex == x,]; a = which.max(table(b$pred)); print(names(a)); return(b[b$pred == names(a), ])}))
      n <- n[,c("pred","obs","rowIndex")]
      n <- n[!duplicated(n), ]
      predictions_to_use <- n
    }
    
    concordant <- predictions_to_use$pred == predictions_to_use$obs
    accuracy <- sum(predictions_to_use$pred == predictions_to_use$obs)/length(predictions_to_use$obs)
    logger.info(msg = sprintf("MAIN - ALGO - ENSEMBLE - %s, accuracy = %f", algorithmL, accuracy))

    discordant <- predictions_to_use$pred != predictions_to_use$obs
    
    list_of_results[[algorithmL]] = discordant
  }
  
  r <- do.call(cbind, list_of_results)
  
  majority_filter <- rowSums(r)/ncol(r) > .5
  consensus_filter <- rowSums(r)/ncol(r) == 1
  
  keep_majority_filter <- !majority_filter
  keep_consensus_filter <- !consensus_filter
  
  return(list("MF" = keep_majority_filter, "CF" = keep_consensus_filter, "full_res" = r) )#, "P" = p))
}

#' Evaluate the performance of the ensemble for identifying known flipped samples
#' 
#' @param keeping array of TRUE/FALSE values indicating whether a sample was kept after filtering
#' @param fdz array of values -1/1 indicating whether a sample has its true label (-1) or was flipped (1)
#' @return A confusionMatrix object
#' @examples
#' flag_flipped_samples(filter$MF, shuffled_fdz)
flag_flipped_samples <- function(keeping, fdz){
  cm <- confusionMatrix(factor(as.integer(keeping)), factor(as.integer(fdz < 0)), positive = "1")
  print(cm)
  plot(cm$table)
  return(cm)
}

make_superlearner <- function(x, y, v, library){
  
  cv_sl = CV.SuperLearner(Y = y, X = x, family = binomial(), V = v,
                          method = "method.AUC",
                          SL.library =  library) #
  
  # create the barplot of learners.
  a <- summary(cv_sl)
  aa <- a$Table
  print(aa)
  min_performance <- mean(aa$Ave) - 2*sd(aa[c(3:nrow(aa)),"Ave"])
  max_performance <- mean(aa$Ave) + 2*sd(aa[c(3:nrow(aa)),"Ave"])
  aa$col <- as.factor(aa$Ave > min_performance)
  
  p<-ggplot(data=aa, aes(x=Algorithm, y=Ave, fill = col )) + ylim(c(0,1)) +
    geom_bar(stat="identity") + scale_fill_brewer(palette="Paired")  + theme_minimal() + ggtitle("AUC for Each SL") +
    geom_abline(slope=0, intercept=0.5,  col = "red",lty=2) +
    geom_abline(slope=0, intercept=min_performance,  col = "gold",lty=2) +
    geom_abline(slope=0, intercept=mean(aa$Ave),  col = "gold",lty=1) +
    geom_abline(slope=0, intercept=max_performance,  col = "gold",lty=2) 
  print(p)
  
  individual_predictions <- cv_sl$library.predict
  print(individual_predictions)
  binary_predictions <- individual_predictions > .5
  print(binary_predictions)
  concordant <- (binary_predictions+1) == as.numeric(y)
  discordant <- (binary_predictions+1) != as.numeric(y)
  
  barplot(table(rowSums(discordant)/ncol(discordant))) # here 1 is ConsensusFiltering, > .5 is Majority Filtering
  abline(v = 2.5, col="red", lty=2)
  
  majority_filter <- rowSums(discordant)/ncol(discordant) > .5
  consensus_filter <- rowSums(discordant)/ncol(discordant) == 1
  
  keep_majority_filter <- !majority_filter
  keep_consensus_filter <- !consensus_filter
  
  return(list("MF" = keep_majority_filter, "CF" = keep_consensus_filter, "P" = p))
}


plot_multiple_learning_curves <- function(list_of_LCs, colors_to_use, learning_curve_iters, success, significance_level = .05, plot_all_points = TRUE, plot_error_bars = FALSE, title="Learning Curve"){
  index <- 1
  for(LC in list_of_LCs){
    
    if(index == 1){
      plot(colMeans(LC[1,success,]), pch = 16, col = colors_to_use[index], ylim = c(0,1), xlim = c(-5, 20), main = title, 
           ylab = "AUC", cex = 1.5, xaxt="n")
      axis(1, at = 1:length(learning_curve_iters), labels = learning_curve_iters)
    }else{
      points(colMeans(LC[1,success,]), pch = 16, col = colors_to_use[index], cex = 1.5)
    }
    
    lines(colMeans(LC[1,success,]), pch = 16, col = colors_to_use[index])
    text(-3, colMeans(LC[1,success,])[1], names(list_of_LCs)[index], cex=.8, col = colors_to_use[index])
    
    if(plot_error_bars){
      sd <- apply(LC[1,success,], 2, sd)
      arrows(seq(1:length(learning_curve_iters)), colMeans(LC[1,success,]) - sd, seq(1:length(learning_curve_iters)), colMeans(auc_full_filt[1,success,]) + sd, length=0.05, angle=90, code=3, col="gray28")
    }
    if(plot_all_points){
      matplot(t(LC[1,success,]),type="p", pch=16, col=alpha(colors_to_use[index],.4), add = TRUE, jitter=1)
    }
    index <- index + 1 
  }
  
  if(length(list_of_LCs) == 2){
    significance <- unlist(lapply(seq(1:length(list_of_LCs[[1]][1,success,1])), function(x){return(wilcox.test(list_of_LCs[[1]][1,success,x], list_of_LCs[[2]][1,success,x])$p.value)}))
    print(significance)
    points(seq(1:length(significance)), as.integer(significance < significance_level)*.1, pch = 8, col=c("white","black")[as.integer(significance < significance_level)+1])
    text(-3, .1, paste("p <",significance_level), cex = .8, col = "black")
  }
}


#' Run feature selection
#' 
#' @param x training data; matrix of features/covariates used for prediction.
#' @param y training labels.
#' @param method The feature selection method used to rank features, per rankfeatures() documentation in matlab (options: "ttest", "entropy", "bhattacharyya", "roc", "wilcoxon")
#' @return The sorted matrix of features.
#' @examples
#' run_FS(x, y, "wilcoxon")
run_FS <- function(x, y, method){
  # method options are "ttest", "entropy", "bhattacharyya", "roc", "wilcoxon"
  setVariable(matlab, x = x)
  setVariable(matlab, y = y)
  setVariable(matlab, method=method)
  evaluate(matlab, "[idx, z] = rankfeatures(x', y, 'Criterion', method);")
  feature_values <- getVariable(matlab, "z")[[1]]
  rownames(feature_values) <- colnames(x)
  ranks <- getVariable(matlab, "idx")[[1]]  
  result <- as.matrix(feature_values[colnames(x)[ranks],])
  return(result)
}


#' Plot the jaccard similarity (percentage of overlapping features) for models built using data with and without mislabelled samples
#' 
#' @param features_selected A list containing the features selected at each iteration - contains list entries for "unflipped" and "flipped_wilcox".
#' @param check_values Vector containing integers specifying the top N genes that should be checked for overlap percentage between flipped and unflipped.
#' @param fileroot Optional parameter specifying the file name to be used for saving .pdf version of the plot.
#' @return The matrix of overlap percentages computed at each check_value(s)
#' @examples
#' overlap <- plot_feature_similarity(features_selected, c(100, 500, 1000, 2500, 5000, 10000))#, "fileroot" = paste(EXPERIMENT_DIR, "features_byflip", sep=""))
plot_feature_similarity <- function(features_selected, check_values, fileroot = NULL){
  for(d in features_selected){
    dims = dim(d$unflipped)
    discrete_pal2 <- colorRampPalette(c("white", "dodgerblue4"))(n = dims[2] * dims[3])
    
    overlap_mat <- array(rep(length(check_values) * 0,dims[2] * dims[3]), c(length(check_values), dims[2], dims[3]))
    for(j in seq(1:dims[2])){
      for(i in seq(1:dims[3])){
        mean_jaccard <- c()
        for(N in seq(1:length(check_values))){
          jaccard <- c();
          for(k in seq(1:dim(d$flipped_wilcox)[1])){
            jaccard <- c(jaccard, length( intersect( head(d$unflipped[1,j,i,], n = check_values[N]), 
                                                     head(d$flipped_wilcox[1,j,i,], n = check_values[N]) )) / check_values[N])
          }
          print(mean(jaccard))
          mean_jaccard <- c(mean_jaccard, jaccard) 
          overlap_mat[N, j, i] <- mean(jaccard)
        }
      }
    }
  }
  
  if(!is.null(fileroot)){
    pdf(paste(fileroot, ".pdf", sep=""), height = 6, width= 6)
  }
  
  plot(overlap[,1,1], ylim = c(0,1), col="black", lwd = 4, cex= 2, xlim = c(-1,length(overlap[,1,1])+2),
       ylab = "% feature overlap - Jaccard Similarity",xaxt = "n", xlab = "overlap considered for top N features")
  axis(1, at = 1:length(overlap[,1,1]), labels = check_values)
  cnt = 0
  for(i in seq(1:dims[3])){
    for(j in seq(1:dims[2])){
      cnt = cnt + 1
      points(overlap[,j,i], col = alpha(discrete_pal2[cnt],.8), pch= 16, cex = 2)
      points(overlap[,j,i], col = alpha(discrete_pal2[cnt],1), cex = 2)
      lines(overlap[,j,i], col = alpha(discrete_pal2[cnt],1), cex = 2)
      text(length(overlap[,j,i]) + 1, overlap[,j,i][length(overlap[,j,i])], paste(c("i: ", EXP_RANGE[i], " j: ",EXP_RANGE_J[j], ", ", overlap[,j,i][length(overlap[,j,i])] * 100, "%"),collapse=""), cex = .5)
      text(0, overlap[,j,i][1], paste(c("i: ", EXP_RANGE[i], " j: ",EXP_RANGE_J[j], ", ", overlap[,j,i][1] * 100, "%"),collapse=""), cex = .5)
    }
  }
  if(!is.null(fileroot)){
    dev.off()
  }
  return(overlap_mat)
}


#' Simulate data from a known dataset (using SimSeq package)
#' 
#' @param original_dataset_name the string identifier for the pre-existing dataset from which data should be sampled. Note: only particular datasets are available/pre-programmed for this method, so string names must match one of the following: "mBAL"
#' @param n_samples the number of iterations of re-sampling to perform; SimSeq will only sample the number of samples in each class of the training dataset, but to increase total samples size, we can iterate running SimSeq.
#' @return A list containing {"x": training data, "y": training class ID, "ff": zeros, "xx": test data, "tt": test class ID, "dd": zeros}.
#' @examples
#' simulate_data("mBAL", 15)
simulate_data <- function(original_dataset_name, n_samples ){
  training_dataset_normalized <- NULL
  test_dataset_normalized <- NULL
  
  if(original_dataset_name == "mBAL"){
    input_data <- filtered_eset     # from mBALPkg
    TRAINING_NAMES <- c("TA.212","TA.225","TA.298","TA.304","TA.314","TA.315","TA.335","TA.337","TA.343","TA.350", 
                        "TA.349","TA.273","TA.331","TA.221","TA.220","TA.215","TA.270","TA.241","TA.211","TA.218")  # same as in mBAL study    
    DEgenes <- read.table("Documents/Research/NMLHC/Exp1_HostBenchmark/data/DEgenes.csv")                           # taken from mBAL study
    
    training_set <- input_data[,TRAINING_NAMES]
    test_set <- input_data[,!(colnames(input_data) %in% TRAINING_NAMES)]
    training_dataset <- generate_simulated_data(training_set[,training_set$effective_group == 4],  
                                                training_set[,training_set$effective_group == 1], n_samples, DEgenes)
    test_dataset <- generate_simulated_data(test_set[,test_set$effective_group == 4], 
                                            test_set[,test_set$effective_group == 1], n_samples, DEgenes)
    training_dataset_normalized <- make_sim_eset(t(t(training_dataset)/colSums(training_dataset)))    # TSS normalize the simulated data
    test_dataset_normalized <- make_sim_eset(t(t(test_dataset)/colSums(test_dataset)))
  }
  
  return(list("x" = t(as.matrix(exprs(training_dataset_normalized))), "y" = as.matrix(as.numeric(training_dataset_normalized$classification == 1) + 1), "ff" = as.matrix(rep(0, length(training_dataset_normalized$classification))),
              "xx" = t(as.matrix(exprs(test_dataset_normalized))), "tt" = as.matrix(as.numeric(test_dataset_normalized$classification == 1) + 1), "dd" = as.matrix(rep(0, length(test_dataset_normalized$classification)))))
}



#' Function called by simulated_data to run SimSeq
#' 
#' @param negative_reference_data ExpressionSet object corresponding to the positive class data.
#' @param positive_reference_data ExpressionSet object corresponding to thenegative class data.
#' @param n The number of iterations of SimSeq sampling to perform.
#' @param DEgenes A list of genes passed to SimSeq which should be differentially expressed in the simulated data.
#' @return An ExpressionSet object containing simulated positive and negative samples, with metadata column "classification" indicating which group they belong to.
#' @examples
#' generate_simulated_data(training_set[,training_set$effective_group == 4], training_set[,training_set$effective_group == 1], 15, DEgenes)
generate_simulated_data <- function(negative_reference_data, positive_reference_data, n = 5, DEgenes){
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
      new_positives <- cbind(new_positives, sd2$counts[,(min_dim + 1):(min_dim*2)])    
    }
  }
  
  return_dataset <- cbind(new_positives, new_negatives)
  colnames(return_dataset) <- c(lapply(seq(1:ncol(new_positives)), function(x){return(paste("Sim",x,"_Group","1",sep=""))}), lapply(seq(1:ncol(new_negatives)), function(x){return(paste("Sim",x,"_Group","0",sep=""))}))
  return(return_dataset)
}


#' Convert matrix of counts with column names to an ExpressionSet with metadata
#' 
#' @param input_data matrix of simulated gene counts, with columns containing next describing which group they came from.
#' @return An ExpressionSet object containing simulated positive and negative samples, with metadata column "classification" indicating which group they belong to.
#' @examples
#' make_sim_eset(test_data)
make_sim_eset <- function(input_data){
  pd <- cbind(as.integer(grepl("Group1", colnames(input_data))), rep("sim", length(colnames(input_data))))
  rownames(pd) <- colnames(input_data)
  colnames(pd) <- c("classification","simulated_status")
  return(ExpressionSet(input_data, phenoData = AnnotatedDataFrame(as.data.frame(pd))))
}


#' Split matrix of counts into a training and a test set
#' 
#' @param input The matrix of counts to be divided into a training and test set
#' @return List containing "training_set" and "test_set" entries, each a matrix of counts. Training set is 80% of the dataset, test set is 20%.
#' @examples
#' split_train_test(input)
split_train_test <- function(input){
  print("inside split_train_test()")
  s <- shuffle(seq(1:dim(input)[2]))
  training_set <- input[,s[1:round(dim(input)[2]*.8)]]
  test_set <- input[,s[(round(dim(input)[2]*.8)+1):dim(input)[2]]]
  return(list("training_set" = training_set, "test_set" = test_set))
}


#' Subset a known dataset from GEO into training and test set.
#' 
#' @param geo_data The geo series object to be split
#' @param pos_regex A regex that can be queried to identify "positive" samples in the GEO series object
#' @param neg_regex A regex that can be queried to identify "negative" samples in the GEO series object
#' @param source_variable The geo field containing the factor of interest
#' @return A list containing {"x": training data, "y": training class ID, "ff": zeros, "xx": test data, "tt": test class ID, "dd": zeros}.
#' @examples
#' subset_known_dataset(geo_data, "BACTERIA", "VIRUS", "characteristics_ch1.2")
subset_known_dataset <- function(geo_data, pos_regex, neg_regex, source_variable){
  print("inside subset_known_dataset()")

  pos_split <- split_train_test(geo_data[, grep(pos_regex, pData(geo_data)[,source_variable])])
  neg_split <- split_train_test( geo_data[, grep(neg_regex, pData(geo_data)[,source_variable])])
  
  print(class(pos_split$training_set))
  print(class(neg_split$training_set))
  print(class(pos_split$test))
  print(class(neg_split$test_set))
  
  full_train <- Biobase::combine(pos_split$training_set, neg_split$training_set)
  full_test <- Biobase::combine(pos_split$test_set, neg_split$test_set)
  
  true_labels_train <- rep(0, length(pData(full_train)[,source_variable]))
  true_labels_train[grep(pos_regex, pData(full_train)[,source_variable])] <- 1
  true_labels_train <- true_labels_train + 1

  true_labels_test <- rep(0, length(pData(full_test)[,source_variable]))
  true_labels_test[grep(pos_regex, pData(full_test)[,source_variable])] <- 1
  true_labels_test <- true_labels_test + 1
  
  print(full_train)
  
  return(list("x" = t(as.matrix(exprs(full_train))), "y" = as.matrix(true_labels_train), "ff" = as.matrix(rep(0, length(true_labels_train))),
              "xx" = t(as.matrix(exprs(full_test))), "tt" = as.matrix(true_labels_test), "dd" = as.matrix(rep(0, length(true_labels_test)))))
  
}


check_params <- function(parameters, required_params){
  if(sum(required_params %in% names(parameters)) == length(required_params)){
    return(TRUE)
  }else{
    logger.info(msg = "WARNING - the following parameters were not found in parameter .json file")
    logger.info(msg = paste("WARNING - PARAMS - ",required_params[!(required_params %in% names(parameters))]))
    return(FALSE)
  }
}

save_data <- function(var, var_name, EXPERIMENT_DIR){
  saveRDS(var, file = paste(c(EXPERIMENT_DIR, var_name,".rds"), collapse=""))
}

get_variable_name <- function(var) {
  return(deparse(substitute(var)))
}

init_log <- function(parameters){
  logger.info(msg = "INIT - Initializing Script")
  logger.info(msg = paste("INIT - PARAMS - ",names(parameters), parameters, sep="\t"))  
}


evaluate_wilcox <- function( list_of_results ){
  
  logger.info(msg="FUNCTION - STATUS - evaluate_wilcox()")
  
  pmat <- list()
  
  for(i in seq(1:length(list_of_results))){
    for(j in seq(1:length(list_of_results))){
      
      if(j > i){
        
        p_values <- list()
        for(k in seq(1:ncol(list_of_results[[1]]))){
          p_values[[k]] <- wilcox.test(list_of_results[[i]][,k], list_of_results[[j]][,k])$p.value
        } 
        pmat[[paste(c("i",i,"_","j",j), collapse="")]] <- p_values
      }
    }
  }
  output_matrix <- matrix(unlist(pmat), nrow=length(pmat[[1]]))
  colnames(output_matrix) <- names(pmat)
  rownames(output_matrix) <- colnames(list_of_results[[1]])
  return(output_matrix)
  
}


plot_relevant_data_withdataset <- function(result_arrays, plot_params, DS_SIZE_RANGE, DIM_RANGE, DATASETS, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette, color_by_pval = TRUE){
  
  logger.info(msg="FUNCTION - STATUS - plot_relevant_data()")
  
  index_1 <- which(DS_SIZE_RANGE == plot_params$DS_SIZE_RANGE)
  index_2 <- which(DIM_RANGE == plot_params$DIM_RANGE)
  index_3 <- which(DATASETS == plot_params$DATASET)
  index_4 <- which(EXP_RANGE_J == plot_params$EXP_RANGE_J)
  index_5 <- which(EXP_RANGE == plot_params$EXP_RANGE)
  
  if(length(index_1) == 0){
    index_1 <- 1:length(DS_SIZE_RANGE)
  }
  if(length(index_2) == 0){
    index_2 <- 1:length(DIM_RANGE)
  }
  if(length(index_3) == 0){
    index_3 <- 1:length(DATASETS)
  }
  if(length(index_4) == 0){
    index_4 <- 1:length(EXP_RANGE_J)
  }
  if(length(index_5) == 0){
    index_5 <- 1:length(EXP_RANGE)
  }
  
  configured_data <- list()
  configured_sd <- list()
  raw_data <- list()
  
  #configure data to plot
  for(i in 1:length(result_arrays)){
    
    # compute the average over the desired part of the matrix
    a = lapply(seq(dim(result_arrays[[i]])[4]), function(x) result_arrays[[i]][index_1, index_2, index_3, x, index_4, index_5])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
    raw_data[[i]] <- do.call(rbind, a)
    print("RAW_DATA")
    print(raw_data)
    plot_data = Reduce("+", a)/length(a)
    print(plot_data)
    configured_data[[i]] <- plot_data
    
    # compute Standard Deviation for each mean value, using Reduce() function and formula: var(x) = E(x^2) - (E(x))^2
    # https://stackoverflow.com/questions/38493741/calculating-standard-deviation-of-variables-in-a-large-list-in-r
    list.squared.mean <- Reduce("+", lapply(a, "^", 2)) / length(a)
    list.variance <- list.squared.mean - plot_data^2
    list.sd <- sqrt(list.variance)
    configured_sd[[i]] <- list.sd
    
  }
  
  
  plot_title <- paste(c("SIZE_RANGE: ", plot_params$DS_SIZE_RANGE, ", DIM_RANGE: ", plot_params$DIM_RANGE, ", DATASET: ", plot_params$DATASET, ", J_RANGE: ", plot_params$EXP_RANGE_J, ", I RANGE: ", plot_params$EXP_RANGE),collapse="")
  logger.info(msg = paste("PLOTTING - PARAMS - ", plot_params$fileroot))
  logger.info(msg = paste("PLOTTING - PARAMS -", plot_title))
  
  #actually generate the plots
  if(class(dim(configured_data[[1]])) == class(NULL)){  # this is a line-plot; plot all lines on top of each other
    
    pdf(paste(plot_params$fileroot,".pdf",sep=""))
    
    print("creating line plot")
    colors <- colorpal 
    
    # GENERATE p-values here
    EW <- evaluate_wilcox(raw_data)
    write.table(EW, file = paste(plot_params$fileroot,"_Ptab.txt",sep=""))
    
    for(i in seq(1:length(configured_data))){
      
      pch_index <- rep(1, nrow(EW))
      if(color_by_pval){
        print("coloring by p-value")
        pch_index <- as.integer(rowSums((EW < .05)[,grepl(i, colnames(EW))]) > 0) + 1
      }
      
      CD = configured_data[[i]]
      CD_sd = configured_sd[[i]]
      xlab_split = strsplit(names(configured_data[[i]]), '_')[[1]]
      xlabel = paste(xlab_split[1:length(xlab_split)-1],collapse=" ")
      real_x_vals <- lapply(strsplit(names(configured_data[[i]]), '_'), function(l){return(as.numeric(l[[3]]))})

      if(i == 1){
        plot(unlist(real_x_vals), as.numeric(CD), xlab = xlabel, ylab = plot_params$performance_metric, ylim = c(0,1), col=colors[i], pch = c(0, 16)[pch_index], main = plot_title)                     # create new plot
        arrows(unlist(real_x_vals), as.numeric(CD) - CD_sd, unlist(real_x_vals), as.numeric(CD) + CD_sd, length=0.05, angle=90, code=3, col="gray28")
      }else{
        points(unlist(real_x_vals), CD, col=colors[i], xlab = xlabel, ylab= plot_params$performance_metric, pch = c(0, 16)[pch_index])                # plot already exists, add new points on top
        arrows(unlist(real_x_vals), as.numeric(CD) - CD_sd, unlist(real_x_vals), as.numeric(CD) + CD_sd, length=0.05, angle=90, code=3, col="gray28")
      }
      
      print(names(result_arrays))
      legend("bottomleft", inset = 0.02, legend = names(result_arrays), col = colors[1:length(result_arrays)], lty=rep(1,length(result_arrays)), cex=.8) #horiz=TRUE, 
      
    }
    dev.off()
    
  }else if(sum(dim(configured_data[[1]]) > 1) == length(dim(configured_data[[1]]))){  # this is a matrix plot; plot multiple matrices next to eachother
    
    overlap_res <- NULL
    
    plot_list = list()
    for(gl in seq(1:length(configured_data))){
      CD = configured_data[[gl]]
      
      # prep plotting data for ggplot 
      if(gl == 1){
        df <- melt(CD)
        df$gl <- rep(names(result_arrays)[gl], nrow(df))
        overlap_res <- df
      }else{
        df <- melt(CD)
        df$gl <- rep(names(result_arrays)[gl], nrow(df))
        overlap_res <- rbind(overlap_res, df)
      }
      
      hmp = pheatmap(CD, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
                     col=colorpal, main = paste(names(result_arrays)[[gl]],plot_params$performance_metric, "\n",plot_title),
                     border_color = "black", fontsize_number = 10,cellheight=25,cellwidth=25, breaks = seq(0,1,by=.001))
      plot_list[[gl]] <- hmp[[4]]
    }

    pdf(paste(plot_params$fileroot, ".pdf",sep=""))
    par(mar=c(5,15,5,5))
    g <- do.call(grid.arrange, plot_list)
    dev.off()
    
    # This can make plots overlapping datasets if you use a command like this:
    #print(df)
    #print(overlap_res)
    ggplot(overlap_res, aes(x=Var2, y=value, colour=Var1, shape = as.factor(gl),
                  group=interaction(Var1, gl))) + geom_point(size = 3) + ylim(0,1) + ggtitle(paste(plot_params$performance_metric, ": ",plot_title)) +
                theme(plot.title = element_text(size = 10))
    ggsave(paste(plot_params$fileroot, "_ggplot.pdf",sep=""), plot = last_plot(), device = NULL, path = NULL,
           scale = 1, width = 8, height = 5, units = c("in"),
           dpi = 300)
    
    
  }
  
}


plot_relevant_data <- function(result_arrays, plot_params, DS_SIZE_RANGE, DIM_RANGE, EXP_RANGE, EXP_RANGE_J, colorpal = my_palette, color_by_pval = TRUE){
  
  logger.info(msg="FUNCTION - STATUS - plot_relevant_data()")
  
  index_1 <- which(DS_SIZE_RANGE == plot_params$DS_SIZE_RANGE)
  index_2 <- which(DIM_RANGE == plot_params$DIM_RANGE)
  index_4 <- which(EXP_RANGE_J == plot_params$EXP_RANGE_J)
  index_5 <- which(EXP_RANGE == plot_params$EXP_RANGE)
  
  if(length(index_1) == 0){
    index_1 <- 1:length(DS_SIZE_RANGE)
  }
  if(length(index_2) == 0){
    index_2 <- 1:length(DIM_RANGE)
  }
  if(length(index_4) == 0){
    index_4 <- 1:length(EXP_RANGE_J)
  }
  if(length(index_5) == 0){
    index_5 <- 1:length(EXP_RANGE)
  }
  
  configured_data <- list()
  configured_sd <- list()
  raw_data <- list()
  
  #configure data to plot
  for(i in 1:length(result_arrays)){
    
    # compute the average over the desired part of the matrix
    a = lapply(seq(dim(result_arrays[[i]])[3]), function(x) result_arrays[[i]][index_1, index_2, x, index_4, index_5])  # get multiple iterations' worth at a contstant error rate (both I and J remain the same), but evaluate all DIM x all DS_SIZEs
    raw_data[[i]] <- do.call(rbind, a)
    plot_data = Reduce("+", a)/length(a)
    print(plot_data)
    configured_data[[i]] <- plot_data
    
    # compute Standard Deviation for each mean value, using Reduce() function and formula: var(x) = E(x^2) - (E(x))^2
    # https://stackoverflow.com/questions/38493741/calculating-standard-deviation-of-variables-in-a-large-list-in-r
    list.squared.mean <- Reduce("+", lapply(a, "^", 2)) / length(a)
    list.variance <- list.squared.mean - plot_data^2
    list.sd <- sqrt(list.variance)
    configured_sd[[i]] <- list.sd
    
  }
  
  
  plot_title <- paste(c("SIZE_RANGE: ", plot_params$DS_SIZE_RANGE, ", DIM_RANGE: ", plot_params$DIM_RANGE, ", J_RANGE: ", plot_params$EXP_RANGE_J, ", I RANGE: ", plot_params$EXP_RANGE),collapse="")
  logger.info(msg = paste("PLOTTING - PARAMS - ", plot_params$fileroot))
  logger.info(msg = paste("PLOTTING - PARAMS -", plot_title))
  
  #actually generate the plots
  if(class(dim(configured_data[[1]])) == class(NULL)){  # this is a line-plot; plot all lines on top of each other
    
    pdf(paste(plot_params$fileroot,".pdf",sep=""))
    
    print("creating line plot")
    colors <- colorpal 
    
    # GENERATE p-values here
    EW <- evaluate_wilcox(raw_data)
    write.table(EW, file = paste(plot_params$fileroot,"_Ptab.txt",sep=""))
    
    for(i in seq(1:length(configured_data))){
      
      pch_index <- rep(1, nrow(EW))
      if(color_by_pval){
        print("coloring by p-value")
        pch_index <- as.integer(rowSums((EW < .05)[,grepl(i, colnames(EW))]) > 0) + 1
      }
      
      CD = configured_data[[i]]
      CD_sd = configured_sd[[i]]
      xlab_split = strsplit(names(configured_data[[i]]), '_')[[1]]
      xlabel = paste(xlab_split[1:length(xlab_split)-1],collapse=" ")
      real_x_vals <- lapply(strsplit(names(configured_data[[i]]), '_'), function(l){return(as.numeric(l[[3]]))})
      
      if(i == 1){
        plot(unlist(real_x_vals), as.numeric(CD), xlab = xlabel, ylab = plot_params$performance_metric, ylim = c(0,1), col=colors[i], pch = c(0, 16)[pch_index], main = plot_title)                     # create new plot
        arrows(unlist(real_x_vals), as.numeric(CD) - CD_sd, unlist(real_x_vals), as.numeric(CD) + CD_sd, length=0.05, angle=90, code=3, col="gray28")
      }else{
        points(unlist(real_x_vals), CD, col=colors[i], xlab = xlabel, ylab= plot_params$performance_metric, pch = c(0, 16)[pch_index])                # plot already exists, add new points on top
        arrows(unlist(real_x_vals), as.numeric(CD) - CD_sd, unlist(real_x_vals), as.numeric(CD) + CD_sd, length=0.05, angle=90, code=3, col="gray28")
      }
      
      print(names(result_arrays))
      legend("bottomleft", inset = 0.02, legend = names(result_arrays), col = colors[1:length(result_arrays)], lty=rep(1,length(result_arrays)), cex=.8) #horiz=TRUE, 
      
    }
    dev.off()
    
  }else if(sum(dim(configured_data[[1]]) > 1) == length(dim(configured_data[[1]]))){  # this is a matrix plot; plot multiple matrices next to eachother
    
    plot_list = list()
    for(gl in seq(1:length(configured_data))){
      CD = configured_data[[gl]]
      hmp = pheatmap(CD, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
                     col=colorpal, main = paste(names(result_arrays)[[gl]],plot_params$performance_metric, "\n",plot_title),
                     border_color = "black", fontsize_number = 10,cellheight=25,cellwidth=25, breaks = seq(0,1,by=.001))
      plot_list[[gl]] <- hmp[[4]]
    }
    pdf(paste(plot_params$fileroot, ".pdf",sep=""))
    par(mar=c(5,15,5,5))
    g <- do.call(grid.arrange, plot_list)
    dev.off()
    
  }
  
}


grab_grob <- function(){
  grid.echo()
  grid.grab()
}


create_new_result_set <- function(DIM_RANGE, DS_SIZE_RANGE, DATASETS, EXP_RANGE_J, EXP_RANGE){
  result_array <- array(
    rep(0, length(DIM_RANGE) * length(DS_SIZE_RANGE) * length(DATASETS) * ITER * length(EXP_RANGE_J) * length(EXP_RANGE)),
    dim = c(length(DIM_RANGE), length(DS_SIZE_RANGE), length(DATASETS), ITER, length(EXP_RANGE_J), length(EXP_RANGE)),
    dimnames = list(
      paste("DIM", DIM_RANGE, sep=("_")),
      paste("DS_SIZE", DS_SIZE_RANGE, sep=("_")),
      paste("DATASETS", DATASETS, sep=("_")),
      paste("ITER", seq(1:ITER), sep=("_")),
      paste("J_FLIP", EXP_RANGE_J, sep=("_")),
      paste("I_FLIP", EXP_RANGE, sep=("_"))
    )
  )
  return(result_array)
}


create_new_result_set_learning_curve <- function(DATASETS, learning_curve_iters){
  logger.info(msg="FUNCTION - STATUS - create_new_result_set()")
  result_array <- array(
    rep(0, length(DATASETS) * ITER * length(learning_curve_iters)),
    dim = c(length(DATASETS), ITER, length(learning_curve_iters)),
    dimnames = list(
      paste("DATASETS", DATASETS, sep=("_")),
      paste("ITER", seq(1:ITER), sep=("_")),
      paste("LC", learning_curve_iters, sep=("_"))
    )
  )
  return(result_array)
}

create_new_feature_set <- function(EXP_RANGE_J, EXP_RANGE, features){
  logger.info(msg="FUNCTION - STATUS - create_new_feature_set()")
  result_array <- array(
    rep(0, ITER * length(EXP_RANGE_J) * length(EXP_RANGE) * length(features)),
    dim = c( ITER, length(EXP_RANGE_J), length(EXP_RANGE), length(features)),
    dimnames = list(
      paste("ITER", seq(1:ITER), sep=("_")),
      paste("J_FLIP", EXP_RANGE_J, sep=("_")),
      paste("I_FLIP", EXP_RANGE, sep=("_")),
      paste("ft", seq(1:length(features)), sep=("_"))
    )
  )
  return(result_array)
}

get_error <- function(weights, test_data, test_labels){
  return(sum(as.numeric((addbias(test_data) %*% weights > 0)) != (as.numeric(castLabel(test_labels,-1)) > 0))/ length(test_labels))
}

castLabel <- function(y, t){
  if (length(y) == 1){
    print('All value of y required to recognise current format')
    return
  }
  if (-1 %in% y){
    # {-1,1} input
    if(t == -1){
      y = y # do nothing, included for clarity
    }else if(t == 0){
      y = (y+1)/2
    }else if(t==2){
      y = (y+3)/2
    }
  }else if(0 %in% y){
    # {0,1} input
    if(t == -1){
      y = y *2 -1
    }else if(t == 0){
      y = y # do nothing, included for clarity
    }else if(t == 2){
      y = y + 1
    }
  }else if (2 %in% y){
    # {1,2} input
    if (t == -1){
      y = y * 2 - 3
    }else if(t == 0){
      y = y -1
    }else if(t == 2){
      y = y # do nothing, included for clarity
    }
  }
  return(y)
}


load_matlab_libraries <- function(){
  logger.info(msg="FUNCTION - STATUS - load_matlab_libraries()")
  #
  # Load all the files required to run the estimation procedures
  #
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/utils/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/kk_utils/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/robust_nlr/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/robust_lr/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/regulariser/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/minFunc');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/netlab');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/libsvm');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/glmnet_matlab');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/minFunc/mex');")
  evaluate(matlab, "addpath('/Users/kkalantar/Documents/Research/NMLHC/robust_logisticregression/third_party_libs/minFunc/compiled');")
}

run_rlr <- function(method, winit, ginit, train_data, train_labels, options, test_data, test_labels){
  
  setVariable(matlab, winit = winit)
  setVariable(matlab, ginit = ginit)
  setVariable(matlab, train_data = train_data)
  setVariable(matlab, train_labels = train_labels)
  setVariable(matlab, common_reg = options$regFunc)
  setVariable(matlab, common_sn = options$sn)
  setVariable(matlab, estG = options$estG)
  setVariable(matlab, test_data = test_data)
  setVariable(matlab, test_labels = test_labels)
  
  # parameters that require True/False argument
  if(options$estG){
    evaluate(matlab, "options.estG = true")
  }else{
    evaluate(matlab, "options.estG = false")
  }
  evaluate(matlab, "options.regFunc = common_reg;")
  evaluate(matlab, "options.verbose = false;")
  evaluate(matlab, "options.sn = common_sn;")
  
  if(method == "rlr"){
    evaluate(matlab, "[w, g, llh] = rlr(winit, ginit, addbias(train_data), train_labels, options);")
  }else if(method == "gammalr"){
    evaluate(matlab, "[w, g, llh] = gammalr(winit, addbias(train_data), train_labels, options);")
  }
  
  w = getVariable(matlab, "w")[[1]]
  g = getVariable(matlab, "g")[[1]]
  llh = getVariable(matlab, "llh")[[1]]

  auc = try(roc(as.numeric(test_labels), as.numeric(addbias(test_data) %*% w))) 
  if(class(auc) == "try-error"){
    auc = list(auc=0)
  }
  
  logger.info(paste("AUC:", auc$auc))
  
  error = get_error(w, test_data, test_labels) 
  
  
  return( list( w = w,
                g = g,
                llh = llh,
                error = error,
                auc = auc) )
  
}




subset_geo <- function(geo_dataset_name, geo_dataset_list, source_variable){
  split <- strsplit(geo_dataset_name, "_")
  print(split)
  pos_regex <- split[[1]][1]
  neg_regex <- split[[1]][2]
  return_value <- subset_known_dataset(geo_dataset_list[[geo_dataset_name]], pos_regex, neg_regex, source_variable)
  return(return_value)
}


generate_data <- function(CLS, DIM, DS_SIZE, N_TEST = 1000, CLS_SEP = 1, DT = 'gen'){
  
  # send variables to matlab server
  setVariable(matlab, CLS = CLS)
  setVariable(matlab, DIM = DIM)
  setVariable(matlab, DS_SIZE = DS_SIZE)
  setVariable(matlab, N_TEST = N_TEST)
  setVariable(matlab, CLS_SEP = CLS_SEP)
  setVariable(matlab, DT = DT)
  
  print("okay")
  
  # use matlab genData function
  evaluate(matlab, "[x, y, ff, xx, tt, dd] = genData(CLS,DIM,DS_SIZE,N_TEST,CLS_SEP,DT);")
  
  x = getVariable(matlab, "x")[[1]]
  colnames(x) <- seq(1:ncol(x))
  
  xx = getVariable(matlab, "xx")[[1]]
  colnames(xx) <- seq(1:ncol(xx))
  
  # get the variables back
  return(list(x = x, 
              y = getVariable(matlab, "y")[[1]], 
              ff = getVariable(matlab, "ff")[[1]], 
              xx = xx, 
              tt = getVariable(matlab, "tt")[[1]], 
              dd = getVariable(matlab, "dd")[[1]]))
}


inject_label_noise <- function(yt, flip_i, flip_j){
  setVariable(matlab, yt = yt)
  setVariable(matlab, flip_i = flip_i)
  setVariable(matlab, flip_j = flip_j)
  
  evaluate(matlab, "[yz, fdz] = injectLabelNoise(yt, [flip_i flip_j]);")
  
  return(list(yz = getVariable(matlab, "yz")[[1]],
              fdz = getVariable(matlab, "fdz")[[1]]))
}


inject_label_noiseR <- function(y, flip_i, flip_j){
  
  fd <- rep(1, length(y)) * -1
  yz <- castLabel(y, -1)
  print("yz")
  print(yz)
  y <- castLabel(y, 2)
  print('y')
  print(y)
  
  flip_rate <- c(flip_i, flip_j)
  
  for(i in c(1,2)){
    print(i)
    prob = rand(length(yz),1)
    idx = intersect(which(y==i), which(prob <= flip_rate[i]))
    print(idx)
    print("yz[idx]")
    print(yz[idx])
    yz[idx] = yz[idx] * -1
    print(yz[idx])
    fd[idx] = fd[idx] * -1
  }
 
  yz = castLabel(yz, 2)
  return(list("yz" = yz, "fd" = fd))
}




standardise <- function(Xt, Xs){
  
  setVariable(matlab, Xt = Xt)
  setVariable(matlab, Xs = Xs)
  
  evaluate(matlab, "[Xt, Xs] = standardise(Xt,Xs);")
  
  return(list(Xt = getVariable(matlab, "Xt")[[1]],
              Xs = getVariable(matlab, "Xs")[[1]]))
}


standardiseR <- function(x, xt){
  DPNT = dim(x[1])
  TPNT = dim(xt[1])
  xa = x;
  
  # standardizing training set uses only training samples
  offset = colMeans(xa)
  var = apply(xa, 2, sd)
  var[var==0] = var[var==0] + 1
  scale = 1/var
  
  X  = x - repmat(offset, DPNT, 1)
  X = X * repmat(scale, DPNT, 1)
  Xt  = xt - repmat(offset, TPNT, 1)
  Xt = xt * repmat(scale, TPNT, 1)
  
  return(list(X,Xt))
}


addbias <- function(x){
  nrow = dim(x)[1]
  x = cbind(matrix(rep(1,nrow), ncol=1), x)
  return(x)  
}



