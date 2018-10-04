
possibly_run_rlr <- possibly(run_rlr, otherwise = NA)

split_train_test <- function(input){
  print("inside split_train_test()")
  s <- shuffle(seq(1:dim(input)[2]))
  training_set <- input[,s[1:round(dim(input)[2]*.8)]]
  test_set <- input[,s[(round(dim(input)[2]*.8)+1):dim(input)[2]]]
  return(list("training_set" = training_set, "test_set" = test_set))
}

subset_known_dataset <- function(geo_data, pos_regex, neg_regex){
  print("inside subset_known_dataset()")
  
  pos_split <- split_train_test(geo_data[, grep(pos_regex, geo_data$source_name_ch1)])
  neg_split <- split_train_test( geo_data[, grep(neg_regex, geo_data$source_name_ch1)])
  
  print(dim(pos_split$training_set))
  print(dim(neg_split$training_set))
  
  print(class(pos_split$training_set))
  
  full_train <- Biobase::combine(pos_split$training_set, neg_split$training_set)
  full_test <- Biobase::combine(pos_split$test_set, neg_split$test_set)
  
  true_labels_train <- rep(0, length(full_train$source_name_ch1))
  true_labels_train[grep(pos_regex, full_train$source_name_ch1)] <- 1
  
  true_labels_test <- rep(0, length(full_test$source_name_ch1))
  true_labels_test[grep(pos_regex, full_test$source_name_ch1)] <- 1
  
  return(list("x" = t(as.matrix(full_train)), "y" = as.matrix(true_labels_train), "ff" = as.matrix(rep(0, length(true_labels_train))),
              "xx" = t(as.matrix(full_test)), "tt" = as.matrix(true_labels_test), "dd" = as.matrix(rep(0, length(true_labels_test)))))
  
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
  logger.info(msg="FUNCTION - STATUS - create_new_result_set()")
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
  
  auc = roc(as.numeric(test_labels), as.numeric(addbias(test_data) %*% w), ci = TRUE)
  error = get_error(w, Xs, ys) 
  
  
  return( list( w = w,
                g = g,
                llh = llh,
                error = error,
                auc = auc) )
  
}


subset_geo <- function(geo_dataset_name, geo_dataset_list){
  split <- strsplit(geo_dataset_name, "_")
  print(split)
  pos_regex <- split[[1]][1]
  neg_regex <- split[[1]][2]
  return_value <- subset_known_dataset(geo_dataset_list[[geo_dataset_name]], pos_regex, neg_regex)
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
  
  # get the variables back
  return(list(x = getVariable(matlab, "x")[[1]], 
              y = getVariable(matlab, "y")[[1]], 
              ff = getVariable(matlab, "ff")[[1]], 
              xx = getVariable(matlab, "xx")[[1]], 
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
