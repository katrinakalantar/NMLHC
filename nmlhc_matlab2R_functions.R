

new_results_matrix <- function(EXP_RANGE, EXP_RANGE_J){
  mat <- matrix(data=NA, nrow=length(EXP_RANGE), ncol=length(EXP_RANGE_J))
  rownames(mat) <- EXP_RANGE
  colnames(mat) <- EXP_RANGE_J
  return(mat)
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
