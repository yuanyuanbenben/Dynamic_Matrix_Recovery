###
# the netflix recommandation system example
###

# environment
library(kernlab)
library(psych)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)
library(npmr)

setwd("/your dictionary/Dynamic_Matrix_Recovery")

# input
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  warning("Using DLR method and no link function version as defaut.")
  method = 'DLR'
  link = 'NoLink'
  save_mode = 'NoSave'
  warning('Results will not be saved.')
} else{
  if (length(args) == 1){
    method = args
    link = 'NoLink'
    save_mode = 'NoSave'
    warning("Using no link function version.")
    warning('Results will not be saved.')
  } else {
    if (length(args) == 2){
      method = args[1]
      link = args[2]
      save_mode = 'NoSave'
      warning('Results will not be saved.')
    } else {
      method = args[1]
      link = args[2]
      save_mode = args[3]
    }
  }
}

if (! method %in% c('DLR','Static','TwoStep','Tensor')){
  warning(paste('The method',method,'has not been implemented in the code yet!'))
} else{
  print(paste('Using the method', method))
}

if (! link %in% c('NoLink','Link')){
  simpleError("No correct link version. Please using 'NoLink' for no link function version and 'Link' for link function version.")
} else{
  print(paste('Using',link,'version'))
}

if (! save_mode %in% c('NoSave','Save')){
  warning(paste('The save mode',save_mode,'has not been implemented in the code yet!'))
} else {
  if (save_mode == 'Save'){
    print('Results will be saved.')
  }
}

# help functions

if (link == 'NoLink'){
  source("code/real_data/netflix_data/netflix_DFISTA.R")
} else{
  if (link == 'Link'){
    source("code/real_data/netflix_data/netflix_DFISTA_link.R")
  }
}
source("code/real_data/netflix_data/netflix_baseline_FISTA.R")
source("code/simulation/help_functions.R")
source("code/real_data/help_functions.R")
load("output/netflixdata.RData")

# parallel computing settings
if (method %in% c('Static','TwoStep')){
  cl.cores = detectCores(logical = F)
  cl <- makeCluster(55)
  registerDoParallel(cl)
}
  
  
  
# method = 'DLR'


p = length(used_id)
q = length(used_movies)
T_ = 100

#
# DLR method
#
if (method == 'DLR'){
  batch_size = as.integer(dim(data)[1]/30*0.8)
  # selected tuning parameters
  if (link == 'NoLink'){
     h = 42
     #lambda = 11.1
    lambda = 1.5
  }
  if (link == 'Link'){
    h = 28 
    lambda = 2.1
  }
  tor=10
  result_mse = rep(0,T_)
  begin <- Sys.time()
  for (t in 1:T_) {
    print(paste('current time point:',t))
    datas <- realdata_createdata_func(t,train_X1_total,train_X2_total,train_Y_total,T_,h)
    X1 <- datas[[1]]
    X2 <- datas[[2]]
    Y <- datas[[3]]
    if (t==1){
      M_in <- matrix(3,p,q)
      len_1 <- length(X1[[t]])
      for (i in 1:len_1) {
        M_in[X1[[t]][i],X2[[t]][i]] <- Y[[t]][i]
      }
    }
    M_input <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=1000,sto=TRUE,batch_size=batch_size,init=FALSE,M_input=M_in,tor=tor)
    # test error 
    result_mse[t] <- test_error_func(M_input,test_X1_total[[t]],test_X2_total[[t]],test_Y_total[[t]])
    print(paste('test mse for t =', t,'is',result_mse[t]))
  }
  if (save_mode == 'Save'){
    write.csv(result_mse,paste("real_data/output/netflix_mse_sample_",link,".csv",sep=""))
  }
  end <- Sys.time()
  print(difftime(end, begin, units = "sec"))
}
#
# Static method
#
if (method=='Static'){
  if (link == 'NoLink'){
    lambda =  5.5 # 15/2
  }
  if (link == 'Link'){
    lambda = 21/2
  }
  batch_size = as.integer(dim(data)[1]/30)
  tor=1
  result_mse = rep(0,T_)
  begin <- Sys.time()
  foreach(t = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
    print(paste('current time point:',t))
    X1 <- train_X1_total[[t]]
    X2 <- train_X2_total[[t]]
    Y <- train_Y_total[[t]]
    M_input = matrix(rnorm(p*q),p,q)*0.2+3
    M_inputs <- base_FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=5000,sto=TRUE,
                        batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)
    if (save_mode == 'Save'){
      write.csv(M_inputs,paste("output/baseline_matrix_",t,".csv",sep=""))
    }
    # test error 
    result_mse[t] <- test_error_func(M_input,test_X1_total[[t]],test_X2_total[[t]],test_Y_total[[t]])
    print(paste('test mse for t =', t,'is',result_mse[t]))
  }
  if (save_mode == 'Save'){
    saved_mse = rep(0,T_)
    for (t in 1:T_){
      M_hat <- read.csv(paste("output/baseline_matrix_",t,".csv",sep=""))[,2:(q+1)]
      saved_mse[t] = test_error_func(M_hat,test_X1_total[[t]],test_X2_total[[t]],test_Y_total[[t]])
    }
    write.csv(saved_mse,"output/baseline_mse.csv")
  }
  end <- Sys.time()
  print(difftime(end, begin, units = "sec"))
}

# for two step method results, we first do the Static method and then smooth its results
if (method == 'TwoStep'){
  mse = rep(0,T_)
  h = 3
  mse_ret = foreach(t = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
    array_ = matrix(0,p,q)
    if (t > h & t <= T_ - h){
      for (s in (t-h):(t+h)){
        a <- read.csv(paste("output/baseline_matrix_",s,".csv",sep=""))[,2:(q+1)]
        for (j in 1:p){
          for (k in 1:q){
            array_[j,k] = array_[j,k] + a[[k]][j]
          }
        }
      }
      array_ = array_/(2*h + 1)
    }
    if (t > T_ - h){
      for (s in (t-h):T_){
        a <- read.csv(paste("output/baseline_matrix_",s,".csv",sep=""))[,2:(q+1)]
        for (j in 1:p){
          for (k in 1:q){
            array_[j,k] = array_[j,k] + a[[k]][j]
          }
        }
      }
      array_ = array_/(h + 1 + T_ - t)
    }
    if (t <= h){
      for (s in 1:(t+h)){
        a <- read.csv(paste("output/baseline_matrix_",s,".csv",sep=""))[,2:(q+1)]
        for (j in 1:p){
          for (k in 1:q){
            array_[j,k] = array_[j,k] + a[[k]][j]
          }
        }
      }
      array_ = array_/(h + t)
    }
  mse_t = test_error_func(array_,test_X1_total[[t]],test_X2_total[[t]],test_Y_total[[t]])
  print(mse_t)
  mse_t
  }
  print(mse_ret)
  if (save_mode == 'Save'){
    write.csv(mse_ret,"output/twostep_mse.csv")
  }
}
#
# Tensor method
#
if (method=='Tensor'){
  X_tensor <- array(FALSE,dim = c(T_,p,q))
  Y_tensor <- array(0,dim = c(T_,p,q))
  for (j in 1:T_) {
    index1 <- train_X1_total[[j]]
    index2 <- train_X2_total[[j]]
    lens <- length(index1)
    for (i in 1:lens) {
      X_tensor[j,index1[i],index2[i]] <- TRUE
      Y_tensor[j,index1[i],index2[i]] <- train_Y_total[[j]][i]
    }
  }
  
  X_tensor_test <- array(FALSE,dim = c(T_,p,q))
  Y_tensor_test <- array(0,dim = c(T_,p,q))
  for (j in 1:T_) {
    index1 <- test_X1_total[[j]]
    index2 <- test_X2_total[[j]]
    lens <- length(index1)
    for (i in 1:lens) {
      X_tensor_test[j,index1[i],index2[i]] <- TRUE
      Y_tensor_test[j,index1[i],index2[i]] <- test_Y_total[[j]][i]
    }
  }
  begin <- Sys.time()
  Y_tensor_hat <- ADM_TR(X_tensor,Y_tensor,T_,p,q,beta = 0.1,lamda=0.3,c_beta=1,c_lamda=1,itertime = 50,netflix=TRUE) #0.8,3
  end <- Sys.time()
  print(difftime(end, begin, units = "sec"))
  if (save_mode == 'Save'){
    for (t in 1:T_){
      write.csv(Y_tensor_hat[t,,],paste("output/tensor_",t,".csv",sep = ""))
    }
  }
  mse = rep(0,T_)
  for (t in 1:T_){
    print(t)
    a = read.csv(paste("output/tensor_",t,".csv",sep=""))[,2:(q+1)]
    Y_tensor_test_hat <- a[X_tensor_test[t,,]]
    Y_test = Y_tensor_test[t,,]
    Y_test = Y_test[X_tensor_test[t,,]]
    mse[t] <- sum((Y_tensor_test_hat - Y_test)*(Y_tensor_test_hat - Y_test))/length(test_X1_total[[t]])
    print(mse[t])
  }
  if (save_mode == 'Save'){
    write.csv(mse,"output/baseline_mse_tensor.csv")
  }
}


# cv used
# batch_size = as.integer(dim(data)[1]/30/42*0.8)
# ret <- foreach(h = c(20,24,28,32,36,42,46,50),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %:% 
#   foreach(lambda=c(5,7,9,11,13,15,17),.packages = c("psych","kernlab","npmr"),.combine="rbind",.verbose=TRUE) %dopar% {
#     tor=h/5
#     loss = c()
#     o=1
#     for (t in c(30,70)) {
#     #  for (l in 1:4) {
#       l=1
#       X1 <- list()
#       X2 <- list()
#       Y <- list()
#       if (l==1){
#         for (s in 1:T_) {
#           X1[[s]] <- c(train_X1_total_cv2[[s]],train_X1_total_cv3[[s]],train_X1_total_cv4[[s]])
#           X2[[s]] <- c(train_X2_total_cv2[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
#           Y[[s]] <- c(train_X2_total_cv2[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
#         }
#       }
#       if (l==2){
#         for (s in 1:T_) {
#           X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv3[[s]],train_X1_total_cv4[[s]])
#           X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
#           Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
#         }
#       }
#       if (l==3){
#         for (s in 1:T_) {
#           X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv2[[s]],train_X1_total_cv4[[s]])
#           X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv4[[s]])
#           Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv4[[s]])
#         }
#       }
#       if (l==4){
#         for (s in 1:T_) {
#           X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv2[[s]],train_X1_total_cv3[[s]])
#           X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv3[[s]])
#           Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv3[[s]])
#         }
#       }
#       datas <- realdata_createdata_func(t,X1,X2,Y,T_,h)
#       X1 <- datas[[1]]
#       X2 <- datas[[2]]
#       Y <- datas[[3]]
#       #if (t==1){
#       M_in <- matrix(3,p,q)
#       len_1 <- length(X1[[t]])
#       for (i in 1:len_1) {
#         M_in[X1[[t]][i],X2[[t]][i]] <- Y[[t]][i]
#       }
#       #}
#       M_input <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=30000,sto=TRUE,
#                             batch_size=4675*3*5,init=FALSE,M_input=M_in,tor=tor)
#       # test error 
#       if (l==1) {
#         loss[o] <- test_error_func(M_input,train_X1_total_cv1[[t]],train_X2_total_cv1[[t]],train_Y_total_cv1[[t]])
#       }
#       if (l==2) {
#         loss <- loss + test_error_func(M_input,train_X1_total_cv2[[t]],train_X2_total_cv2[[t]],train_Y_total_cv2[[t]])
#       }
#       if (l==3) {
#         loss <- loss + test_error_func(M_input,train_X1_total_cv3[[t]],train_X2_total_cv3[[t]],train_Y_total_cv3[[t]])
#       }
#       if (l==4) {
#         loss <- loss + test_error_func(M_input,train_X1_total_cv4[[t]],train_X2_total_cv4[[t]],train_Y_total_cv4[[t]])
#       }
#       o = o+1
#      # }
#       loss
#     }
# }
if (method %in% c('Static','TwoStep')){
  stopCluster(cl) 
}
