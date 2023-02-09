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

# help functions
source("~/Dynamic_Matrix_Recovery/code/simulation/DFISTA.R")
source("~/Dynamic_Matrix_Recovery/code/simulation/baseline_FISTA.R")
source("~/Dynamic_Matrix_Recovery/code/real_data/help_functions.R")
load("~/Dynamic_Matrix_Recovery/code/real_data/netflixdata.RData")

# paraller computing settings
cl.cores = detectCores(logical = F)
cl <- makeCluster(58)
registerDoParallel(cl)


p = length(used_id)
q = length(used_movies)
T_ = 100

#
# DLR method
#
batch_size = as.integer(dim(train_data)[1]/30)
h=42
lambda=11
tor=10
result_mse = rep(0,T_)
for (t in 1:T_) {
#   print(t)
  datas <- realdata_createdata_func(t,X1_total,X2_total,Y_total,T_,h)
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
  M_input <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=30000,sto=TRUE,batch_size=batch_size,init=FALSE,M_input=M_in,tor=tor)
  # test error 
  result_mse[t] <- test_error_func(M_input,test_X1_total[[t]],test_X2_total[[t]],test_Y[[t]])
#   print(result_mse[t])
  write.csv(result_mse,paste("~/Dynamic_Matrix_Recovery/data/mse_",t,".csv",sep=""))
}

#
# Static method
#
batch_size = as.integer(dim(train_data)[1]/30/42)
lambda = 0.55
tor=10
result_mse = rep(0,T_)
for (t in 1:T_) {
  X1 <- X1_total[[t]]
  X2 <- X2_total[[t]]
  Y <- Y_total[[t]]
  M_input = matrix(3,p,q)
  M_input <- base_FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=3000,sto=TRUE,
                        batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)
  # test error 
  result_mse[t] <- test_error_func(M_input,test_X1_total[[t]],test_X2_total[[t]],test_Y[[t]])
#   print(result_mse[t])
  write.csv(result_mse,paste("~/Dynamic_Matrix_Recovery/data/baseline_mse_",t,".csv",sep=""))
}


# cv used
# batch_size = as.integer(dim(train_data)[1]/30/42)
# ret <- foreach(h = c(20,24,28,32,36,42,46,50),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %:% 
#   foreach(lambda=c(5,7,9,11,13,15,17),.packages = c("matlab","psych","kernlab","npmr"),.combine="rbind",.verbose=TRUE) %dopar% {
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
#     write.csv(loss,paste("~/final_version/datas/cv_",h,"_",lambda,".csv",sep = ""))
# }

#
# Tensor method
#
X_tensor <- array(FALSE,dim = c(T_,p,q))
Y_tensor <- array(0,dim = c(T_,p,q))
for (j in 1:T_) {
  index1 <- X1_total[[j]]
  index2 <- X2_total[[j]]
  lens <- length(index1)
  for (i in 1:lens) {
    X_tensor[j,index1[i],index2[i]] <- TRUE
    Y_tensor[j,index1[i],index2[i]] <- Y[[j]][i]
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
    Y_tensor_test[j,index1[i],index2[i]] <- test_Y[[j]][i]
  }
}

Y_tensor_hat <- ADM_TR(X_tensor,Y_tensor,T_,p,q,beta = 0.1,lamda=0.3,c_beta=1,c_lamda=1,itertime = 1000)
Y_tensor_test_hat <- Y_tensor_hat[X_tensor_test]
tensor_mse <- mean_time_mse(Y_tensor_test,X_tensor_test,Y_tensor_test_hat,T_,p,q)
# write.csv(tensor_mse,"~/Dynamic_Matrix_Recovery/output/baseline_mse_tensor.R")

stopCluster(cl) 
