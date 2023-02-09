###
# some simulations
###

# environment
library(kernlab)
# library(CVXR)
library(psych)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)
library(npmr)

source("~/project_dmc/real_data/large_baseline_FISTA.R")
source("~/project_dmc/real_data/large_FISTA.R")

cl.cores = detectCores(logical = F)
cl <- makeCluster(58)
registerDoParallel(cl)

# real data 1
# netflix recommandation system data

# error functions
test_error_func <- function(M,X1,X2,Y){
  # X1,X2,Y: vector
  len = length(X1)
  err = 0
  for (i in 1:len) {
    err <- err + (M[X1[i],X2[i]] - Y[i])^2
  }
  return(err/len)
}

realdata_createdata_func <- function(t,X1_total,X2_total,Y,T_,h){
  h_ = as.integer(h/2)
  
  if (t - h_ <= 1){
    X1 <- X1_total[1:(h_+t)]
    X2 <- X2_total[1:(h_+t)]
    Y <- Y[1:(h_+t)]
  }
  if (t + h_ > T_){
    X1 <- X1_total[(t-h_):T_]
    X2 <- X2_total[(t-h_):T_]
    Y <- Y[(t-h_):T_]
  }
  if (t - h_ > 1 & t + h_ <= T_) {
    X1 <- X1_total[(t-h_):(t+h_)]
    X2 <- X2_total[(t-h_):(t+h_)]
    Y <- Y[(t-h_):(t+h_)]
  }
  return(list(X1,X2,Y))
}

simulation_func_baseline <- function(U_0,U_1,V_0,V_1,D_0,D_1,T_=100,n_t=100,p=30,q=25,r=5,t=50,eps_sd=0.1,lambda=1.5){
  M <- (cos(pi*t/2/T_)*U_0 + sin(pi*t/2/T_)*U_1)%*%(D_1 + 
                                                      t/T_*D_0)%*%t(cos(pi*t/2/T_)*V_0 + sin(pi*t/2/T_)*V_1)
  X1 <- sample(1:p,n_t,replace = TRUE)
  X2 <- sample(1:q,n_t,replace = TRUE)
  Y <- rep(0,n_t)
  for (i in 1:n_t) {
    Y[i] <- M[X1[i],X2[i]]
  }
  Y <- Y + rnorm(n_t,sd = eps_sd)
  result_baseline <- baseline_FISTA_func(X1,X2,Y,n_t,p,q,lambda=lambda)
  #print(M-result_baseline)
  return(compare_matrix_func(result_baseline,M,p,q)[[1]])
}




p = length(used_id)
q = length(used_movies)
T_ = 100
h = 10
batch_size = as.integer(dim(train_data)[1]/30/42)
lambda = 0.55
tor=10
result_mse = rep(0,T_)
for (t in 1:T_) {
  print(t)
  #datas <- realdata_createdata_func(t,X1_total,X2_total,Y_total,T_,h)
  #X1 <- datas[[1]]
  #X2 <- datas[[2]]
  #Y <- datas[[3]]
  X1 <- X1_total[[t]]
  X2 <- X2_total[[t]]
  Y <- Y_total[[t]]
  #if (t==1){
   # M_input <- matrix(0,p,q)
   # len_1 <- length(X1[[1]])
   # for (i in 1:len_1) {
   #   M_input[X1[[1]][i],X2[[1]][i]] <- Y[[1]][i]
   # }
 # }
  M_input = matrix(3,p,q)
  M_input <- base_FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=3000,sto=TRUE,
                        batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)
  # test error 
  result_mse[t] <- test_error_func(M_input,test_X1_total[[t]],test_X2_total[[t]],test_Y[[t]])
  print(result_mse[t])
}


h=42
lambda=11
batch_size = as.integer(dim(train_data)[1]/30/42)
ret <- foreach(h = c(20,24,28,32,36,42,46,50),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %:% 
  foreach(lambda=c(5,7,9,11,13,15,17),.packages = c("matlab","psych","kernlab","npmr"),.combine="rbind",.verbose=TRUE) %dopar% {
    tor=h/5
    loss = c()
    o=1
    for (t in c(30,70)) {
    #  for (l in 1:4) {
      l=1
      X1 <- list()
      X2 <- list()
      Y <- list()
      if (l==1){
        for (s in 1:T_) {
          X1[[s]] <- c(train_X1_total_cv2[[s]],train_X1_total_cv3[[s]],train_X1_total_cv4[[s]])
          X2[[s]] <- c(train_X2_total_cv2[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
          Y[[s]] <- c(train_X2_total_cv2[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
        }
      }
      if (l==2){
        for (s in 1:T_) {
          X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv3[[s]],train_X1_total_cv4[[s]])
          X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
          Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv3[[s]],train_X2_total_cv4[[s]])
        }
      }
      if (l==3){
        for (s in 1:T_) {
          X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv2[[s]],train_X1_total_cv4[[s]])
          X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv4[[s]])
          Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv4[[s]])
        }
      }
      if (l==4){
        for (s in 1:T_) {
          X1[[s]] <- c(train_X1_total_cv1[[s]],train_X1_total_cv2[[s]],train_X1_total_cv3[[s]])
          X2[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv3[[s]])
          Y[[s]] <- c(train_X2_total_cv1[[s]],train_X2_total_cv2[[s]],train_X2_total_cv3[[s]])
        }
      }
      datas <- realdata_createdata_func(t,X1,X2,Y,T_,h)
      X1 <- datas[[1]]
      X2 <- datas[[2]]
      Y <- datas[[3]]
      #if (t==1){
      M_in <- matrix(3,p,q)
      len_1 <- length(X1[[t]])
      for (i in 1:len_1) {
        M_in[X1[[t]][i],X2[[t]][i]] <- Y[[t]][i]
      }
      #}
      M_input <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=30000,sto=TRUE,
                            batch_size=4675*3*5,init=FALSE,M_input=M_in,tor=tor)
      # test error 
      if (l==1) {
        loss[o] <- test_error_func(M_input,train_X1_total_cv1[[t]],train_X2_total_cv1[[t]],train_Y_total_cv1[[t]])
      }
      if (l==2) {
        loss <- loss + test_error_func(M_input,train_X1_total_cv2[[t]],train_X2_total_cv2[[t]],train_Y_total_cv2[[t]])
      }
      if (l==3) {
        loss <- loss + test_error_func(M_input,train_X1_total_cv3[[t]],train_X2_total_cv3[[t]],train_Y_total_cv3[[t]])
      }
      if (l==4) {
        loss <- loss + test_error_func(M_input,train_X1_total_cv4[[t]],train_X2_total_cv4[[t]],train_Y_total_cv4[[t]])
      }
      o = o+1
     # }
      loss
    }
    #result_mse
     # trick_result <- seq(0,101)
     # for (i in 0:100) {
     #   M_input_ = M_input
     #   trun <- i*0.5/100
     #   M_input_[M_input_ < (1+trun)] = 1
     #   M_input_[(2-trun) < M_input_ & M_input_ < (2+trun)] = 2
     #   M_input_[(3-trun) < M_input_ & M_input_ < (3-trun)] = 3
     #   M_input_[(4-trun) < M_input_ & M_input_ < (4-trun)] = 4
     #   M_input_[(5-trun) < M_input_] = 5
     #   trick_result[i+1] <-  test_error_func(M_input_,test_X1_total[[t]],test_X2_total[[t]],test_Y[[t]])
     # }
     # write.csv(trick_result,paste("~/final_version/datas/mse_",t,".csv",sep = ""))
    write.csv(loss,paste("~/final_version/datas/cv_",h,"_",lambda,".csv",sep = ""))
}
datas <- array(0,c(100,66))
for (i in 1:6) {
  for (j in 1:11) {
    datas[,(i-1)*11+j] <- tuning_para[[i]][j,]
  }
}
datas <- data.frame(datas)
write.csv(datas,"~/final_version/datas/dependent_X_mc.csv")




stopCluster(cl) 
