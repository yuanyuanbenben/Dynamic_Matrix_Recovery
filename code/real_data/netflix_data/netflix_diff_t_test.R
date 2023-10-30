###############################################################################
# train for different time number T
###############################################################################

# environment
library(kernlab)
library(psych)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)
library(npmr)

setwd("/your dictionary/Dynamic_Matrix_Recovery")

source("code/real_data/netflix_data/netflix_DFISTA.R")
source("code/real_data/netflix_data/netflix_baseline_FISTA.R")
source("code/simulation/help_functions.R")
source("code/real_data/help_functions.R")
load("output/diff_t.RData")
args = commandArgs(trailingOnly = TRUE)


p = length(used_id)
q = length(used_movies)
# T = 10,20,50,100,200,500,1000
T_ = as.integer(args[1])
start_time = as.integer(args[2])
print(paste('use time number:',T_))
print(paste('start from:',start_time))
# train data and test data


set.seed(1278451)
if (T_ < 1000){
  train_X1_total_ <- list()
  train_X2_total_ <- list()
  train_Y_total_ <- list()
  test_X1_total_ <- list()
  test_X2_total_ <- list()
  test_Y_total_ <- list()
  for (i in 1:T_) {
    train_X1_total_[[i]] <- train_X1_total[[(i-1)*round(1000/T_)+1]]
    test_X1_total_[[i]] <- test_X1_total[[(i-1)*round(1000/T_)+1]]
    train_X2_total_[[i]] <- train_X2_total[[(i-1)*round(1000/T_)+1]]
    test_X2_total_[[i]] <- test_X2_total[[(i-1)*round(1000/T_)+1]]
    train_Y_total_[[i]] <- train_Y_total[[(i-1)*round(1000/T_)+1]]
    test_Y_total_[[i]] <- test_Y_total[[(i-1)*round(1000/T_)+1]]
    for (j in 2:round(1000/T_)){
      train_X1_total_[[i]] <- c(train_X1_total_[[i]],train_X1_total[[(i-1)*round(1000/T_)+j]])
      test_X1_total_[[i]] <- c(test_X1_total_[[i]],test_X1_total[[(i-1)*round(1000/T_)+j]])
      train_X2_total_[[i]] <- c(train_X2_total_[[i]],train_X2_total[[(i-1)*round(1000/T_)+j]])
      test_X2_total_[[i]] <- c(test_X2_total_[[i]],test_X2_total[[(i-1)*round(1000/T_)+j]])
      train_Y_total_[[i]] <- c(train_Y_total_[[i]],train_Y_total[[(i-1)*round(1000/T_)+j]])
      test_Y_total_[[i]] <- c(test_Y_total_[[i]],test_Y_total[[(i-1)*round(1000/T_)+j]])
    }
  }
}
if (T_ == 1000){
  train_X1_total_ <- train_X1_total
  train_X2_total_ <- train_X2_total
  train_Y_total_ <- train_Y_total
  test_X1_total_ <- test_X1_total
  test_X2_total_ <- test_X2_total
  test_Y_total_ <- test_Y_total
}



batch_size = as.integer(dim(data)[1]/30*0.8)
# selected tuning parameters

h = round(42/100*T_)
lambda = 1.5*100/T_
size = round(1000/T_)
tor=10*100/T_
result_mse = rep(0,1000)
begin <- Sys.time()

for (t in 1:T_) {
  print(paste('current time point:',t))
  datas <- realdata_createdata_func(t,train_X1_total_,train_X2_total_,train_Y_total_,T_,h)
  X1 <- datas[[1]]
  X2 <- datas[[2]]
  Y <- datas[[3]]
  if (t==1){
    M_in <- matrix(3,p,q)
    # len_1 <- length(X1[[t-start_time]])
    # for (i in 1:len_1) {
    #   M_in[X1[[t-start_time]][i],X2[[t-start_time]][i]] <- Y[[t-start_time]][i]
    # }
  }
  M_input <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,itertime=1000,sto=TRUE,batch_size=batch_size,init=FALSE,M_input=M_in,tor=tor)
  # write.csv(M_input,paste("real_data/output/matrix_",t,".csv",sep=""))
  # test error 
  for (index in ((t-1)*size+1):(t*size)){
    result_mse[index] <- test_error_func(M_input,test_X1_total[[index]],test_X2_total[[index]],test_Y_total[[index]])
    print(paste('test mse for t =', index,'is',result_mse[index]))
  }
}
write.csv(result_mse,paste("output/netflix_mse_diff_t_",T_,".csv",sep=""))
end <- Sys.time()
print(difftime(end, begin, units = "sec"))

