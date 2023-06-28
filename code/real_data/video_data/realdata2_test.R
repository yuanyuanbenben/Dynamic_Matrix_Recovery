# real data 2
# DIVIS 2017 videos

# environment
library(kernlab)
library(psych)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)
library(npmr)


setwd('/your dictionary/Dynamic_Matrix_Recovery')


source("code/real_data/video_data/cs_DFISTA.R")
source("code/real_data/video_data/cs_baseline_FISTA.R")
source("code/real_data/help_functions.R")
source("code/real_data/video_data/robust_pca.R")

load("data/lions_video.RData")

cl.cores = detectCores(logical = F)
cl <- makeCluster(21)
registerDoParallel(cl)

conv_ker = 1/16*matrix(c(1,2,1,2,4,2,1,2,1),3,3)

set.seed(378461087)

# compress step
n = 60000
h = 5
X1_total <- array(sample(2:(p-1),n*T_,replace = TRUE),dim = c(n,T_))
X2_total <- array(sample(2:(q-1),n*T_,replace = TRUE),dim = c(n,T_))
Y_total <- array(0,dim = c(3,n,T_))
for (i in 1:T_) {
  for (j in 1:n) {
    s = X1_total[j,i]
    t = X2_total[j,i]
    Y_total[1,j,i] <- sum(L_total[1,i,c(s-1,s,s+1),c(t-1,t,t+1)]*conv_ker)
    Y_total[2,j,i] <- sum(L_total[2,i,c(s-1,s,s+1),c(t-1,t,t+1)]*conv_ker)
    Y_total[3,j,i] <- sum(L_total[3,i,c(s-1,s,s+1),c(t-1,t,t+1)]*conv_ker)
  }
}

#
# DLR method
#
batch_size = as.integer(n)
lambda = 2.75
tor = 0
label = c("red","green","blue")
result_mse <- foreach(mode=c(1,2,3),.packages = c("psych","kernlab","npmr"),.verbose=TRUE)%:%
  foreach(t = c(5,25,45,65,85),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
    datas <- cs_realdata_createdata_func(t,X1_total,X2_total,Y_total,T_,h,mode)
    X1 <- datas[[1]]
    X2 <- datas[[2]]
    Y <- datas[[3]]
    M_input <- matrix(rnorm(p*q),p,q)
    M_input <- cs_FISTA_func(X1,X2,Y,T_,t,p,q,lambda,conv_ker,itertime=40000,sto=TRUE,
                           batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)
    
  write.csv(construc_pic(t,M_input,mode,S_total),paste("output/lions/lions_",label[mode],"_",t,".csv",sep = ""))
  cs_test_error_func(M_input,L_total[mode,t,,],p,q)
}

# 
# Static method
#
batch_size = as.integer(n/5)
lambda = 0.224
tor=0
label = c("red","green","blue")
result_mse <- foreach(mode=c(1,2,3),.packages = c("psych","kernlab","npmr"),.verbose=TRUE)%:%
  foreach(t = c(5,25,45,65,85),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
    X1 <- X1_total[,t]
    X2 <- X2_total[,t]
    Y <- Y_total[mode,,t]
    M_input <- matrix(rnorm(p*q),p,q)
    M_input <- baseline_cs_FISTA_func(X1,X2,Y,T_,t,p,q,lambda,conv_ker,itertime=40000,sto=TRUE,
                           batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)

  write.csv(construc_pic(t,M_input,mode,S_total),paste("output/lions/baseline_lions_",label[mode],"_",t,".csv",sep = ""))
  cs_test_error_func(M_input,L_total[mode,t,,],p,q)
}

stopCluster(cl)
