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

source("~/Dynamic_Matrix_Recovery/code/real_data/cs_DFISTA.R")
source("~/Dynamic_Matrix_Recovery/code/real_data/cs_baseline_FISTA.R")
source("~/Dynamic_Matrix_Recovery/code/real_data/help_functions.R")
source("~/Dynamic_Matrix_Recovery/code/real_data/robust_pca.R")
load("~/Dynamic_Matrix_Recovery/code/real_data/vedio_data.RData")

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
batch_size = as.integer(n/5)
lambda = 0.224#2.75
tor=0
label = c("red","green","blue")
#result_mse = rep(0,T_)
#result_pic = array(0,dim = c(T_,p,q))
result_mse <- foreach(mode=c(1,2,3),.packages = c("psych","kernlab","npmr"),.verbose=TRUE)%:%
  foreach(t = c(5,25,45,65,85),.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
 # datas <- cs_realdata_createdata_func(t,X1_total,X2_total,Y_total,T_,h,mode)
 # X1 <- datas[[1]]
 # X2 <- datas[[2]]
 # Y <- datas[[3]]
    X1 <- X1_total[,t]
    X2 <- X2_total[,t]
    Y <- Y_total[mode,,t]
  M_input <- matrix(rnorm(p*q),p,q)
  M_input <- baseline_cs_FISTA_func(X1,X2,Y,T_,t,p,q,lambda,conv_ker,itertime=40000,sto=TRUE,
                           batch_size=batch_size,init=FALSE,M_input=M_input,tor=tor)
  # test error 
  #result_mse[t] <- 
  #result_pic[t,,] <- M_input
  #print(result_mse[t])
  write.csv(construc_pic(t,M_input,mode),paste("~/pic_data/outdata/baseline_lions_",label[mode],"_",t,".csv",sep = ""))
  cs_test_error_func(M_input,L_total[mode,t,,],p,q)
}

baseline_red = array(0,dim=c(5,480,854))
baseline_green = array(0,dim=c(5,480,854))
baseline_blue = array(0,dim=c(5,480,854))

for (i in 1:5) {
  t = i*20 - 15
  ret = construc_func(t)
  baseline_red[i,,] = ret[[1]]
  baseline_green[i,,] = ret[[2]]
  baseline_blue[i,,] = ret[[3]]
}

baseline_sparse_red = array(0,dim=c(5,480,854))
baseline_sparse_green = array(0,dim=c(5,480,854))
baseline_sparse_blue = array(0,dim=c(5,480,854))

for (i in 1:5) {
  t = i*20 - 15
  ret = construc_func(t)
  baseline_sparse_red[i,,] = ret[[1]]
  baseline_sparse_green[i,,] = ret[[2]]
  baseline_sparse_blue[i,,] = ret[[3]]
}
baseline_low_red = baseline_red - baseline_sparse_red
baseline_low_green = baseline_green - baseline_sparse_green 
baseline_low_blue = baseline_blue - baseline_sparse_blue

twostep_low_red = array(0,dim=c(5,480,854))
twostep_low_green = array(0,dim=c(5,480,854))
twostep_low_blue = array(0,dim=c(5,480,854))

twostep_low_red[1,,] = (5*baseline_low_red[1,,] + baseline_low_red[2,,])/6
twostep_low_red[2,,] = (baseline_low_red[1,,] + 5*baseline_low_red[2,,] + baseline_low_red[3,,])/7
twostep_low_red[3,,] = (baseline_low_red[2,,] + 5*baseline_low_red[3,,] + baseline_low_red[4,,])/7
twostep_low_red[4,,] = (baseline_low_red[3,,] + 5*baseline_low_red[4,,] + baseline_low_red[5,,])/7
twostep_low_red[5,,] = (5*baseline_low_red[5,,] + baseline_low_red[4,,])/6

twostep_low_green[1,,] = (5*baseline_low_green[1,,] + baseline_low_green[2,,])/6
twostep_low_green[2,,] = (baseline_low_green[1,,] + 5*baseline_low_green[2,,] + baseline_low_green[3,,])/7
twostep_low_green[3,,] = (baseline_low_green[2,,] + 5*baseline_low_green[3,,] + baseline_low_green[4,,])/7
twostep_low_green[4,,] = (baseline_low_green[3,,] + 5*baseline_low_green[4,,] + baseline_low_green[5,,])/7
twostep_low_green[5,,] = (5*baseline_low_green[5,,] + baseline_low_green[4,,])/6

twostep_low_blue[1,,] = (5*baseline_low_blue[1,,] + baseline_low_blue[2,,])/6
twostep_low_blue[2,,] = (baseline_low_blue[1,,] + 5*baseline_low_blue[2,,] + baseline_low_blue[3,,])/7
twostep_low_blue[3,,] = (baseline_low_blue[2,,] + 5*baseline_low_blue[3,,] + baseline_low_blue[4,,])/7
twostep_low_blue[4,,] = (baseline_low_blue[3,,] + 5*baseline_low_blue[4,,] + baseline_low_blue[5,,])/7
twostep_low_blue[5,,] = (5*baseline_low_blue[5,,] + baseline_low_blue[4,,])/6

twostep_red = array(0,dim=c(5,480,854))
twostep_green = array(0,dim=c(5,480,854))
twostep_blue = array(0,dim=c(5,480,854))

for (i in 1:5) {
  twostep_red[i,,] = twostep_low_red[i,,] + baseline_sparse_red[i,,]
  twostep_green[i,,] = twostep_low_green[i,,] + baseline_sparse_green[i,,]
  twostep_blue[i,,] = twostep_low_blue[i,,] + baseline_sparse_blue[i,,]
}

twostep_red[twostep_red>1]=1
twostep_green[twostep_green>1]=1
twostep_blue[twostep_blue>1]=1
twostep_red[twostep_red<0]=0
twostep_green[twostep_green<0]=0
twostep_blue[twostep_blue<0]=0
baseline_low_red = array(0,dim=c(5,480,854))
baseline_low_green = array(0,dim=c(5,480,854))
baseline_low_blue = array(0,dim=c(5,480,854))
for (i in 1:5) {
  t = i*20 - 15
  ret = construc_func(t)
  baseline_low_red[i,,] = ret[[1]]
  baseline_low_green[i,,] = ret[[2]]
  baseline_low_blue[i,,] = ret[[3]]
}
mse = 0
for (i in 1:5) {
  mse = mse + cs_test_error_func(twostep_low_red[i,,],baseline_low_red[i,,],p,q) + 
    cs_test_error_func(twostep_low_green[i,,],baseline_low_green[i,,],p,q)+
    cs_test_error_func(twostep_low_blue[i,,],baseline_low_blue[i,,],p,q)
}

construc_pic <- function(t,l,mode){
#  l = result_pic[t,,]
  s = S_total[mode,t,,]
#  s[abs(s)<0.06] = 0
  m = s+l
  m[m<0]=0
  m[m>1]=1
  return(m)
}
construc_func <- function(t){
  m_data_red = read.csv(paste("~/pic_data/outdata/lions_red_lowrank_",t,".csv",sep = ""))[,2:855]
  m_data_green = read.csv(paste("~/pic_data/outdata/lions_green_lowrank_",t,".csv",sep = ""))[,2:855]
  m_data_blue = read.csv(paste("~/pic_data/outdata/lions_blue_lowrank_",t,".csv",sep = ""))[,2:855]
  m_red = matrix(0,p,q)
  m_green = matrix(0,p,q)
  m_blue = matrix(0,p,q)
  for (i in 1:p) {
    for (j in 1:q) {
      m_red[i,j] = m_data_red[i,j]
      m_green[i,j] = m_data_green[i,j]
      m_blue[i,j] = m_data_blue[i,j]
    }
  }
  return(list(m_red,m_green,m_blue))
}
a =0
for (i in 1:3) {
  for (j in 1:5) {
    a = a +norm(S_total[i,5+(j-1)*20,,]+S_total[i,5+(j-1)*20,,],"F")
  }
}
b =0
for (i in 1:3) {
  for (j in 1:5) {
    b = b + result_mse[[i]][[j]]
  }
}

stopCluster(cl)
