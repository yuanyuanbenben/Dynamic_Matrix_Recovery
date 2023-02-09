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

source("~/pic_data/cs_FISTA_1.R")
source("~/pic_data/cs_baseline_FISTA.R")
source("~/pic_data/robust_pca.R")
load("~/Dynamic_Matrix_Recovery/code/real_data/vedio_data.RData")

cl.cores = detectCores(logical = F)
cl <- makeCluster(21)
registerDoParallel(cl)

conv_ker = 1/16*matrix(c(1,2,1,2,4,2,1,2,1),3,3)

n = 60000
h = 5

T_ = dim(img_total_b)[1]
p = dim(img_total_b)[2]
q = dim(img_total_b)[3]

set.seed(378461087)
X1_total <- array(sample(2:(p-1),n*T_,replace = TRUE),dim = c(n,T_))
X2_total <- array(sample(2:(q-1),n*T_,replace = TRUE),dim = c(n,T_))



#M_smooth <- function(M1,M2,M3,T_){
 # Y <- array(0,dim = c(2*T_+1,p,q))
  #svd1 <- svd(M1)
  #svd2 <- svd(M2)
  #svd3 <- svd(M3)
  #for (i in (9*T_+1):(10*T_)) {
  #  Y[i-9*T_,,] <- (svd1$u*cos((i-1)*pi/20/T_)+svd2$u*sin((i-1)*pi/20/T_))%*%diag(svd1$d*(10*T_-i+1)/10/T_ + svd2$d*(i-1)/10/T_)%*%t(svd1$v*cos((i-1)*pi/20/T_)+svd2$v*sin((i-1)*pi/20/T_))
  #}
  #Y[T_+1,,] <- M2
  #for (i in 1:T_) {
  #  Y[T_+1+i,,] <-(svd2$u*cos(i*pi/20/T_)+svd3$u*sin(i*pi/20/T_))%*%diag(svd2$d*(10*T_-i)/10/T_ + svd3$d*i/10/T_)%*%t(svd2$v*cos(i*pi/20/T_)+svd3$v*sin(i*pi/20/T_))
  #}
  #return(Y)
#}
#M_total <-  M_smooth(img_total_r[2,,],img_total_r[3,,],img_total_r[4,,],T_)
#T_ = 20
#M_total = M_total[(T_+1-20):(T_+1+20),,]
L_total = array(0,dim = c(3,T_,p,q))
S_total = array(0,dim = c(3,T_,p,q))
foreach (i = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
  #print(i)
  l = rpca_func(img_total_b[i,,],0.25,0.0125)
  write.csv(l[[1]],paste("~/pic_data/outdata/lions_blue_lowrank_",i,".csv",sep = ""))
  write.csv(l[[2]],paste("~/pic_data/outdata/lions_blue_sparse_",i,".csv",sep = ""))
}

for (t in 1:T_) {
  print(t)
  m1 = matrix(0,p,q)
  m1_data = read.csv(paste("~/pic_data/outdata/lions_red_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m1[i,j] = m1_data[i,j]
    }
  }
  L_total[1,t,,] = m1
  m2 = matrix(0,p,q)
  m2_data = read.csv(paste("~/pic_data/outdata/lions_green_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m2[i,j] = m2_data[i,j]
    }
  }
  L_total[2,t,,] = m2
  m3 = matrix(0,p,q)
  m3_data = read.csv(paste("~/pic_data/outdata/lions_blue_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m3[i,j] = m3_data[i,j]
    }
  }
  L_total[3,t,,] = m3
}
for (t in 1:T_) {
  print(t)
  m1 = matrix(0,p,q)
  m1_data = read.csv(paste("~/pic_data/outdata/lions_red_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m1[i,j] = m1_data[i,j]
    }
  }
  S_total[1,t,,] = m1
  m2 = matrix(0,p,q)
  m2_data = read.csv(paste("~/pic_data/outdata/lions_green_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m2[i,j] = m2_data[i,j]
    }
  }
  S_total[2,t,,] = m2
  m3 = matrix(0,p,q)
  m3_data = read.csv(paste("~/pic_data/outdata/lions_blue_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m3[i,j] = m3_data[i,j]
    }
  }
  S_total[3,t,,] = m3
}
# compressing step
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


# data in each epoch t 
cs_realdata_createdata_func <- function(t,X1_total,X2_total,Y_total,T_,h,mode){
  h_ = as.integer(h/2)
  
  if (t - h_ <= 1){
    X1 <- X1_total[,1:(h_+t)]
    X2 <- X2_total[,1:(h_+t)]
    Y <- Y_total[mode,,1:(h_+t)]
  }
  if (t + h_ > T_){
    X1 <- X1_total[,(t-h_):T_]
    X2 <- X2_total[,(t-h_):T_]
    Y <- Y_total[mode,,(t-h_):T_]
  }
  if (t - h_ > 1 & t + h_ <= T_) {
    X1 <- X1_total[,(t-h_):(t+h_)]
    X2 <- X2_total[,(t-h_):(t+h_)]
    Y <- Y_total[mode,,(t-h_):(t+h_)]
  }
  return(list(X1,X2,Y))
}

cs_test_error_func <- function(M,N,p,q){
  return(norm(M-N,type="F")^2/p/q)
}


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
