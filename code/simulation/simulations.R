#
# The following code is used for simulations in paper "Dynamic Matrix Recovery"
#

# environment
library(kernlab)
library(psych)
library(ggplot2)
library(reshape2)
# if parallel computing need
library(foreach)
library(doParallel)

# help functions
source("~/Dynamic_Matrix_Recovery/code/simulation/DFISTA.R")
source("~/Dynamic_Matrix_Recovery/code/simulation/baseline_FISTA.R")
source("~/Dynamic_Matrix_Recovery/code/simulation/help_functions.R")

# paraller computing settings
cl.cores = detectCores(logical = F)
cl <- makeCluster(24)
registerDoParallel(cl)

#
# simulation 1:
# independent case
# T = 100, n_t = 100, M(t) is rank 10, 500*300 matrix
# the structure of M(t) is M(t) = U(t)D(t)V(t)^T,
# where U(t) = cos(pi*t/2)*U_0 +sin(pi*t/2)*U_1, V(t) = cos(pi*t/2)*V_0 +sin(pi*t/2)*V_1, D(t) = D_0 + t*D_1
# U_0,U_1 are orthogonal basis and V_0,V_1 are orthogonal basis, D_0 = 10*diag(100,81,...,1), D_1 = diag(10,9,...,1)
# observations Y = Tr(M,X) + eps, eps ~ N(0,1)
#

# for fixed n_t and all t from 1 to T_=100
T_ = 100
n_t = 30000 # 120000 for benchmark Static and TwoStep 
p = 500
q = 300
r = 10
# no use
X1 = 0
X2 = 0
Y = 0
M_input = 0
eps_gene = 0
# fixed seed
set.seed(1246326361)
# construct U(t), D(t) and V(t)
N = matrix(rnorm(p*q),p,q)
svd_N <- svd(N)
U_0 <- svd_N$u[,1:r]
U_1 <- svd_N$u[,(r+1):(2*r)]
V_0 <- svd_N$v[,1:r]
V_1 <- svd_N$v[,(r+1):(2*r)]
D_0 <- diag(r:1)
D_1 <- D_0*D_0*10

#
# DLR:our proposed algorithm 1
#
# calculate bandwidth h and tuning parameter lambda
h=22
lambda = 13
# record mse results and conputational time 
result_mse_DLR1 = rep(0,T_)
time_DLR1 = rep(0,T_)
for (t in 1:T_) {
  start <- Sys.time()
  ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,
                                X2,Y,T_,n_t,p,q,r,t,h,eps_sd=1, M_input = M_input,
                                lambda = lambda)
  result_mse_DLR1[t] <- ret[[1]]
  M_input <- ret[[2]]
  X1 <- ret[[3]]
  X2 <- ret[[4]]
  Y <- ret[[5]]
  end <- Sys.time()
  #print(difftime(end, start, units = "sec"))
  time_DLR1[t] <- difftime(end1, start1, units = "sec")
}

# find the influences of rho and tau on the estimation results
# change settings for n_t from 5000 to 30000 and T from 100 to 300 
# combine those results mse as combine_mse and save it 
# write.csv(combine_mse,"~/Dynamic_Matrix_Recovery/output/dmc_5000_30000.csv")
tuning_para <- foreach(n_t = seq(5000,30000,1000),
                       .packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE) %dopar% {
  result_mse = rep(0,T_)
  for (t in 1:T_) {
    ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,
                                  X2,Y,T_,n_t,p,q,r,t,h,eps_sd = 1, M_input = M_input,lambda = lambda,tor=tor)
    result_mse[t] <- ret[[1]][[1]]
    M_input <- ret[[2]]
    X1 <- ret[[3]]
    X2 <- ret[[4]]
    Y <- ret[[5]]
  }
  result_mse
  }
datas <- array(0,c(100,26))
for (i in 1:100) {
  for (j in 1:26) {
    datas[i,j] <- tuning_para[j,i]
  }
}
# save
# datas <- data.frame(datas)
# write.csv(datas,"~/Dynamic_Matrix_Recovery/output/phase_transition_precise.csv")

#
# DLR: random initialization
#
h=22
lambda = 13
result_mse_DLR2 <- rep(0, T_)
time_DLR2 = rep(0,T_)
# paraller version
foreach(t=1:T_,.packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE)%dopar% {
  start <- Sys.time()
  ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,X2,Y,T_,
                                n_t,p,q,r,t,h,eps_sd=1, M_input = M_input,lambda = lambda)
  result_mse_DLR2[t] <- ret[[1]]
  end <- Sys.time()
  time_DLR2[t] <- difftime(end1, start1, units = "sec")
}

#
# Static
#
# calculate bandwidth h and tuning parameter lambda
lambda=58
t=20
result_mse_Static <- rep(0, T_)
time_Static = rep(0,T_)
M_estimate_Static <- array(0,dim=c(T_,p,q)
foreach(t=1:T_,.packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE)%dopar% {
  begin <- Sys.time()
  ret <- simulation_func_baseline(U_0,U_1,V_0,V_1,D_0,D_1,T_,n_t,p,q,r,t,eps_sd = 1,lambda = 58)
  result_mse_Static[t] <- ret[[1]]
  M_estimate_Static[t,,] <- ret[[2]]
  end <- Sys.time()
  difftime(end1, begin1, units = "sec")
}
#save                        
#write.csv(result_mse_Static,"~/Dynamic_Matrix_Recovery/output/baseline_120000.csv")
#write.csv(M_estimate_Static,"~/Dynamic_Matrix_Recovery/output/baseline_120000_matrix.csv")

#
# TwoStep
#                  
matrix_data <- read.csv("~/Dynamic_MAtrix_Recovery/output/baseline_120000_matrix.csv")[,2:301]
matrix_array <- array(0,dim = c(T_,p,q))
for (i in 1:T_) {
  print(i)
  for (j in 1:p) {
    for (k in 1:q) {
      matrix_array[i,j,k] <- matrix_data[(i-1)*p+j,k]
    }
  }
}
result_mse <- rep(0,T_)
h=24
lambda = 1
for (t in 1:T_) {
  M_in <- createdata_func(matrix_array,t,T_,h)
  M <- local_smooth(M_in,T_,t,h,p,q,lambda = lambda)
  M_ =  (cos(pi*t/2/T_)*U_0 + sin(pi*t/2/T_)*U_1)%*%(D_1 + 
                                                       t/T_*D_0)%*%t(cos(pi*t/2/T_)*V_0 + sin(pi*t/2/T_)*V_1)
  result_mse[t] <- compare_matrix_func(M,M_,p,q)[[1]]
  print(result_mse[t])
}
#write.csv(result_mse,"~/Dynamic_MAtrix_Recovery/output/localsmooth_120000_2.csv")                                                   
                           
#
# Tensor
#
M <- array(0,dim = c(T_,p,q))
for (j in 1:T_) {
  M[j,,]<- (cos(pi*j/2/T_)*U_0 + sin(pi*j/2/T_)*U_1)%*%(D_1 + 
                                                          j/T_*D_0)%*%t(cos(pi*j/2/T_)*V_0 + sin(pi*j/2/T_)*V_1)
}

X <- array(FALSE,dim = c(T_,p,q))
Y <- array(0,dim = c(T_,p,q))
for (j in 1:T_) {
  index1 <- sample(1:p,n_t,replace = TRUE)
  index2 <- sample(1:q,n_t,replace = TRUE)
  for (i in 1:n_t) {
    X[j,index1[i],index2[i]] <- TRUE
    Y[j,index1[i],index2[i]] <- M[j,index1[i],index2[i]] + rnorm(1,sd=1)
  }
}

M_hat <- ADM_TR(X,Y,T_,p,q,beta = 0.1,lamda=0.3,c_beta=1,c_lamda=1,itertime = 3000)
# a = c(1/3,1/3,1/3)
# rho=10
# M_hat <- HaLRTC(X,Y,a,T_,p,q,rho=rho,itertime = 10000)
# write.csv(result_mse,"~/Dynamic_MAtrix_Recovery/output/tensor_30000.csv")      

                           
#
# simulation2
# dependent case
# other settings are the same as in simulation1
#
tuning_para <- foreach(eps_sd=c(0.1,0.5,1,1.5,2,5),.packages = c("matlab","psych","kernlab"),.verbose=TRUE) %:% 
foreach(rela_para=seq(0,1,0.1),.packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE) %dopar% {
  result_mse = rep(0,T_)
  for (t in 1:T_) {
  ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,
                                  X2,Y,T_,n_t,p,q,r,t,h,eps_sd=eps_sd, M_input = M_input,
                                  lambda = lambda,X_rela = TRUE,rela_para = rela_para)
  result_mse[t] <- ret[[1]][[1]]
  M_input <- ret[[2]]
  X1 <- ret[[3]]
  X2 <- ret[[4]]
  Y <- ret[[5]]
  }
  result_mse
}                
datas <- array(0,c(100,66))
for (i in 1:6) {
  for (j in 1:11) {
    datas[,(i-1)*11+j] <- tuning_para[[i]][j,]
  }
}
#datas <- data.frame(datas)
#write.csv(datas,"~/Dynamic_Matrix_Recovery/output/dependent_X_mc.csv")

stopCluster(cl) 
