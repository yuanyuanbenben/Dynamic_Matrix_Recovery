#
# The following code is used for do simulation in paper
#

# environment
library(kernlab)
library(psych)
library(ggplot2)
library(reshape2)
# if parallel computing need
library(foreach)
library(doParallel)

# DLR help functions and benchmark Static help functions
source("~/project_dmc/final_version/large_FISTA.R")
source("~/project_dmc/final_version/large_baseline_FISTA.R")
source("~/Dynamic_Matrix_Recovery/code/simulation/help_functions.R")

# paralle computing settings
cl.cores = detectCores(logical = F)
cl <- makeCluster(24)
registerDoParallel(cl)

# simulation 1:
# T = 100, n_t = 100, M(t) is rank 5, 30*25 matrix
# X[j,i,,] is 30*25 matrix in which only one entry is one and other entries are zero
# the structure of M is M(t) = U(t)D(t)V(t)^T,
# where U(t) = cos(pi*t/200)*U_0 +sin(pi*t/200)*U_1, V(t) = cos(pi*t/200)*V_0 +sin(pi*t/200)*V_1, D(t) = D_0 + t/100*D_1
# U_0,U_1 are orthogonal basis and V_0,V_1 are orthogonal basis, D_0 = diag(25,16,9,4,1), D_1 = diag(5,4,3,2,1)
# observations Y[j,i] = Tr(M,X[j,i,,]) + eps, eps ~ N(0,0.1)

# for fixed n_t and all t from 1 to T_=100
T_ = 100
n_t = 30000
p = 500
q = 300
r = 10
h = 40
result_mse <- rep(0,T_)
result_mse_baseline <- rep(0,T_)
X1=0
X2=0
Y=0
M_input = 0
eps_gene = 0
set.seed(1246326361)
N = matrix(rnorm(p*q),p,q)
svd_N <- svd(N)
U_0 <- svd_N$u[,1:r]
U_1 <- svd_N$u[,(r+1):(2*r)]
V_0 <- svd_N$v[,1:r]
V_1 <- svd_N$v[,(r+1):(2*r)]
D_0 <- diag(r:1)
D_1 <- D_0*D_0*10
h=22
tuning_para1 <- foreach(lambda = seq(10,15,0.25),.packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE)%dopar% {
    result_mse = rep(0,6)
    for (t in 1:6) {
      ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,
                                    X2,Y,T_,n_t,p,q,r,t,h,eps_sd=1, M_input = M_input,
                                    lambda = lambda)
      result_mse[t] <- ret[[1]]
      M_input <- ret[[2]]
      X1 <- ret[[3]]
      X2 <- ret[[4]]
      Y <- ret[[5]]
    }
    result_mse
}
print(tuning_para1)
lambda = 13
result_mse = rep(0,T_)
time = rep(0,T_)
for (t in 1:T_) {
  start1 <- Sys.time()
  ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,
                                X2,Y,T_,n_t,p,q,r,t,h,eps_sd=1, M_input = M_input,
                                lambda = lambda)
  result_mse[t] <- ret[[1]]
  M_input <- ret[[2]]
  X1 <- ret[[3]]
  X2 <- ret[[4]]
  Y <- ret[[5]]
  end1 <- Sys.time()
  print(difftime(end1, start1, units = "sec"))
  time[t] <- difftime(end1, start1, units = "sec")
}
datas <- data.frame(tuning_para1)
write.csv(datas,"~/final_version/datas/cv_result.csv")
for (t in 1:T_) {
  print(t)
  ret <- simulation_func_single(U_0,U_1,V_0,V_1,D_0,D_1,X1,X2,Y,T_,n_t,p,q,r,t,h,eps_sd = 1, 
                                M_input = M_input,lambda = lambda,eps_rela = FALSE,X_rela = TRUE,rela_para = 1)
  result_mse[t] <- ret[[1]][[1]]
  M_input = ret[[2]]
  print(result_mse[t])
  X1 = ret[[3]]
  X2 = ret[[4]]
  Y = ret[[5]]
}
lambda=58
t=20
tuning_para_baseline <- foreach(t=1:T_,.packages = c("matlab","psych","kernlab"),.combine="rbind",.verbose=TRUE) %dopar% {
  simulation_func_baseline(U_0,U_1,V_0,V_1,D_0,D_1,T_,n_t,p,q,r,t,eps_sd = 1,lambda = 58)
}
# time test
begin1 <- Sys.time()
simulation_func_baseline(U_0,U_1,V_0,V_1,D_0,D_1,T_,n_t,p,q,r,t,eps_sd = 1,lambda = 58)
end1 <- Sys.time()
print(difftime(end1, begin1, units = "sec"))
#mse_baseline <- rep(0,T_)
#for (i in 1:T_) {
#  mse_baseline[i] <- tuning_para_baseline[i]
#}
write.csv(tuning_para_baseline,"~/final_version/datas/baseline_120000_matrix.csv")
#m <- read.csv("~/final_version/datas/baseline_30000.csv")[,2]
# h=22, lambda=13
#ret = list(0,0,0,0,0)
h=22
lambda=13
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
datas <- data.frame(datas)
write.csv(datas,"~/final_version/datas/dependent_X_mc.csv")
#mse_results <- data.frame("5000"=tuning_para[1,],"10000"=tuning_para[2,],"15000"=tuning_para[3,],"20000"=tuning_para[4,],"25000"=tuning_para[5,],"30000"=tuning_para[6,])
#write.csv(mse_results,"~/final_version/datas/dmc_5000_30000.csv")
#m <- read.csv("~/final_version/datas/dmc_5000_30000.csv")

# phase transition
h=22
tuning_para <- foreach(n_t = seq(5000,30000,1000),
                       lambda = seq(13/6,13,13/30),tor=seq(10,60,2),
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
datas <- data.frame(datas)
write.csv(datas,"~/final_version/datas/phase_transition_precise.csv")


stopCluster(cl) 
