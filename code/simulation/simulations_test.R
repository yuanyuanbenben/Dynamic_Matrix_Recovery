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
source("~/Dynamic-Matrix-Recovery/code/simulation/help_functions.R")

# paralle computing settings
cl.cores = detectCores(logical = F)
cl <- makeCluster(24)
registerDoParallel(cl)

# kernel function used in DLR method
kernel_weight<- function(t,h,T_){
  epan_ker <- function(x,h){
    if (abs(x) > as.integer(h/2)){
      return(0)
    }
    else {
      return(3/4*(1-4*x^2/h^2))
    }
  }
  if (t - as.integer(h/2) <= 1){
    return(sapply((-t+1):as.integer(h/2), epan_ker,h))
  }
  if (t + as.integer(h/2) >= T_){
    return(sapply(as.integer(-h/2):(T_-t), epan_ker,h))
  }
  if (t -  as.integer(h/2) > 1 & t + as.integer(h/2) < T_) {
    return(sapply(as.integer(-h/2):as.integer(h/2), epan_ker,h))
  }
}

# compare matrix difference
compare_matrix_func <- function(M_,M,p,q){
  err_matrix <- M_-M
  mse <- norm(err_matrix,"F")^2/(p*q)
  max_err <- max(abs(err_matrix))
  mean_err <- mean(abs(err_matrix))
  max_err_rate <- max_err/mean(abs(M))
  mean_err_rate <- mean_err/mean(abs(M))
  return(list("mse"=mse,"max_err"=max_err,"mean_err"=mean_err,
              "max_err_rate"=max_err_rate,"mean_err_rate"=mean_err_rate))
}

# compare using cv
compare_cv_func <- function(M_,X1,X2,Y){
  l = length(X1)
  err = 0
  for (i in 1:l) {
    err = err + (M_[X1[i],X2[i]] - Y[i])^2
  }
  return(err/l)
}
# simulation function
simulation_func_single <- function(U_0,U_1,V_0,V_1,D_0,D_1,
                                   X1=FALSE,X2=FALSE,Y=FALSE,T_=100,n_t=100,p=30,q=25,r=5,t=50,
                                   h=9,X_input=FALSE,X_rela=FALSE,rela_para=FALSE,
                                   eps_sd=0.1,
                                   M_input=FALSE,
                                   eps_rela=FALSE,auto_rela=FALSE,eps_gene=FALSE,
                                   lambda=1.5,tor=30){
  # generate simulation data
  # generate matrices M
  
  M <- (cos(pi*t/2/T_)*U_0 + sin(pi*t/2/T_)*U_1)%*%(D_1 + 
                                                      t/T_*D_0)%*%t(cos(pi*t/2/T_)*V_0 + sin(pi*t/2/T_)*V_1)
  
  
  # generate treatment X
  h_ = 2*as.integer(h/2)+1
  if (X_rela){
    repl_num <- as.integer(rela_para*n_t)
    generate_rela_X_func <- function(index1,index2,repl_num){
      nondele_index <- sample(1:n_t,n_t-repl_num,replace = FALSE)
      repl_index1 <- sample(1:p,repl_num,replace = TRUE)
      repl_index2 <- sample(1:q,repl_num,replace = TRUE)
      index1 <- c(index1[nondele_index],repl_index1)
      index2 <- c(index2[nondele_index],repl_index2)
      return(list(index1,index2))
    }
    if (t==1) {
      X1 = array(0,dim = c(as.integer(h/2)+1,n_t))
      X2 = array(0,dim = c(as.integer(h/2)+1,n_t))
      index1 <- sample(1:p,n_t,replace = TRUE)
      index2 <- sample(1:q,n_t,replace = TRUE)
      X1[1,] <- index1
      X2[1,] <- index2
      for (j in 2:(as.integer(h/2)+1)) {
        index <- generate_rela_X_func(X1[j-1,],X2[j-1,],repl_num)
        X1[j,] <- index[[1]]
        X2[j,] <- index[[2]]
      }
    }
    else {
      if (t - as.integer(h/2) <= 1){
        index <- generate_rela_X_func(X1[t+as.integer(h/2)-1,],X2[t+as.integer(h/2)-1,],repl_num)
        X1 <- rbind(X1,index[[1]])
        X2 <- rbind(X2,index[[2]])
      }
      if (t + as.integer(h/2) > T_){
        X1 = X1[2:(as.integer(h/2) + T_ - t + 2),]
        X2 = X2[2:(as.integer(h/2) + T_ - t + 2),]
      }
      if (t - as.integer(h/2) > 1 & t + as.integer(h/2) <= T_) {
        index <- generate_rela_X_func(X1[h_,],X2[h_,],repl_num)
        X1 <- rbind(X1[2:h_,],index[[1]])
        X2 <- rbind(X2[2:h_,],index[[2]])
      }
    }
  }
  else{
    if (t==1) {
      X1 = array(0,dim = c(as.integer(h/2)+1,n_t))
      X2 = array(0,dim = c(as.integer(h/2)+1,n_t))
      for (j in 1:(as.integer(h/2)+1)) {
        index1 <- sample(1:p,n_t,replace = TRUE)
        index2 <- sample(1:q,n_t,replace = TRUE)
        X1[j,] <- index1
        X2[j,] <- index2
      }
    }
    else {
      if (t - as.integer(h/2) <= 1){
        index1 <- sample(1:p,n_t,replace = TRUE)
        index2 <- sample(1:q,n_t,replace = TRUE)
        X1 <- rbind(X1,index1)
        X2 <- rbind(X2,index2)
      }
      if (t + as.integer(h/2) > T_){
        X1 = X1[2:(as.integer(h/2) + T_ - t + 2),]
        X2 = X2[2:(as.integer(h/2) + T_ - t + 2),]
      }
      if (t - as.integer(h/2) > 1 & t + as.integer(h/2) <= T_) {
        index1 <- sample(1:p,n_t,replace = TRUE)
        index2 <- sample(1:q,n_t,replace = TRUE)
        X1 <- rbind(X1[2:h_,],index1)
        X2 <- rbind(X2[2:h_,],index2)
      }
    }
  }
  # generate Y
  if (t == 1) {
    Y = matrix(0,as.integer(h/2)+1,n_t)
    for (j in 1:(as.integer(h/2)+1)) {
      M_ <- (cos(pi*j/2/T_)*U_0 + sin(pi*j/2/T_)*U_1)%*%(D_1 + 
                                                           j/T_*D_0)%*%t(cos(pi*j/2/T_)*V_0 + sin(pi*j/2/T_)*V_1)
      for (i in 1:n_t) {
        Y[j,i] <- M_[X1[j,i],X2[j,i]]
      }
    }
  }
  else {
    if (t - as.integer(h/2) <= 1){
      y = matrix(0,1,n_t)
      M_ <- (cos(pi*(t + as.integer(h/2))/2/T_)*U_0 + sin(pi*(t + as.integer(h/2))/2/T_)*U_1)%*%(D_1 + 
                                                                                                   (t + as.integer(h/2))/T_*D_0)%*%t(cos(pi*(t + as.integer(h/2))/2/T_)*V_0 + sin(pi*(t + as.integer(h/2))/2/T_)*V_1)
      for (i in 1:n_t) {
        y[1,i] <- M_[X1[(t + as.integer(h/2)),i],X2[(t + as.integer(h/2)),i]]
      }
      Y <- rbind(Y,y)
    }
    if (t + as.integer(h/2) > T_){
      Y = Y[2:(as.integer(h/2) + T_ - t + 2),]
    }
    if (t - as.integer(h/2) > 1 & t + as.integer(h/2) <= T_){
      Y[1:(h_-1),] = Y[2:h_,]
      M_ <- (cos(pi*(t + as.integer(h/2))/2/T_)*U_0 + sin(pi*(t + as.integer(h/2))/2/T_)*U_1)%*%(D_1 + 
                                                                                                   (t + as.integer(h/2))/T_*D_0)%*%t(cos(pi*(t + as.integer(h/2))/2/T_)*V_0 + sin(pi*(t + as.integer(h/2))/2/T_)*V_1)
      for (i in 1:n_t) {
        Y[h_,i] <- M_[X1[h_,i],X2[h_,i]]
      }
    }
  }
  
  if (eps_rela){
    resid_para <- sqrt(1-auto_rela**2)
    if (t==1){
      eps = matrix(0,1+as.integer(h/2),n_t)
      eps_gene = matrix(rnorm(p*q,sd=eps_sd),p,q)
      for (j in 1:(1+as.integer(h/2))) {
        for (i in 1:n_t) {
          eps[j,i] <- eps_gene[X1[j,i],X2[j,i]]
        }
        #update eps_gene
        eps_gene <- auto_rela*eps_gene + resid_para*matrix(rnorm(p*q,sd=eps_sd),p,q)
      }
      Y <- Y + eps
    }
    if (t + as.integer(h/2) <= T_){
      for (i in 1:n_t) {
        Y[min(t+as.integer(h/2),h_),i] <- Y[min(t+as.integer(h/2),h_),i] +
          eps_gene[X1[min(t + as.integer(h/2),h_),i],X2[min(t + as.integer(h/2),h_),i]]
      }
      eps_gene <- auto_rela*eps_gene + resid_para*matrix(rnorm(p*q,sd=eps_sd),p,q)
    }
  }
  else{
    if (t==1){
      eps = matrix(rnorm((1+as.integer(h/2))*n_t),1+as.integer(h/2),n_t)
      Y <- Y + eps
    }
    if (t + as.integer(h/2) <= T_){
      Y[min(t+as.integer(h/2),h_),] <- Y[min(t+as.integer(h/2),h_),] + rnorm(n_t,sd=eps_sd)
    }
  }
  if (t==1){
    M_init <- matrix(0,p,q)
    for (i in 1:n_t) {
      M_init[X1[1,i],X2[1,i]] <- Y[1,i]
    }
    ret_FISTA <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,sto = TRUE,batch_size = as.integer(n_t*h_/10),init = FALSE,M_input=M_init,tor=tor)
  }
  else{
    ret_FISTA <- FISTA_func(X1,X2,Y,T_,t,h,p,q,lambda,sto = TRUE,batch_size = as.integer(n_t*h_/10),init=FALSE,M_input=M_input,tor=tor)
  }
  #svd_ret_FISTA <- svd(ret_FISTA)
  #M_FISTA <- svd_ret_FISTA$u[,1:r]%*%diag(svd_ret_FISTA$d[1:r])%*%t(svd_ret_FISTA$v[,1:r])
  
  # cv_x
  # x1 <- sample(1:p,round(n_t/5),replace = TRUE)
  # x2 <- sample(1:q,round(n_t/5),replace = TRUE)
  # y <- rep(0,round(n_t/5))
  # for (i in 1:round(n_t/5)){
  #  y[i] = M[x1[i],x2[i]] + rnorm(1,sd=eps_sd)
  #}
  #result_FISTA <- compare_cv_func(ret_FISTA,x1,x2,y)
  result_FISTA <- compare_matrix_func(ret_FISTA,M,p,q)
  if (eps_rela) {
    return(list(result_FISTA,ret_FISTA,X1,X2,Y,eps_gene))
  }
  return(list(result_FISTA,ret_FISTA,X1,X2,Y))
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
  M_init <- matrix(0,p,q)
  #for (i in 1:n_t) {
  #  M_init[X1[i],X2[i]] <- Y[i]
  #}
  result_baseline <- baseline_FISTA_func(X1,X2,Y,n_t,p,q,lambda=lambda,init = TRUE,M_input = M_init)
  #print(M-result_baseline)
  #return(result_baseline)
  return(compare_matrix_func(result_baseline,M,p,q))
}
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
