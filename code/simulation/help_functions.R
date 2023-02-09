# some help functions

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

# compare matrix difference for matrix M and M_ with dimension p*q
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

# compare using cv for estimated M with sample X1,X2,Y
compare_cv_func <- function(M_,X1,X2,Y){
  l = length(X1)
  err = 0
  for (i in 1:l) {
    err = err + (M_[X1[i],X2[i]] - Y[i])^2
  }
  return(err/l)
}

#' DLR simulation function at a single time point t
#' @param U_0,U_1,V_0,V_1,D_0,D_1: matrix->construct M(t)
#' @param X1,X2,Y: array->samples with Y = M[X1,X2] + xi or FALSE-> random samples
#' @param T_,n_t,t: int->total time points,  sample size in each time, current time
#' @param p,q,r: int->dimensions and rank of underground truth M(t)
#' @param h,lambda: num->bandwidth and tuning parameter
#' @param M_input: matrix->initial matrix or FALSE->random initialization 
#' @param X_rela,eps_rela: bool-> whether dependence exist for X and xi
#' @param rela_para,auto_rela: num-> control dependent coefficients
#' @param eps_sd:num->standard variance of xi
#' @param eps_gene: matrix->dependent xi construction of FALSE-> random sturcture
#' @param tor: num-> stop condition
#' @return list->result_FISTA:MSE,ret_FISTA:estimated M(t):X1,X2,Y:sample used currently
simulation_func_single <- function(U_0,U_1,V_0,V_1,D_0,D_1,
                                   X1=FALSE,X2=FALSE,Y=FALSE,T_=100,n_t=100,p=30,q=25,r=5,t=50,
                                   h=9,X_rela=FALSE,rela_para=FALSE,
                                   eps_sd=0.1,
                                   M_input=FALSE,
                                   eps_rela=FALSE,auto_rela=FALSE,eps_gene=FALSE,
                                   lambda=1.5,tor=30){
  
  # generate simulation data
  
  # generate matrices M(t)
  M <- (cos(pi*t/2/T_)*U_0 + sin(pi*t/2/T_)*U_1)%*%(D_1 + 
                                                      t/T_*D_0)%*%t(cos(pi*t/2/T_)*V_0 + sin(pi*t/2/T_)*V_1)
  
  
  # generate treatment X
  h_ = 2*as.integer(h/2)+1
  # dependence of X exist
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
    # start time t=1
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
      # time points that near begining
      if (t - as.integer(h/2) <= 1){
        index <- generate_rela_X_func(X1[t+as.integer(h/2)-1,],X2[t+as.integer(h/2)-1,],repl_num)
        X1 <- rbind(X1,index[[1]])
        X2 <- rbind(X2,index[[2]])
      }
      # time points that near ending
      if (t + as.integer(h/2) > T_){
        X1 = X1[2:(as.integer(h/2) + T_ - t + 2),]
        X2 = X2[2:(as.integer(h/2) + T_ - t + 2),]
      }
      # regular time points
      if (t - as.integer(h/2) > 1 & t + as.integer(h/2) <= T_) {
        index <- generate_rela_X_func(X1[h_,],X2[h_,],repl_num)
        X1 <- rbind(X1[2:h_,],index[[1]])
        X2 <- rbind(X2[2:h_,],index[[2]])
      }
    }
  }
  # independent case for X 
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
  
  # generate Y = M[X1,X2]
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
  
  # add noise xi for Y
  # dependent case for xi
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
  #independent case for xi
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

#' simulation functon for Static
#' @param U_0,U_1,V_0,V_1,D_0,D_1: matrix->construct M(t)
#' @param X1,X2,Y: array->samples with Y = M[X1,X2] + xi or FALSE-> random samples
#' @param T_,n_t,t: int->total time points,  sample size in each time, current time
#' @param p,q,r: int->dimensions and rank of underground truth M(t)
#' @param lambda: num-> tuning parameter
#' @param eps_sd:num->standard variance of xi
#' @return list-> MSE, estimated M(t)
simulation_func_baseline <- function(U_0,U_1,V_0,V_1,D_0,D_1,X1=FALSE,X2=FALSE,Y=FALSE,T_=100,n_t=100,p=30,q=25,r=5,t=50,eps_sd=0.1,lambda=1.5){
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
  result_baseline <- baseline_FISTA_func(X1,X2,Y,n_t,p,q,lambda=lambda,init = TRUE,M_input = M_init)
  return(list(compare_matrix_func(result_baseline,M,p,q),result_baseline))
}

# from results of Static method to TwoStep method
# using an additional local smooth

createdata_func <- function(matrix_array,t,T_,h){
  h_ = as.integer(h/2)
  
  if (t - h_ <= 1){
    m <- matrix_array[1:(h_+t),,]
  }
  if (t + h_ > T_){
    m <- matrix_array[(t-h_):T_,,]
  }
  if (t - h_ > 1 & t + h_ <= T_) {
    m <- matrix_array[(t-h_):(t+h_),,]
  }
  return(m)
}

threshold_func <- function(x) sapply(x, function(z) max(0,z))

# local smooth 
local_lowrank_smooth <- function(matrix_array,T_,t,h,p,q,lambda){
  weight <- kernel_weight(t,h,T_)
  weight <- weight/sum(weight)
  M <- matrix(0,p,q)
  len <- dim(matrix_array)[1]
  for (i in 1:len) {
    M <- M + weight[i]*matrix_array[i,,]
  }
  svd_M <- svd(M)
  N <- svd_M$u%*%diag(threshold_func(svd_M$d-lambda))%*%t(svd_M$v)
  return(N)
}

local_smooth <- function(matrix_array,T_,t,h,p,q,lambda){
  weight <- kernel_weight(t,h,T_)
  weight <- weight/sum(weight)
  M <- matrix(0,p,q)
  len <- dim(matrix_array)[1]
  for (i in 1:len) {
    M <- M + weight[i]*matrix_array[i,,]
  }
  return(M)
}
                                     
