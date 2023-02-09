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

mean_time_mse <- function(Y_tensor_test,X_tensor_test,Y_tensor_test_hat,T_,p,q){
  err <- rep(0,T_)
  for (i in 1:T_) {
    err[i] = sum((Y_tensor_test[i,,]-Y_tensor_test_hat[i,,])^2)/sum(X_tensor_test[i,,])
  }
  return(err)
}
threshold_func <- function(x) sapply(x, function(z) max(0,z))

shrinkage_function <- function(X,rho){
  svd_X <- svd(X)
  return(svd_X$u%*%(diag(threshold_func(svd_X$d-rho)))%*%t(svd_X$v))
}

unfold_tensor <- function(X,mode){
  dim_X <- dim(X)
  if (mode==1){
    ret_mat <- matrix(0,nrow = dim_X[1],ncol = dim_X[2]*dim_X[3])
    for (i in 1:dim_X[2]) {
      for (j in 1:dim_X[3]) {
        ret_mat[,(i-1)*dim_X[3]+j] <- X[,i,j]
      }
    }
    return(ret_mat)
  }
  if (mode==2){
    ret_mat <- matrix(0,nrow = dim_X[2],ncol = dim_X[3]*dim_X[1])
    for (i in 1:dim_X[3]) {
      for (j in 1:dim_X[1]) {
        ret_mat[,(i-1)*dim_X[1]+j] <- X[j,,i]
      }
    }
    return(ret_mat)
  }
  if (mode==3){
    ret_mat <- matrix(0,nrow = dim_X[3],ncol = dim_X[1]*dim_X[2])
    for (i in 1:dim_X[1]) {
      for (j in 1:dim_X[2]) {
        ret_mat[,(i-1)*dim_X[2]+j] <- X[i,j,]
      }
    }
    return(ret_mat)
  }
}

fold_tensor <- function(X,mode,shape){
  ret_tensor <- array(0,dim = shape)
  if (mode==1){
    for (i in 1:shape[2]){
      for (j in 1:shape[3]) {
        ret_tensor[,i,j] <- X[,(i-1)*shape[3]+j]
      }
    }
    return(ret_tensor)
  }
  if (mode==2){
    for (i in 1:shape[3]){
      for (j in 1:shape[1]) {
        ret_tensor[j,,i] <- X[,(i-1)*shape[1]+j]
      }
    }
    return(ret_tensor)
  }
  if (mode==3){
    for (i in 1:shape[1]){
      for (j in 1:shape[2]) {
        ret_tensor[i,j,] <- X[,(i-1)*shape[2]+j]
      }
    }
    return(ret_tensor)
  }
}

ADM_TR <- function(X,Y,T_,p,q,beta=0.1,lamda=0.1,c_beta=1,c_lamda=1,itertime=100){
  # initial
  Z = Y
  Y1 = Y
  Y2 = Y
  Y3 = Y
  W1 = array(0,dim = c(T_,p,q))
  W2 = array(0,dim = c(T_,p,q))
  W3 = array(0,dim = c(T_,p,q))
  
  for (iter in 1:itertime) {
    # update Z
    # begin <- Sys.time()
    Z[X] = 1/(lamda + 3*beta)*(W1+W2+W3+beta*(Y1+Y2+Y3)+lamda*Y)[X]
    Z[!X] = 1/(3*beta)*(W1+W2+W3+beta*(Y1+Y2+Y3))[!X]
    # update stepsize
    # update Y and W
    Y1 = fold_tensor(shrinkage_function(unfold_tensor(Z,mode = 1) - 
                                          1/beta*unfold_tensor(W1,mode = 1),rho = 1/beta),mode = 1,shape = c(T_,p,q))
    Y2 = fold_tensor(shrinkage_function(unfold_tensor(Z,mode = 2) - 
                                          1/beta*unfold_tensor(W2,mode = 2),rho = 1/beta),mode = 2,shape = c(T_,p,q))
    Y3 = fold_tensor(shrinkage_function(unfold_tensor(Z,mode = 3) - 
                                          1/beta*unfold_tensor(W3,mode = 3),rho = 1/beta),mode = 3,shape = c(T_,p,q))
    
    W1 = W1 - beta*(Z - Y1)
    W2 = W2 - beta*(Z - Y2)
    W3 = W3 - beta*(Z - Y3)
    beta = beta*c_beta
    lamda = lamda*c_lamda
    end <- Sys.time()
    # print(iter)
    # print(difftime(end, begin, units = "sec"))
  }
  return(Z)
}
