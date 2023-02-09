# tensor completion using HaLRTC 

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
    #print(iter)
    begin <- Sys.time()
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
    print(iter)
    print(difftime(end, begin, units = "sec"))
    #print(mean_time_mse(M,Z,T_,p,q))
  }
  return(Z)
}



HaLRTC <- function(X,Y,a,T_,p,q,rho=1,itertime=1000){
  # Y target tensor
  # X treatment index 
  # a weight
  S1 = array(0,dim = c(T_,p,q))
  S2 = array(0,dim = c(T_,p,q))
  S3 = array(0,dim = c(T_,p,q))
  for (iter in 1:itertime) {
    print(iter)
    M1 = fold_tensor(shrinkage_function(unfold_tensor(Y,mode=1)+unfold_tensor(S1,mode=1)/rho,a[1]/rho),mode=1,shape=c(T_,p,q))
    M2 = fold_tensor(shrinkage_function(unfold_tensor(Y,mode=2)+unfold_tensor(S1,mode=2)/rho,a[2]/rho),mode=2,shape=c(T_,p,q))
    M3 = fold_tensor(shrinkage_function(unfold_tensor(Y,mode=3)+unfold_tensor(S1,mode=3)/rho,a[3]/rho),mode=3,shape=c(T_,p,q))
    Y_hat = (1-X)*(M1+M2+M3-(S1+S2+S3)/rho)/3+Y
    S1 = S1-rho*(M1-Y_hat)
    S2 = S2-rho*(M2-Y_hat)
    S3 = S3-rho*(M3-Y_hat)
  }
  return(Y_hat)
}

write.csv(train_data,'tensor_netflix_traindata.csv')
write.csv(test_data,'tensor_netflix_testdata.csv')

# X1_total,X2_total,Y
X_tensor <- array(FALSE,dim = c(T_,p,q))
Y_tensor <- array(0,dim = c(T_,p,q))
for (j in 1:T_) {
  index1 <- X1_total[[j]]
  index2 <- X2_total[[j]]
  lens <- length(index1)
  for (i in 1:lens) {
    X_tensor[j,index1[i],index2[i]] <- TRUE
    Y_tensor[j,index1[i],index2[i]] <- Y[[j]][i]
  }
}

X_tensor_test <- array(FALSE,dim = c(T_,p,q))
Y_tensor_test <- array(0,dim = c(T_,p,q))
for (j in 1:T_) {
  index1 <- test_X1_total[[j]]
  index2 <- test_X2_total[[j]]
  lens <- length(index1)
  for (i in 1:lens) {
    X_tensor_test[j,index1[i],index2[i]] <- TRUE
    Y_tensor_test[j,index1[i],index2[i]] <- test_Y[[j]][i]
  }
}
#a = c(1/3,1/3,1/3)
#rho=10
Y_tensor_hat <- ADM_TR(X_tensor,Y_tensor,T_,p,q,beta = 0.1,lamda=0.3,c_beta=1,c_lamda=1,itertime = 1000)
#M_hat <- HaLRTC(X,Y,a,T_,p,q,rho=rho,itertime = 10000)
write.csv(Y_tensor_hat,'~/final_version/datas/tensor_netflix.csv')
Y_tensor_test_hat <- Y_tensor_hat[X_tensor_test]
mean_time_mse <- function(Y_tensor_test,X_tensor_test,Y_tensor_test_hat,T_,p,q){
  err <- rep(0,T_)
  for (i in 1:T_) {
    err[i] = sum((Y_tensor_test[i,,]-Y_tensor_test_hat[i,,])^2)/sum(X_tensor_test[i,,])
  }
  return(err)
}

mean_time_mse(Y_tensor_test,X_tensor_test,Y_tensor_test_hat,T_,p,q)


