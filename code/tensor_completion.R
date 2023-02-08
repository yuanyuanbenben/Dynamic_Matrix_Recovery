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
    print(iter)
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
    print(difftime(end, begin, units = "sec"))
    print(mean_time_mse(M,Z,T_,p,q))
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


T_ = 100
n_t = 30000
p = 500
q = 300
r = 10
set.seed(1246326361)
N = matrix(rnorm(p*q),p,q)
svd_N <- svd(N)
U_0 <- svd_N$u[,1:r]
U_1 <- svd_N$u[,(r+1):(2*r)]
V_0 <- svd_N$v[,1:r]
V_1 <- svd_N$v[,(r+1):(2*r)]
D_0 <- diag(r:1)
D_1 <- D_0*D_0*10
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
#a = c(1/3,1/3,1/3)
#rho=10
M_hat <- ADM_TR(X,Y,T_,p,q,beta = 0.1,lamda=0.3,c_beta=1,c_lamda=1,itertime = 3000)
#M_hat <- HaLRTC(X,Y,a,T_,p,q,rho=rho,itertime = 10000)

mean_time_mse <- function(M,M_hat,T_,p,q){
  return(sum((M-M_hat)*(M-M_hat))/T_/p/q)
}

mean_time_mse(M,M_hat,T_,p,q)



# test

l1 = 20
l2 = 30
l3 = 40
r = 2
set.seed(12108432)
ker <- array(rnorm(r^3),dim = c(r,r,r))*10
s1 <- matrix(rnorm(l1*l2),l1,l2)
psi1 <- svd(s1)$u[,1:r]
psi2 <- svd(s1)$v[,1:r]
s2 <- matrix(rnorm(l2*l3),l2,l3)
psi3 <- svd(s2)$v[,1:r]
M <- array(0,dim = c(l1,l2,l3))

compos_func <- function(v1,v2,v3,ker){
  shape = dim(ker)
  m <- array(0,dim = shape[2:3])
  for (j in 1:shape[2]) {
    for (k in 1:shape[3]) {
      m[j,k] <- sum(ker[,j,k]*v1)
    }
  }
  n <- rep(0,shape[2])
  for (k in 1:shape[2]) {
    n[k] <- sum(m[,k]*v2)
  }
  return(sum(n*v3))
}

for (i in 1:l1) {
  for (j in 1:l2) {
    for (k in 1:l3) {
      M[i,j,k] <- compos_func(psi1[i,],psi2[j,],psi3[k,],ker)
    }
  }
}

X <- array(FALSE,dim = c(l1,l2,l3))
Y <- array(0,dim = c(l1,l2,l3))
rho <- 0.6
n_t <- round(l1*l2*l3*rho)
index1 <- sample(1:l1,n_t,replace = TRUE)
index2 <- sample(1:l2,n_t,replace = TRUE)
index3 <- sample(1:l3,n_t,replace = TRUE)
for (i in 1:n_t) {
  X[index1[i],index2[i],index3[i]] <- TRUE
  Y[index1[i],index2[i],index3[i]] <- M[index1[i],index2[i],index3[i]] + rnorm(1,sd=1)
}
M_hat <- ADM_TR(X,Y,l1,l2,l3,beta = 1,lamda=3,c_beta=1,c_lamda=1,itertime = 1000)
M_hat <- HaLRTC(X,Y,a,l1,l2,l3,rho=3,itertime = 1000)
mean_time_mse(M,M_hat,l1,l2,l3)
mean_time_mse(M,Y,l1,l2,l3)



