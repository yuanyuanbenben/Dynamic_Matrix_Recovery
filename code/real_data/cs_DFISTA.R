# using FISTA algorithm for compressed sensing problem
# real data example 2: compress and recover the lion video from Davis dataset 

# kernel 
kernel_weight<- function(t,h,T_){
  epan_ker <- function(x,h){
    if (abs(x) > as.integer(h/2)){
      return(0)
    }
    else {
      #return(3/4*(1-4*x^2/h^2))
      return((as.integer(h/2)+1-abs(x))^3*4.95/45)
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


# object function
inner_pro_func_M <- function(M,X1,X2,Y,len,n,d,v,minibatch,batch_size,conv_ker){
  # X T*n_t*p*q
  # M p*q
  ret <- rep(0,batch_size)
  for (index in 1:batch_size) {
    j <- minibatch[1,index]
    i <- minibatch[2,index]
    s <- X1[i,j]
    t <- X2[i,j]
    ret[index] <- sum(M[c(s-1,s,s+1),c(t-1,t,t+1)]*conv_ker) - Y[i,j]
  }
  return(ret)
}

obj_func  <- function(M,inner_pro,minibatch,batch_size,weight,lambda){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  obj <- 0
  for (index in 1:batch_size) {
    obj <- obj + weight[minibatch[1,index]]*(inner_pro[index]^2)
  }
  return(obj + lambda*nuclear(M))
}

# gradient function
grad_func <- function(X1,X2,inner_pro,minibatch,batch_size,d,v,weight,p,q,conv_ker){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  # X T*n_t*p*q
  grad_mat <- matrix(0,p,q)
  for (index in 1:batch_size) {
    j <- minibatch[1,index]
    i <- minibatch[2,index]
    s <- X1[i,j]
    t <- X2[i,j]
    grad_mat[c(s-1,s,s+1),c(t-1,t,t+1)] <- grad_mat[c(s-1,s,s+1),c(t-1,t,t+1)] + weight[j]*inner_pro[index]*conv_ker
  }
  return(2*grad_mat)
}


sample_func <- function(batch_size,weight,len,n){
  size1 <- sample(1:len,batch_size,replace = TRUE)#,prob = diag(weight_mat))
  size2 <- sample(1:n,batch_size,replace = TRUE)
  return(rbind(size1,size2))
}

# compute Lipschitz constant
Lipschitz_func <- function(X1,X2,len,n,weight,p,q,conv_ker){
  L_mat <- matrix(0,p,q)
  for (j in 1:len) {
    for (i in 1:n) {
      s = X1[i,j]
      t = X2[i,j]
      L_mat[c(s-1,s,s+1),c(t-1,t,t+1)] <- L_mat[c(s-1,s,s+1),c(t-1,t,t+1)] + weight[j]*conv_ker
    }
  }
  return(2*6*norm(L_mat,type = "F"))
}

# main function
cs_FISTA_func <- function(X1,X2,Y,T_,t_,h,p,q,lambda,conv_ker,itertime=35000,sto=FALSE,
                          batch_size=FALSE,init=TRUE,M_input=FALSE,
                          eta=1.2,constant=TRUE,tor=100){
  # initial
  if (t_ == 1) {
    itertime = 20000
  }
  # if (t_ < h/2) {
  #   tor = tor/2
 #  }
  n = dim(X1)[1]
  len = dim(X1)[2]
  weight <- kernel_weight(t_,h,T_)
  L <- Lipschitz_func(X1,X2,len,n,weight,p,q,conv_ker)
  print(L)
  t = 1
  if (init){
    M = matrix(rnorm(p*q),p,q)
  }
  else{
    M = M_input
  }
  N <- M
  #iteration
  obj_value_before <- Inf
  for (iter in 1:itertime) {
    # print(iter)
    minibatch <- sample_func(batch_size,weight,len,n)
    inner_pro <- inner_pro_func_M(N,X1,X2,Y,len,n,d,v,minibatch,batch_size,conv_ker)
    grad_N <- grad_func(X1,X2,inner_pro,minibatch,batch_size,d,v,weight,p,q,conv_ker)
    svd_G <- svd(N - 1/L*grad_N)
    M_ <- svd_G$u%*%diag(threshold_func(svd_G$d,lambda/L))%*%t(svd_G$v)
    t_ <- (1+sqrt(1+4*t^2))/2
    N <- M + (t-1)/t_*(M_-M)
    M <- M_
    t <- t_
    if (as.integer(iter/20)*20==iter){
      print(iter)
      obj_value <- obj_func(M,inner_pro,minibatch,batch_size,weight,lambda)
      print(obj_value)
      if (abs(obj_value-obj_value_before) < tor) break
      obj_value_before <- obj_value
      N[N>1]=0.98
      N[N<0]=0.02
      M[M>1]=0.98
      M[M<0]=0.02
    }
  }
  return(M)
}
