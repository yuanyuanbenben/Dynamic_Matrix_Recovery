# using FISTA algorithm


# object function
baseline_inner_pro_func_M <- function(M,X1,X2,Y,n,d,v,minibatch,batch_size,conv_ker){
  # X T*n_t*p*q
  # M p*q
  ret <- rep(0,batch_size)
  for (index in 1:batch_size) {
    i <- minibatch[index]
    s <- X1[i]
    t <- X2[i]
    ret[index] <- sum(M[c(s-1,s,s+1),c(t-1,t,t+1)]*conv_ker) - Y[i]
  }
  return(ret)
}

baseline_obj_func  <- function(M,inner_pro,lambda){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  obj <- sum(inner_pro*inner_pro)
  return(obj + lambda*nuclear(M))
}

# gradient function
baseline_grad_func <- function(X1,X2,inner_pro,minibatch,batch_size,d,v,p,q,conv_ker){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  # X T*n_t*p*q
  grad_mat <- matrix(0,p,q)
  for (index in 1:batch_size) {
    i <- minibatch[index]
    s <- X1[i]
    t <- X2[i]
    grad_mat[c(s-1,s,s+1),c(t-1,t,t+1)] <- grad_mat[c(s-1,s,s+1),c(t-1,t,t+1)] + inner_pro[index]*conv_ker
  }
  return(2*grad_mat)
}


baseline_sample_func <- function(batch_size,n){
  size <- sample(1:n,batch_size,replace = TRUE)
  return(size)
}

# compute Lipschitz constant
baseline_Lipschitz_func <- function(X1,X2,n,p,q,conv_ker){
  L_mat <- matrix(0,p,q)
  for (i in 1:n) {
    s = X1[i]
    t = X2[i]
    L_mat[c(s-1,s,s+1),c(t-1,t,t+1)] <- L_mat[c(s-1,s,s+1),c(t-1,t,t+1)] + conv_ker
  }
  return(2*6*norm(L_mat,type = "F"))
}

# main function
baseline_cs_FISTA_func <- function(X1,X2,Y,T_,t_,p,q,lambda,conv_ker,itertime=35000,sto=FALSE,
                          batch_size=FALSE,init=TRUE,M_input=FALSE,
                          eta=1.2,constant=TRUE,tor=100){
  # initial
  if (t_ == 1) {
    itertime = 20000
  }
  # if (t_ < h/2) {
  #   tor = tor/2
  #  }
  n = length(X1)
  L <- baseline_Lipschitz_func(X1,X2,n,p,q,conv_ker)
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
    minibatch <- baseline_sample_func(batch_size,n)
    inner_pro <- baseline_inner_pro_func_M(N,X1,X2,Y,n,d,v,minibatch,batch_size,conv_ker)
    grad_N <- baseline_grad_func(X1,X2,inner_pro,minibatch,batch_size,d,v,p,q,conv_ker)
    svd_G <- svd(N - 1/L*grad_N)
    M_ <- svd_G$u%*%diag(threshold_func(svd_G$d,lambda/L))%*%t(svd_G$v)
    t_ <- (1+sqrt(1+4*t^2))/2
    N <- M + (t-1)/t_*(M_-M)
    M <- M_
    t <- t_
    if (as.integer(iter/20)*20==iter){
      print(iter)
      obj_value <- baseline_obj_func(M,inner_pro,lambda)
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
