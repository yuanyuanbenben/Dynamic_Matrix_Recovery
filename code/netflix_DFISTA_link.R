# using FISTA algorithm

# kernel 
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

# object functions

batch_func <- function(X1,X2,Y,len,batch_size){
  index_time <- c(rmultinom(n=1,prob = rep(1/len,len),size=batch_size))
  for (k in 1:len) {
    len_k = length(X1[[k]])
    index <- sample(1:len_k,index_time[k],replace = TRUE)
    X1[[k]] <- X1[[k]][index]
    X2[[k]] <- X2[[k]][index]
    Y[[k]] <- Y[[k]][index]
  }
  return(list(X1,X2,Y))
}

variance_func <- function(a){
  return((a - floor(a))*(ceiling(a) - a))
}
grad_variance_func <- function(a){
  return(1 - 2*(a - floor(a)))
}

obj_func  <- function(X1,X2,Y,weight,M,len,lambda){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  err = 0
  for (k in 1:len) {
    len_k <- length(X1[[k]])
    err_k = 0
    for (s in 1:len_k) {
      m_i = M[X1[[k]][s],X2[[k]][s]]
      err_k <- err_k + (m_i-Y[[k]][s])^2 + variance_func(m_i)
    }
    err <- err + weight[k]*err_k
  }
  err <- err + lambda*nuclear(M)
  return(err)
}

# gradient function
grad_func <- function(minibatch,weight,M,len,p,q){
  X1 <- minibatch[[1]]
  X2 <- minibatch[[2]]
  Y <- minibatch[[3]]
  grad_mat <- matrix(0,p,q)
  for (k in 1:len) {
    len_k <- length(X1[[k]])
    for (s in 1:len_k) {
      m_i = M[X1[[k]][s],X2[[k]][s]]
      grad_mat[X1[[k]][s],X2[[k]][s]] <- weight[k]*(m_i - Y[[k]][s] + grad_variance_func(m_i))
    }
  }
  return(2*grad_mat)
}


# compute Lipschitz constant
Lipschitz_func <- function(X1,X2,len,weight,p,q){
  # L_mat <- matrix(0,p,q)
  l = 0
  w = weight^2
  for (j in 1:len) {
    l = l + length(X1[[j]])*w[j]
  }
  return(2*sqrt(l))
}

threshold_func <- function(x) sapply(x, function(z) max(0,z))
# main function
FISTA_func <- function(X1,X2,Y,T_,t_,h,p,q,lambda,itertime=1000,sto=FALSE,
                       batch_size=FALSE,init=TRUE,M_input=FALSE,
                       eta=1.2,constant=TRUE,tor=30){
  # input
  # X1,X2: row and col index for observations, Y: score
  # M_input: the initial matrix with row customers, col movies
  # p,q: row and col dimensions, T_:time points
  # h: bandwidth, lambda: truncation parameter
  
  # initial
  if (t_ == 1) {
    itertime = 3000
  }
  if (t_ < h/2) {
    tor = tor/2
  }
  len <- length(X1)
  weight <- kernel_weight(t_,h,T_)
  L <- Lipschitz_func(X1,X2,len,weight,p,q)*batch_size/len/length(X1[[1]])/10/3
  t = 1
  if (init){
    M = matrix(rnorm(p*q),p,q)
  }
  else{
    M = M_input
  }
  N <- M
  M_old <- M
  #iteration
  obj_value_before <- Inf
  for (iter in 1:itertime) {
    minibatch <- batch_func(X1,X2,Y,len,batch_size)
    grad_N <- grad_func(minibatch,weight,N,len,p,q)
    svd_G <- svd(N - 1/L*grad_N)
    M_ <- svd_G$u%*%diag(threshold_func(svd_G$d-lambda/L))%*%t(svd_G$v)
    t_ <- (1+sqrt(1+4*t^2))/2
    N <- M + (t-1)/t_*(M_-M)
    M <- M_
    t <- t_
    if (as.integer(iter/20)*20==iter){
      print(paste('iteration',iter))
      obj_value <- obj_func(X1,X2,Y,weight,M,len,lambda)
      print(paste('loss',obj_value))
      if (abs(obj_value-obj_value_before) < tor || obj_value-obj_value_before > 0) break
      obj_value_before <- obj_value
      M[M>5]=5
      M[M<1]=1
      N[N>5]=5
      N[N<1]=1
      M_old = M
    }
  }
  return(M_old)
}
