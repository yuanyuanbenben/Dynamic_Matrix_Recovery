# using FISTA algorithm


# object functions
base_batch_func <- function(X1,X2,Y,batch_size){
  len_k = length(X1)
  index <- sample(1:len_k,batch_size,replace = TRUE)
  X1 <- X1[index]
  X2 <- X2[index]
  Y <- Y[index]
  return(list(X1,X2,Y))
}


base_obj_func  <- function(X1,X2,Y,M,lambda){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  err = 0
  len_k <- length(X1)
  for (s in 1:len_k) {
    err <- err + (M[X1[s],X2[s]]-Y[s])^2
  }
  err <- err + lambda*nuclear(M)
  return(err)
}

# gradient function
base_grad_func <- function(minibatch,M,p,q){
  X1 <- minibatch[[1]]
  X2 <- minibatch[[2]]
  Y <- minibatch[[3]]
  grad_mat <- matrix(0,p,q)
  len_k <- length(X1)
  for (s in 1:len_k) {
    grad_mat[X1[s],X2[s]] <- M[X1[s],X2[s]]-Y[s]
  }
  return(2*grad_mat)
}


# compute Lipschitz constant
base_Lipschitz_func <- function(X1,X2,p,q){
  l =  length(X1)
  return(2*sqrt(l))
}

threshold_func <- function(x) sapply(x, function(z) max(0,z))
# main function
base_FISTA_func <- function(X1,X2,Y,T_,t_,h,p,q,lambda,itertime=30000,sto=FALSE,
                            batch_size=FALSE,init=TRUE,M_input=FALSE,
                            eta=1.2,constant=TRUE,tor=30){
  # input
  # X1,X2: row and col index for observations, Y: score
  # M_input: the initial matrix with row customers, col movies
  # p,q: row and col dimensions, T_:time points
  # h: bandwidth, lambda: truncation parameter
  
  
  L <- base_Lipschitz_func(X1,X2,p,q)*batch_size/length(X1[[1]])
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
    # print(iter)
    minibatch <- base_batch_func(X1,X2,Y,batch_size)
    grad_N <- base_grad_func(minibatch,N,p,q)
    svd_G <- svd(N - 1/L*grad_N)
    M_ <- svd_G$u%*%diag(threshold_func(svd_G$d-lambda/L))%*%t(svd_G$v)
    t_ <- (1+sqrt(1+4*t^2))/2
    N <- M + (t-1)/t_*(M_-M)
    M <- M_
    t <- t_
    if (as.integer(iter/20)*20==iter){
      print(paste('iteration',iter))
      obj_value <- base_obj_func(X1,X2,Y,M,lambda)
      print(paste('loss',obj_value))
      if (abs(obj_value-obj_value_before) < tor) break
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
