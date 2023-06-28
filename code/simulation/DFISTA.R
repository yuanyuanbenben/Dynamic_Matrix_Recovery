# using FISTA algorithm

# compute M{x1,x2]
inner_pro_func_M <- function(X1,X2,M,len,n_t,sto=FALSE,minibatch=FALSE,batch_size=FALSE){
  # X T*n_t*p*q
  # M p*q
  ret_mat <- matrix(0,len,n_t)
  if (sto){
    for (index in 1:batch_size) {
      j <- minibatch[1,index]
      i <- minibatch[2,index]
      ret_mat[j,i] <- M[X1[j,i],X2[j,i]]
    }
  }
  else{
    for (j in 1:len) {
      for (i in 1:n_t) {
        ret_mat[j,i] <- M[X1[j,i],X2[j,i]]
      }
    }
  }
  return(ret_mat)
}

obj_func  <- function(weight_mat,Y,inner_pro){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  obj_mat <- sqrt(weight_mat)%*%(Y-inner_pro)
  return(norm(obj_mat,"F")^2)
}

# gradient function
grad_func <- function(weight_mat,Y,inner_pro,X1,X2,len,n_t,p,q,sto=FALSE,minibatch=FALSE,batch_size=FALSE){
  # weight_mat T*T
  # Y, inner_pro T*n_t
  # X T*n_t*p*q
  grad_mat <- matrix(0,p,q)
  tem_mat <- weight_mat%*%(inner_pro-Y)
  if (sto){
    for (index in 1:batch_size) {
      j <- minibatch[1,index]
      i <- minibatch[2,index]
      grad_mat[X1[j,i],X2[j,i]] <- grad_mat[X1[j,i],X2[j,i]] + tem_mat[j,i]
    }
    return(2*grad_mat)
  }
  else{
    for (j in 1:len) {
      for (i in 1:n_t) {
        grad_mat[X1[j,i],X2[j,i]] <- grad_mat[X1[j,i],X2[j,i]] + tem_mat[j,i]
      }
    }
    return(2*grad_mat)
  }
}


sample_func <- function(batch_size,weight_mat,len,n_t){
  size1 <- sample(1:len,batch_size,replace = TRUE)#,prob = diag(weight_mat))
  size2 <- sample(1:n_t,batch_size,replace = TRUE)
  return(rbind(size1,size2))
}

# compute Lipschitz constant
Lipschitz_func <- function(X1,X2,len,n_t,weight_mat,p,q){
  L_mat <- matrix(0,p,q)
  w <- diag(weight_mat)
  for (j in 1:len) {
    for (i in 1:n_t) {
      L_mat[X1[j,i],X2[j,i]] <- L_mat[X1[j,i],X2[j,i]] + w[j]
    }
  }
  return(2*norm(L_mat,type = "F"))
}

threshold_func <- function(x) sapply(x, function(z) max(0,z))
# main function
FISTA_func <- function(X1,X2,Y,T_,t_,h,p,q,lambda,itertime=3000,sto=FALSE,
                       batch_size=FALSE,init=TRUE,M_input=FALSE,
                       eta=1.2,constant=TRUE,tor=30){
  # initial
  if (t_==1) {
    itertime = 30000
  }
  if (t_<h/2) {
    tor = tor/2
  }
  len = dim(X1)[1]
  n_t = dim(X1)[2]
  weight_mat <- diag(kernel_weight(t_,h,T_))
  L <- Lipschitz_func(X1,X2,len,n_t,weight_mat,p,q)
  if (sto){
    L = L*batch_size/n_t/len
  }
  #print(L)
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
  if (sto){
    for (iter in 1:itertime) {
      #start <- Sys.time()
      minibatch <- sample_func(batch_size,weight_mat,len,n_t)
      inner_pro <- inner_pro_func_M(X1,X2,N,len,n_t,sto = sto,minibatch = minibatch,batch_size = batch_size)
      grad_N <- grad_func(weight_mat,Y,inner_pro,X1,X2,len,n_t,p,q,sto = sto,minibatch = minibatch,batch_size = batch_size)
      svd_G <- svd(N - 1/L*grad_N)
      M_ <- svd_G$u%*%diag(threshold_func(svd_G$d-lambda/L))%*%t(svd_G$v)
      t_ <- (1+sqrt(1+4*t^2))/2
      N <- M + (t-1)/t_*(M_-M)
      M <- M_
      t <- t_
      if (as.integer(iter/20)*20==iter){
        inner_pro_bas <- inner_pro_func_M(X1,X2,M,len,n_t)
        print(paste('iteration',iter))
        obj_value <- obj_func(weight_mat,Y,inner_pro_bas)
        print(paste('loss',obj_value))
        if (abs(obj_value-obj_value_before) < tor || obj_value-obj_value_before > 0) break
        obj_value_before <- obj_value
        M_old = M
      }
      #end <- Sys.time()
      #print(difftime(end, start, units = "sec"))
    }
    return(M)
  }
  else{
    for (iter in 1:itertime) {
      inner_pro <- inner_pro_func_M(X1,X2,N,len,n_t)
      grad_N <- grad_func(weight_mat,Y,inner_pro,X1,X2,len,n_t,p,q)
      svd_G <- svd(N - 1/L*grad_N)
      M_ <- svd_G$u%*%diag(threshold_func(svd_G$d-lambda/L))%*%t(svd_G$v)
      t_ <- (1+sqrt(1+4*t^2))/2
      N <- M + (t-1)/t_*(M_-M)
      M <- M_
      t <- t_
      if (as.integer(iter/20)*20==iter){
        inner_pro_bas <- inner_pro_func_M(X1,X2,M,len,n_t)
        print(paste('iteration',iter))
        obj_value <- obj_func(weight_mat,Y,inner_pro_bas)
        print(paste('loss',obj_value))
        if (abs(obj_value-obj_value_before) < tor || obj_value-obj_value_before > 0) break
        obj_value_before <- obj_value
        M_old = M
      }
    }
    return(M_old)
  }
}
