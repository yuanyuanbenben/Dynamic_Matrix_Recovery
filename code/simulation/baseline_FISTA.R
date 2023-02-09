# using FISTA algorithm

# compute Lipschitz constant
Lipschitz_base_func <- function(X1,X2,n_t,p,q){
  L_mat <- matrix(0,p,q)
  for (i in 1:n_t) {
    L_mat[X1[i],X2[i]] <- L_mat[X1[i],X2[i]] + 1
  }
  return(2*norm(L_mat,type = "F"))
}

threshold_func <- function(x) sapply(x, function(z) max(0,z))
                                     
# compute M[x1,x2]
baseline_inner_pro_func_M <- function(X1,X2,M,n_t){
  ret_vec <- 1:n_t
  for (i in 1:n_t) {
    ret_vec[i] <- M[X1[i],X2[i]]
  }
  return(ret_vec)
}
                                     
baseline_obj_func <- function(Y,inner_pro){
  return(sum((Y-inner_pro)*(Y-inner_pro)))
}

# compute gradient
baseline_grad_func <- function(Y,inner_pro,X1,X2,n_t,p,q){
  # Y, inner_pro T*n_t
  # X T*n_t*p*q
  grad_mat <- matrix(0,p,q)
  tem_mat <- (inner_pro-Y)
  for (i in 1:n_t) {
    grad_mat[X1[i],X2[i]] <- grad_mat[X1[i],X2[i]] + tem_mat[i]
  }
  return(2*grad_mat)
}

# main function
baseline_FISTA_func <- function(X1,X2,Y,n_t,p,q,lambda,itertime=30000,constant=TRUE,tor=1,init=FALSE,M_input=FALSE){
  # initial
  L <- Lipschitz_base_func(X1,X2,n_t,p,q)
  #print(L)
  t = 1
  if (init){
    M = M_input
  }
  else{
    M = matrix(rnorm(p*q),p,q)
  }
  N <- M
  obj_value_before <- Inf
  # iteration
  for (iter in 1:itertime) {
    start <- Sys.time()
    inner_pro <- baseline_inner_pro_func_M(X1,X2,N,n_t)
    grad_N <- baseline_grad_func(Y,inner_pro,X1,X2,n_t,p,q)
    svd_G <- svd(N - 1/L*grad_N)
    M_ <- svd_G$u%*%diag(threshold_func(svd_G$d-lambda/L))%*%t(svd_G$v)
    t_ <- (1+sqrt(1+4*t^2))/2
    N <- M + (t-1)/t_*(M_-M)
    M <- M_
    t <- t_
    if (as.integer(iter/20)*20==iter){
      print(iter)
      obj_value <- baseline_obj_func(Y,inner_pro)
      print(obj_value)
      if (abs(obj_value-obj_value_before) < tor) break
      obj_value_before <- obj_value
    }
    end <- Sys.time()
    print(difftime(end, start, units = "sec"))
  }
  return(M)
}

