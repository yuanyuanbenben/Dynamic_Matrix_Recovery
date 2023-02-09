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
                                     
