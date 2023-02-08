# input data
matrix_data <- read.csv("~/project_dmc/final_version/datas/baseline_120000_matrix.csv")[,2:301]
matrix_array <- array(0,dim = c(100,500,300))
for (i in 1:100) {
  print(i)
  for (j in 1:500) {
    for (k in 1:300) {
      matrix_array[i,j,k] <- matrix_data[(i-1)*500+j,k]
    }
  }
}

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
  #svd_M <- svd(M)
  #N <- svd_M$u%*%diag(threshold_func(svd_M$d-lambda))%*%t(svd_M$v)
  return(M)
}

result_mse <- rep(0,T_)
h=24
lambda = 1
for (t in 1:T_) {
  M_in <- createdata_func(matrix_array,t,T_,h)
  M <- local_smooth(M_in,T_,t,h,p,q,lambda = lambda)
  #M <- matrix_array[t,,]
  M_ =  (cos(pi*t/2/T_)*U_0 + sin(pi*t/2/T_)*U_1)%*%(D_1 + 
                                                       t/T_*D_0)%*%t(cos(pi*t/2/T_)*V_0 + sin(pi*t/2/T_)*V_1)
  result_mse[t] <- compare_matrix_func(M,M_,p,q)[[1]]
  print(result_mse[t])
}
write.csv(result_mse,"~/project_dmc/final_version/datas/localsmooth_120000_2.csv")
write.csv(result_mse,"~/final_version/datas/baseline_120000.csv")
