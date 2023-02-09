cs_realdata_createdata_func <- function(t,X1_total,X2_total,Y_total,T_,h,mode){
  h_ = as.integer(h/2)
  
  if (t - h_ <= 1){
    X1 <- X1_total[,1:(h_+t)]
    X2 <- X2_total[,1:(h_+t)]
    Y <- Y_total[mode,,1:(h_+t)]
  }
  if (t + h_ > T_){
    X1 <- X1_total[,(t-h_):T_]
    X2 <- X2_total[,(t-h_):T_]
    Y <- Y_total[mode,,(t-h_):T_]
  }
  if (t - h_ > 1 & t + h_ <= T_) {
    X1 <- X1_total[,(t-h_):(t+h_)]
    X2 <- X2_total[,(t-h_):(t+h_)]
    Y <- Y_total[mode,,(t-h_):(t+h_)]
  }
  return(list(X1,X2,Y))
}

cs_test_error_func <- function(M,N,p,q){
  return(norm(M-N,type="F")^2/p/q)
}
