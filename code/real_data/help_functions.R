# error functions
test_error_func <- function(M,X1,X2,Y){
  # X1,X2,Y: vector
  len = length(X1)
  err = 0
  for (i in 1:len) {
    err <- err + (M[X1[i],X2[i]] - Y[i])^2
  }
  return(err/len)
}

realdata_createdata_func <- function(t,X1_total,X2_total,Y,T_,h){
  h_ = as.integer(h/2)
  
  if (t - h_ <= 1){
    X1 <- X1_total[1:(h_+t)]
    X2 <- X2_total[1:(h_+t)]
    Y <- Y[1:(h_+t)]
  }
  if (t + h_ > T_){
    X1 <- X1_total[(t-h_):T_]
    X2 <- X2_total[(t-h_):T_]
    Y <- Y[(t-h_):T_]
  }
  if (t - h_ > 1 & t + h_ <= T_) {
    X1 <- X1_total[(t-h_):(t+h_)]
    X2 <- X2_total[(t-h_):(t+h_)]
    Y <- Y[(t-h_):(t+h_)]
  }
  return(list(X1,X2,Y))
}

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


construc_pic <- function(t,l,mode,S_total){
  s = S_total[mode,t,,]
  m = s+l
  m[m<0]=0
  m[m>1]=1
  return(m)
}

construc_func <- function(t,versions){
  if (versions == "static"){
    m_data_red = read.csv(paste("~/Dynamic_Matrix_Recovery/output/baseline_lions_red_",t,".csv",sep = ""))[,2:855]
    m_data_green = read.csv(paste("~/Dynamic_Matrix_Recovery/output/baseline_lions_green_",t,".csv",sep = ""))[,2:855]
    m_data_blue = read.csv(paste("~/Dynamic_Matrix_Recovery/output/baseline_lions_blue_",t,".csv",sep = ""))[,2:855]
  }
  else if (versions == "twostep"){
    m_data_red = read.csv(paste("~/Dynamic_Matrix_Recovery/output/twostep_lions_red_",t,".csv",sep = ""))[,2:855]
    m_data_green = read.csv(paste("~/Dynamic_Matrix_Recovery/output/twostep_lions_reen_",t,".csv",sep = ""))[,2:855]
    m_data_blue = read.csv(paste("~/Dynamic_Matrix_Recovery/output/twostep_lions_blue_",t,".csv",sep = ""))[,2:855]
  }
  else {
    m_data_red = read.csv(paste("~/Dynamic_Matrix_Recovery/output/lions_red_",t,".csv",sep = ""))[,2:855]
    m_data_green = read.csv(paste("~/Dynamic_Matrix_Recovery/output/lions_green_",t,".csv",sep = ""))[,2:855]
    m_data_blue = read.csv(paste("~/Dynamic_Matrix_Recovery/output/lions_blue_",t,".csv",sep = ""))[,2:855]
  }
  m_red = matrix(0,p,q)
  m_green = matrix(0,p,q)
  m_blue = matrix(0,p,q)
  for (i in 1:p) {
    for (j in 1:q) {
      m_red[i,j] = m_data_red[i,j]
      m_green[i,j] = m_data_green[i,j]
      m_blue[i,j] = m_data_blue[i,j]
    }
  }
  return(list(m_red,m_green,m_blue))
}
