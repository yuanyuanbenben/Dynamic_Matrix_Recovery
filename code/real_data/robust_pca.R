# implement robust pca 
# M = L + S

# shrinkage function D
threshold_func <- function(x,mu) sapply(x, function(z) max(0,z-mu),simplify = "array")

shrinkage_function <- function(X,rho){
  svd_X <- svd(X)
  return(svd_X$u%*%(diag(threshold_func(svd_X$d,rho)))%*%t(svd_X$v))
}

# truncation function S
truncation_function <- function(X,mu,p,q) {
  matrix(sapply(X, function(z) sign(z)*max(0,abs(z)-mu),simplify = "array"),p,q)
}

# robust pca main function
rpca_func <- function(M,mu=FALSE,lambda=FALSE,itertime=30){
  # init
  p = dim(M)[1]
  q = dim(M)[2]
  S = matrix(0,p,q)
  Y = matrix(0,p,q)
  if (!mu){
    mu = p*q/4/norm(M,"1")
  }
  if (!lambda){
    lambda = 1/sqrt(max(p,q))
  }
  for (i in 1:itertime) {
    L <- shrinkage_function(M - S + 1/mu*Y,1/mu)
    S <- truncation_function(M - L + 1/mu*Y,lambda/mu,p,q)
    Y <- Y + mu*(M-L-S)
   # print(i)
   # print(norm(Y))
  }
  return(list(L,S,Y))
}
