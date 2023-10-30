library(jpeg)
library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)

setwd('/your dictionary/Dynamic_Matrix_Recovery')

source("code/real_data/video_data/robust_pca.R")

# read pic
img_total_r <- array(0,dim = c(96,480,854))
img_total_g <- array(0,dim = c(96,480,854))
img_total_b <- array(0,dim = c(96,480,854))
for (index in 0:95) {
  if (index < 10){
    dic <- paste("data/lions/0000",index,".jpg",sep = "")
  }
  else{
    dic <- paste("data/lions/000",index,".jpg",sep = "")
  }
  img <- readJPEG(dic)
  img_total_r[index+1,,] <- img[,,1]
  img_total_g[index+1,,] <- img[,,2]
  img_total_b[index+1,,] <- img[,,3]
}

T_ = dim(img_total_b)[1]
p = dim(img_total_b)[2]
q = dim(img_total_b)[3]

# do RPCA to decomposition M to L and S
L_total = array(0,dim = c(3,T_,p,q))
S_total = array(0,dim = c(3,T_,p,q))


cl.cores = detectCores(logical = F)
cl <- makeCluster(21)
registerDoParallel(cl)

save = TRUE
foreach (i = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
  l = rpca_func(img_total_r[i,,],0.25,0.0125)
  if (save == TRUE){
    write.csv(l[[1]],paste("output/lions/lions_red_lowrank_",i,".csv",sep = ""))
    write.csv(l[[2]],paste("output/lions/lions_red_sparse_",i,".csv",sep = ""))
  }
}
foreach (i = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
  l = rpca_func(img_total_g[i,,],0.25,0.0125)
  if (save == TRUE){
    write.csv(l[[1]],paste("output/lions/lions_green_lowrank_",i,".csv",sep = ""))
    write.csv(l[[2]],paste("output/lions/lions_green_sparse_",i,".csv",sep = ""))
  }
}
foreach (i = 1:T_,.packages = c("psych","kernlab","npmr"),.verbose=TRUE) %dopar% {
  l = rpca_func(img_total_b[i,,],0.25,0.0125)
  if (save == TRUE){
    write.csv(l[[1]],paste("output/lions/lions_blue_lowrank_",i,".csv",sep = ""))
    write.csv(l[[2]],paste("output/lions/lions_blue_sparse_",i,".csv",sep = ""))
  }
}

stopCluster(cl)

for (t in 1:T_) {
#   print(t)
  m1 = matrix(0,p,q)
  m1_data = read.csv(paste("output/lions/lions_red_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m1[i,j] = m1_data[i,j]
    }
  }
  L_total[1,t,,] = m1
  m2 = matrix(0,p,q)
  m2_data = read.csv(paste("output/lions/lions_green_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m2[i,j] = m2_data[i,j]
    }
  }
  L_total[2,t,,] = m2
  m3 = matrix(0,p,q)
  m3_data = read.csv(paste("output/lions/lions_blue_lowrank_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m3[i,j] = m3_data[i,j]
    }
  }
  L_total[3,t,,] = m3
}

for (t in 1:T_) {
  # print(t)
  m1 = matrix(0,p,q)
  m1_data = read.csv(paste("output/lions/lions_red_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m1[i,j] = m1_data[i,j]
    }
  }
  S_total[1,t,,] = m1
  m2 = matrix(0,p,q)
  m2_data = read.csv(paste("output/lions/lions_green_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m2[i,j] = m2_data[i,j]
    }
  }
  S_total[2,t,,] = m2
  m3 = matrix(0,p,q)
  m3_data = read.csv(paste("output/lions/lions_blue_sparse_",t,".csv",sep = ""))[,2:855]
  for (i in 1:p) {
    for (j in 1:q) {
      m3[i,j] = m3_data[i,j]
    }
  }
  S_total[3,t,,] = m3
}

save.image("output/lions_video.RData")
