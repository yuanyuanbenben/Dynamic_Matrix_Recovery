# pre process 

# read data from original dataset
setwd("/home/yuanyuanbenben/project_dmc/dmr_github/Dynamic_Matrix_Recovery-main/code")
data <- read.csv("real_data/netflix_data/netflix_data.csv")
names(data) <- c("movies","id","score","time")
dims <- dim(data)[1]
data$time <- as.Date(data$time)
movies <- as.data.frame(table(data$movies))

# select movies for recommendation systems (movies that rated more than 25000)
used_movies <- as.numeric(levels(droplevels(movies[movies$Freq>25000,1])))
data <- data[data$movies %in% used_movies,]
ids <- as.data.frame(table(data$id))

# select a subgroup users for completion (randomly select 3000 users)
# used_id <- as.numeric(levels(droplevels(ids[ids$Freq>700,1])))

used_id <- as.numeric(levels(droplevels(ids[ids$Freq>30,1])))
data <- data[data$id %in% used_id,]
set.seed(1228716052)
used_id <- sample(used_id,3000,replace=FALSE)
data <- data[data$id %in% used_id,]
 


#index_train <- sample(1:dim(data)[1],round(0.8*dim(data)[1]))
#train_data <- data[index_train,]
#train_data <- train_data[order(train_data$time),]
#test_data <- data[-index_train,]
#test_data <- test_data[order(test_data$time),]
#data <- rbind(train_data,test_data)
data <- data[order(data$time),]
data <- data[1:526500,]
# data <- data[1:2337500,]
#rm(data)
#rm(index_train)
#gc()
#
origin = "1999-12-30"
finish = "2005-12-30"
#finish_ = "2006-3-01"
#date_split = 22
T_ = 100
# input X1,X2,Y
len <- dim(data)[1]
len_ <- round(len/T_)
X1_total <- list()
X2_total <- list()
Y_total <- list()

index <- seq(0,len,len_)

# data splitting
for (i in 1:T_) {
  x1 <- rep(0,index[i+1] - index[i])
  x2 <- rep(0,index[i+1] - index[i])
  y <- rep(0,index[i+1] - index[i])
  j = 1
  for (k in (index[i]+1):index[i+1]) {
    x1[j] <- which(used_id == data$id[k])
    x2[j] <- which(used_movies == data$movies[k])
    y[j] <- data$score[k]
    j = j + 1
  }
  if (i==1){
    X1_total <- list(x1)
    X2_total <- list(x2)
    Y_total <- list(y)
  }
  else{
    X1_total <- c(X1_total,list(x1))
    X2_total <- c(X2_total,list(x2))
    Y_total <- c(Y_total,list(y))
  }
}

# train data and test data
train_X1_total <- list()
train_X2_total <- list()
train_Y_total <- list()
test_X1_total <- list()
test_X2_total <- list()
test_Y_total <- list()

set.seed(1278451)
for (i in 1:T_) {
  #index_train <- sample(1:23375,round(0.8*23375))
  index_train <- sample(1:5265,round(0.8*5265))
  train_X1_total[[i]] <- X1_total[[i]][index_train]
  test_X1_total[[i]] <- X1_total[[i]][-index_train]
  train_X2_total[[i]] <- X2_total[[i]][index_train]
  test_X2_total[[i]] <- X2_total[[i]][-index_train]
  train_Y_total[[i]] <- Y_total[[i]][index_train]
  test_Y_total[[i]] <- Y_total[[i]][-index_train]
}

# 5-fold cv

# set.seed(12471124)
# train_X1_total_cv1 <- list()
# train_X2_total_cv1 <- list()
# train_Y_total_cv1 <- list()

# train_X1_total_cv2 <- list()
# train_X2_total_cv2 <- list()
# train_Y_total_cv2 <- list()

# train_X1_total_cv3 <- list()
# train_X2_total_cv3 <- list()
# train_Y_total_cv3 <- list()

# train_X1_total_cv4 <- list()
# train_X2_total_cv4 <- list()
# train_Y_total_cv4 <- list()

# for (i in 1:T_) {
#   index_1 <- sample(1:18700,round(0.25*18700))
#   index_2 <- sample(1:14025,round(14025/3))
#   index_3 <- sample(1:9350,round(9350/2))
#   c0 <- train_X1_total[[i]] 
#   train_X1_total_cv1[[i]] <- c0[index_1]
#   c1 <- c0[-index_1]
#   train_X1_total_cv2[[i]] <- c1[index_2]
#   c2 <- c1[-index_2]
#   train_X1_total_cv3[[i]] <- c2[index_3]
#   c3 <- c2[-index_3]
#   train_X1_total_cv4[[i]] <- c3
  
#   c0 <- train_X2_total[[i]] 
#   train_X2_total_cv1[[i]] <- c0[index_1]
#   c1 <- c0[-index_1]
#   train_X2_total_cv2[[i]] <- c1[index_2]
#   c2 <- c1[-index_2]
#   train_X2_total_cv3[[i]] <- c2[index_3]
#   c3 <- c2[-index_3]
#   train_X2_total_cv4[[i]] <- c3
  
#   c0 <- train_Y_total[[i]] 
#   train_Y_total_cv1[[i]] <- c0[index_1]
#   c1 <- c0[-index_1]
#   train_Y_total_cv2[[i]] <- c1[index_2]
#   c2 <- c1[-index_2]
#   train_Y_total_cv3[[i]] <- c2[index_3]
#   c3 <- c2[-index_3]
#   train_Y_total_cv4[[i]] <- c3
# }
save.image("netflixdata.RData")
