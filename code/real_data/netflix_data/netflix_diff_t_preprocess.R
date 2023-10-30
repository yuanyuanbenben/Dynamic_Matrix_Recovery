# pre process 

# read data from original dataset
setwd("/your dictionary/Dynamic_Matrix_Recovery")
# download original data
data <- read.csv("data/netflix_data.csv")


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
data <- data[1:500000,]
# data <- data[1:2337500,]
#rm(data)
#rm(index_train)
#gc()
#
origin = "1999-12-30"
finish = "2005-11-04"
# T = 10,20,50,100,200,500,1000
T_ = 1000
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
  index_train <- sample(1:500,round(0.8*500))
  train_X1_total[[i]] <- X1_total[[i]][index_train]
  test_X1_total[[i]] <- X1_total[[i]][-index_train]
  train_X2_total[[i]] <- X2_total[[i]][index_train]
  test_X2_total[[i]] <- X2_total[[i]][-index_train]
  train_Y_total[[i]] <- Y_total[[i]][index_train]
  test_Y_total[[i]] <- Y_total[[i]][-index_train]
}



save.image("data/diff_t.RData")
