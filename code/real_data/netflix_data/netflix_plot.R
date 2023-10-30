# plots 
library(latex2exp)
library(ggplot2)
library(reshape2)
library(viridis)
setwd("/your dictionary/Dynamic_Matrix_Recovery")

source("real_data/help_functions.R")


origin = "1999-12-30"
finish = "2005-09-08"
finish_ = "2005-11-04"
T_ = 1000

mse_netflix_t_1 <- read.csv("output/netflix_mse_diff_t_1.csv")[,2]
mse_netflix_t_2 <- read.csv("output/netflix_mse_diff_t_2.csv")[,2]
mse_netflix_t_5 <- read.csv("output/netflix_mse_diff_t_5.csv")[,2]
mse_netflix_t_10 <- read.csv("output/netflix_mse_diff_t_10.csv")[,2]
mse_netflix_t_20 <- read.csv("output/netflix_mse_diff_t_20.csv")[,2]
mse_netflix_t_50 <- read.csv("output/netflix_mse_diff_t_50.csv")[,2]
mse_netflix_t_100 <- read.csv("output/netflix_mse_diff_t_100.csv")[,2]
mse_netflix_t_200 <- read.csv("output/netflix_mse_diff_t_200.csv")[,2]
mse_netflix_t_500 <- read.csv("output/netflix_mse_diff_t_500.csv")[,2]
mse_netflix_t_1000 <- read.csv("output/netflix_mse_diff_t_1000.csv")[,2]
for (i in 1:100){
  mse_netflix_t_1[i]<- mean(mse_netflix_t_1[((i-1)*10+1):(i*10)])
  mse_netflix_t_2[i]<- mean(mse_netflix_t_2[((i-1)*10+1):(i*10)])
  mse_netflix_t_5[i]<- mean(mse_netflix_t_5[((i-1)*10+1):(i*10)])
  mse_netflix_t_10[i]<- mean(mse_netflix_t_10[((i-1)*10+1):(i*10)])
  mse_netflix_t_20[i]<- mean(mse_netflix_t_20[((i-1)*10+1):(i*10)])
  mse_netflix_t_50[i]<- mean(mse_netflix_t_50[((i-1)*10+1):(i*10)])
  mse_netflix_t_100[i]<- mean(mse_netflix_t_100[((i-1)*10+1):(i*10)])
  mse_netflix_t_200[i]<- mean(mse_netflix_t_200[((i-1)*10+1):(i*10)])
  mse_netflix_t_500[i]<- mean(mse_netflix_t_500[((i-1)*10+1):(i*10)])
  mse_netflix_t_1000[i]<- mean(mse_netflix_t_1000[((i-1)*10+1):(i*10)])
}
mse_netflix_t_1<- mse_netflix_t_1[1:100]
mse_netflix_t_2<- mse_netflix_t_2[1:100]
mse_netflix_t_5<- mse_netflix_t_5[1:100]
mse_netflix_t_10<- mse_netflix_t_10[1:100]
mse_netflix_t_20<- mse_netflix_t_20[1:100]
mse_netflix_t_50<- mse_netflix_t_50[1:100]
mse_netflix_t_100<- mse_netflix_t_100[1:100]
mse_netflix_t_200<- mse_netflix_t_200[1:100]
mse_netflix_t_500<- mse_netflix_t_500[1:100]
mse_netflix_t_1000<- mse_netflix_t_1000[1:100]

mse_data <- data.frame("t"=seq.Date(as.Date(origin),as.Date(finish),by=21),"t1"=mse_netflix_t_1,"t2"=mse_netflix_t_2,"t5"=mse_netflix_t_5,
                       "t10"=mse_netflix_t_10,"t20"=mse_netflix_t_20,"t50"=mse_netflix_t_50,
                       "t100"=mse_netflix_t_100,"t200"=mse_netflix_t_200,"t500"=mse_netflix_t_500,"t1000"=mse_netflix_t_1000)
mydata <- melt(mse_data,id="t")
cc_ <- viridis::viridis(15)
cc = rep(0,10)
for (i in 1:10){
  cc[i] <- cc_[13-i]
}

a<-ggplot(data = mydata,aes(x=t,y=value,group=variable,
                         color=variable))+
  geom_line(size=0.5)+scale_colour_manual(values=cc,name="size of T",
                                          labels=c("T=1","T=2","T=5","T=10","T=20","T=50","T=100","T=200","T=500","T=1000"))+
  scale_x_continuous(name="Date",limits = c(as.Date(origin),as.Date(finish_)),breaks =seq.Date(as.Date(origin),as.Date(finish_),by=172))+
  scale_y_continuous(name="Mean Square Error",limits=c(0.6,1.6),breaks = seq(0.6,1.6,0.2))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.55,.9),legend.direction = "horizontal",
        legend.box.background = element_rect(color="black"),
        axis.text.x = element_text(angle = 70, hjust = 0, vjust = 0,size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))
ggsave("manuscript/netflix_diff_t.png",plot=a,width=12,height=8)
