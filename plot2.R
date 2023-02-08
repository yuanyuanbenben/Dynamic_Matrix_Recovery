# plots 
library(latex2exp)
library(ggplot2)
library(reshape2)
# read data
mse_baseline <- read.csv("~/project_dmc/final_version/datas/baseline_120000.csv")[,2]
mse_tensor <- read.csv("~/project_dmc/final_version/datas/tensor_30000.csv")[,2]
mse_dynamic <- read.csv("~/project_dmc/final_version/datas/dmc_5000_30000.csv")[,7]
mse_local_baseline <- read.csv("~/project_dmc/final_version/datas/localsmooth_120000.csv")[,2]
mse_data <- data.frame("t"=seq(0.01,1,0.01)*100,"benchmark1"=mse_baseline, "benchmark2" = mse_local_baseline,
                       "benchmark3"=mse_tensor, "DLRTR"=mse_dynamic)
mydata <- melt(mse_data,id="t")
ggplot(data = mydata,aes(x=t,y=value,group = variable,
                         color=variable))+
  geom_line(size=1)+
  #scale_size_manual(values=c(10,10,10,10))+
  scale_colour_discrete(name="Estimator",
                        labels=c(expression(paste("Static: ",rho,"=0.8",sep="")),
                                 expression(paste("TwoStep:",rho,"=0.8",sep="")),
                                 expression(paste("Tensor: ",rho,"=0.2",sep="")),
                                 expression(paste("DLR:  ",rho,"=0.2",sep=""))))+
  scale_x_continuous(name="t",limits = c(0,100),breaks = seq(0,100,10))+
  scale_y_continuous(name="Mean Square Error",limits=c(0,0.7),breaks = seq(0,0.7,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.55,.85),
        legend.box.background = element_rect(color="black"),
              axis.title.x=element_text(size=18),
              axis.title.y=element_text(size=18),
              axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))


#
mse = read.csv("~/project_dmc/final_version/datas/phase_transition.csv")
mse_10000_300 = mse[1,seq(2,301,6)]
mse_15000_200 = mse[2,seq(2,201,4)]
mse_30000_100 = mse[3,seq(2,101,2)]
mse_5000_300 = mse[4,seq(2,301,6)]
mse_7500_200 = mse[5,seq(2,201,4)]
mse_15000_100 = mse[6,seq(2,101,2)]
mse_5000_150 = mse[7,seq(2,151,3)]
mse_7500_100 = mse[8,seq(2,101,2)]
mse_15000_50 = mse[9,seq(2,51,1)]
mse1 = rep(0,50)
mse2 = rep(0,50)
mse3 = rep(0,50)
mse4 = rep(0,50)
mse5 = rep(0,50)
mse6 = rep(0,50)
mse7 = rep(0,50)
mse8 = rep(0,50)
mse9 = rep(0,50)
for (i in 1:50) {
  mse1[i] = mse_10000_300[[i]]
  mse2[i] = mse_15000_200[[i]]
  mse3[i] = mse_30000_100[[i]]
  mse4[i] = mse_5000_300[[i]]
  mse5[i] = mse_7500_200[[i]]
  mse6[i] = mse_15000_100[[i]]
  mse7[i] = mse_5000_150[[i]]
  mse8[i] = mse_7500_100[[i]]
  mse9[i] = mse_15000_50[[i]]
}

mse_data = data.frame("t"=seq(0,0.98,0.02)*100,mse1,mse2,mse3,mse4,mse5,mse6,mse7,mse8)
mydata <- subset(melt(mse_data,id="t"),value<0.7)
ggplot(data = mydata,aes(x=t,y=value,group = variable,
                         color=variable))+
  geom_line(size=1)+
  scale_colour_discrete(name="Settings",labels=c(
    expression(paste(rho,"=0.067,",tau,"=0.0033",sep = "")),
    expression(paste(rho,"=0.100,",tau,"=0.0050",sep = "")),
    expression(paste(rho,"=0.200,",tau,"=0.0100",sep = "")),
    expression(paste(rho,"=0.033,",tau,"=0.0033",sep = "")),
    expression(paste(rho,"=0.050,",tau,"=0.0050",sep = "")),
    expression(paste(rho,"=0.100,",tau,"=0.0100",sep = "")),
    expression(paste(rho,"=0.033,",tau,"=0.0067",sep = "")),
    expression(paste(rho,"=0.050,",tau,"=0.0100",sep = ""))))+
  scale_x_continuous(name="t",limits = c(0,100),breaks = seq(0,100,10))+
  scale_y_continuous(name="Mean Square Error",limits=c(0,0.7),breaks = seq(0,0.7,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.85,.8),
        legend.box.background = element_rect(color="black"), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))

# 
mse_dependent <- read.csv("~/project_dmc/final_version/datas/dependent_mc.csv")[,2:67]
mse_beta_0 <- mse_dependent[,23]
mse_beta_0.3 <- mse_dependent[,26]
mse_beta_0.6 <- mse_dependent[,29]
mse_beta_0.9 <- mse_dependent[,32]
mse_beta_0_2 <- mse_dependent[,45]
mse_beta_0.3_2 <- mse_dependent[,48]
mse_beta_0.6_2 <- mse_dependent[,51]
mse_beta_0.9_2 <- mse_dependent[,54]
mse_data <- data.frame("t"=seq(0.01,1,0.01)*100,mse_beta_0,mse_beta_0.3,mse_beta_0.6,mse_beta_0.9,
                       mse_beta_0_2,mse_beta_0.3_2,mse_beta_0.6_2,mse_beta_0.9_2)
mydata <- subset(melt(mse_data,id="t"),value<0.7)
ggplot(data = mydata,aes(x=t,y=value,group = variable,
                         color=variable))+
  geom_line(size=1)+
  scale_colour_discrete(name="Settings",
                        labels=c(expression(paste(sigma,"=1,",beta,"=0.0",sep="")),
                                 expression(paste(sigma,"=1,",beta,"=0.3",sep="")),
                                 expression(paste(sigma,"=1,",beta,"=0.6",sep="")),
                                 expression(paste(sigma,"=1,",beta,"=0.9",sep="")),
                                 expression(paste(sigma,"=2,",beta,"=0.0",sep="")),
                                 expression(paste(sigma,"=2,",beta,"=0.3",sep="")),
                                 expression(paste(sigma,"=2,",beta,"=0.6",sep="")),
                                 expression(paste(sigma,"=2,",beta,"=0.9",sep=""))))+
  scale_x_continuous(name="t",limits = c(0,100),breaks = seq(0,100,10))+
  scale_y_continuous(name="Mean Square Error",limits=c(0,0.7),breaks = seq(0,0.7,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.85,.8),
        legend.box.background = element_rect(color="black"), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))

dependent_xi_1 <- rep(0,9)
dependent_xi_2 <- rep(0,9)
dependent_xi_3 <- rep(0,9)
for (i in 1:9) {
  dependent_xi_1[i] <- mean(mse_dependent[30:70,22+i])
  dependent_xi_2[i] <- mean(mse_dependent[30:70,33+i])
  dependent_xi_3[i] <- mean(mse_dependent[30:70,44+i])
}
phi_xi_func <- function(beta){
  ret <- rep(0,length(beta))
  for (i in 0:300) {
    ret <- ret + sqrt(1-sqrt(1-beta^(2*i)))
  }
  return(ret)
}
beta = seq(0,0.8,0.1)
phi_xi <- phi_xi_func(beta)^(1.6)
xy_data <- data.frame(phi=phi_xi,x1=dependent_xi_1,x2=dependent_xi_2,x3=dependent_xi_3)
mydata <- melt(xy_data,id="phi")
model_1 <- lm(x1~phi,xy_data)
model_2 <- lm(x2~phi,xy_data)
model_3 <- lm(x3~phi,xy_data)
ggplot(data = mydata,aes(x=phi,y=value,group = variable,
                         color=variable))+
  geom_point(size=1.5)+
  geom_abline(slope=1.591e-03,intercept = 8.624e-02,linetype=3,color="red")+
  geom_abline(slope= 0.0085866,intercept = 0.0958381,linetype=3,color="green")+
  geom_abline(slope=0.028762,intercept = 0.130730,linetype=3,color="blue")+
  scale_colour_discrete(name="Settings",
                        labels=c(expression(paste(sigma,"=1.0",sep="")),
                                 expression(paste(sigma,"=1.5",sep="")),
                                 expression(paste(sigma,"=2.0",sep=""))))+
  scale_x_continuous(name=TeX("$\\Phi_{Y}^{8/5}$"),limits = c(1,10),breaks = seq(1,10,1))+
  scale_y_continuous(name="Average Mean Square Error",limits=c(0,0.4),breaks = seq(0,0.4,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.95,.9),
        legend.box.background = element_rect(color="black"), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))


# 
mse_X_dependent <- read.csv("~/project_dmc/final_version/datas/dependent_X_mc.csv")[,2:67]
mse_X_beta_0 <- mse_X_dependent[,33]
mse_X_beta_0.3 <- mse_X_dependent[,30]
mse_X_beta_0.6 <- mse_X_dependent[,27]
mse_X_beta_0.9 <- mse_X_dependent[,24]
mse_X_beta_0_2 <- mse_X_dependent[,55]
mse_X_beta_0.3_2 <- mse_X_dependent[,52]
mse_X_beta_0.6_2 <- mse_X_dependent[,49]
mse_X_beta_0.9_2 <- mse_X_dependent[,46]
mse_X_data <- data.frame("t"=seq(0.01,1,0.01)*100,mse_X_beta_0,mse_X_beta_0.3,mse_X_beta_0.6,mse_X_beta_0.9,
                         mse_X_beta_0_2,mse_X_beta_0.3_2,mse_X_beta_0.6_2,mse_X_beta_0.9_2)
mydata <- subset(melt(mse_X_data,id="t"),value<0.7)
ggplot(data = mydata,aes(x=t,y=value,group = variable,
                         color=variable))+
  geom_line(size=1)+
  scale_colour_discrete(name="Settings",
                        labels=c(expression(paste(sigma,"=1,",alpha,"=0.0",sep="")),
                                 expression(paste(sigma,"=1,",alpha,"=0.3",sep="")),
                                 expression(paste(sigma,"=1,",alpha,"=0.6",sep="")),
                                 expression(paste(sigma,"=1,",alpha,"=0.9",sep="")),
                                 expression(paste(sigma,"=2,",alpha,"=0.0",sep="")),
                                 expression(paste(sigma,"=2,",alpha,"=0.3",sep="")),
                                 expression(paste(sigma,"=2,",alpha,"=0.6",sep="")),
                                 expression(paste(sigma,"=2,",alpha,"=0.9",sep=""))))+
  scale_x_continuous(name="t",limits = c(0,100),breaks = seq(0,100,10))+
  scale_y_continuous(name="Mean Square Error",limits=c(0,0.7),breaks = seq(0,0.7,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.85,.8),
        legend.box.background = element_rect(color="black"), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))

dependent_X_1 <- rep(0,9)
dependent_X_2 <- rep(0,9)
dependent_X_3 <- rep(0,9)
for (i in 1:9) {
  dependent_X_1[i] <- mean(mse_X_dependent[30:70,34-i])
  dependent_X_2[i] <- mean(mse_X_dependent[30:70,45-i])
  dependent_X_3[i] <- mean(mse_X_dependent[30:70,56-i])
}
phi_X_func <- function(alpha){
  return(1/(1-alpha^(1/2)))
}
alpha = seq(0,0.8,0.1)
phi_X <- phi_X_func(alpha)^1.6
xy_data <- data.frame(phi=phi_X,x1=dependent_X_1,x2=dependent_X_2,x3=dependent_X_3)
mydata <- melt(xy_data,id="phi")
model_1 <- lm(x1~phi,xy_data)
model_2 <- lm(x2~phi,xy_data)
model_3 <- lm(x3~phi,xy_data)
ggplot(data = mydata,aes(x=phi,y=value,group = variable,
                         color=variable))+
  geom_point(size=1.5)+
  geom_abline(slope=1.281e-03,intercept = 8.550e-02,linetype=3,color="red")+
  geom_abline(slope=0.0015180,intercept = 0.0982554,linetype=3,color="green")+
  geom_abline(slope=1.980e-03,intercept = 1.315e-01,linetype=3,color="blue")+
  scale_colour_discrete(name="Settings",
                        labels=c(expression(paste(sigma,"=1.0",sep="")),
                                 expression(paste(sigma,"=1.5",sep="")),
                                 expression(paste(sigma,"=2.0",sep=""))))+
  scale_x_continuous(name=TeX("$\\Phi_{X}^{8/5}$"),limits = c(1,40),breaks = seq(1,40,4))+
  scale_y_continuous(name="Average Mean Square Error",limits=c(0,0.4),breaks = seq(0,0.4,0.1))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.95,.9),
        legend.box.background = element_rect(color="black"),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))

#
mse_precise <- read.csv("~/project_dmc/final_version/datas/phase_transition_precise.csv")[,2:27]
plot_result <- rep(0,26)
for (i in 1:26) {
  plot_result[i] <- mean(mse_precise[30:70,i])
}
index <- seq(1/3,2,1/15)*100/10
mydata = data.frame(index=log(index),mse=log(plot_result))
model <- lm(mse~index,data = mydata[1:11,])
ggplot(data = mydata[1:11,],aes(x=index,y=mse))+
  geom_point(size=1.5)+
  geom_abline(slope=-0.82587,intercept = -0.34902,linetype=3)+
  scale_x_continuous(name=TeX("$\\log \\frac{\\rho}{\\tau}$"),limits = c(1,2.5),breaks=seq(1,2.5,0.5))+
  scale_y_continuous(name="Log of Average Mean Square Error",limits=c(-2.5,-1),breaks = seq(-2.5,-1,0.5))+
  theme( panel.grid.minor = element_blank(),legend.position = c(.85,.9),
         legend.box.background = element_rect(color="black"), 
         axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))


# 
origin = "1999-12-30"
finish = "2005-12-18"
finish_ = "2006-3-01"
T_=100
mse_total = matrix(0,T_,101)
for (i in 1:T_) {
  mse_total[i,] <- read.csv(paste("~/project_dmc/final_version/datas/mse_",i,".csv",sep = ""))[1:101,2]
}
mse_baseline = matrix(0,T_,101)
for (i in 1:T_) {
  mse_baseline[i,] <- read.csv(paste("~/project_dmc/final_version/datas/baseline_mse_",i,".csv",sep = ""))[1:101,2]
}
mse_netflix = mse_total[,73]
mse_netflix_baseline = mse_baseline[,73]
mse_netflix_twostep = read.csv("~/project_dmc/final_version/datas/baseline_mse_twostep")[,2]
mse_netflix_tensor = read.csv("~/project_dmc/final_version/datas/baseline_mse_tensor")[,2]
mse_data <- data.frame("t"=seq.Date(as.Date(origin),as.Date(finish),by=22),"Static"=mse_netflix_baseline,"TwoStep"=mse_netflix_twostep,
                       "Tensor"=mse_netflix_tensor,"DLRTR"=mse_netflix)
mydata <- melt(mse_data,id="t")
ggplot(data = mydata,aes(x=t,y=value,group=variable,
                         color=variable))+
  geom_line(size=1)+
  scale_colour_discrete(name="Estimator",
                        labels=c("Static","TwoStep","Tensor","DLR"))+
  scale_x_continuous(name="Date",limits = c(as.Date(origin),as.Date(finish_)),breaks =seq.Date(as.Date(origin),as.Date(finish_),by=172))+
  scale_y_continuous(name="Mean Square Error",limits=c(0,1.5),breaks = seq(0,1.5,0.3))+
  theme(panel.grid.minor = element_blank(),legend.position = c(.85,.25),
        legend.box.background = element_rect(color="black"),
        axis.text.x = element_text(angle = 70, hjust = 0, vjust = 0,size=16),
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18),
        legend.title =element_text(size=18),
        legend.text = element_text(size=18))



