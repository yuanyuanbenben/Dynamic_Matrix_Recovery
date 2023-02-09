library(jpeg)
library(ggplot2)
library(reshape2)

# read pic
img_total_r <- array(0,dim = c(96,480,854))
img_total_g <- array(0,dim = c(96,480,854))
img_total_b <- array(0,dim = c(96,480,854))
for (index in 0:95) {
  if (index < 10){
    dic <- paste("~/Dynamic_Matrix_Recovery/data/lions/0000",index,".jpg",sep = "")
  }
  else{
    dic <- paste("~/Dynamic_Matrix_Recovery/data/lions/000",index,".jpg",sep = "")
  }
  img <- readJPEG(dic)
  img_total_r[index+1,,] <- img[,,1]
  img_total_g[index+1,,] <- img[,,2]
  img_total_b[index+1,,] <- img[,,3]
}

# show jpeg
# test 
# par(mar=c(0,1,0,1))
# img <- array(0,dim = c(480,854,3))
# img[,,1] <- twostep_red[1,,]#img_total_r[85,,]
# img[,,2] <- twostep_green[1,,]#img_total_g[85,,]
# img[,,3] <- twostep_blue[1,,]#img_total_b[85,,]
# longImage<-melt(img)
# rgbImage<-reshape(longImage,timevar='Var3',idvar=c('Var1','Var2'),direction='wide')
# rgbImage$Var1<- -rgbImage$Var1
# colorColumns<- rgbImage[, substr(colnames(rgbImage), 1, 5)== "value"]
# with(rgbImage,plot(Var2, Var1, col = rgb(colorColumns), asp = 1, pch =".",axes=F,xlab='',ylab=''))

save.image("~/Dynamic_Matrix_Recovery/code/real_data/vedio_data.RData")
