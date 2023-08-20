####导入刘老师的数据
Soil <- read.csv("data/HSD_Soil_Nutrition_Liu.csv",header = T,fileEncoding = "GBK")
####把数据和坐标合并一下
soil_to_liu <- read.csv("data/Ru Wang/HSD_Soil_To Liu.csv",header = T,fileEncoding = "GBK")
Soil <- left_join(Soil,soil_to_liu,by = "SampleID")
Soil <- na.omit(Soil)

###画一下
par(mfrow = c(2,4))
apply(Soil[,2:9], 2, hist)
##都挺正态的
dev.off()

library(automap)
library(sp)

root_qrl_na <- root_qrl[!is.na(root_qrl$GX), ]
prd.loc <- root_qrl_na[,c("GX","GY","TagNew")]
prd.loc$GX <- as.numeric(prd.loc$GX)
prd.loc$GY <- as.numeric(prd.loc$GY)
coordinates(Soil) = ~ GX + GY
coordinates(prd.loc)= ~ GX + GY
dat.krnw <- root_qrl_na
c <-  automap::autoKrige(pH ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "pH"
c <-  automap::autoKrige(AP..mg.kg. ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AP..mg.kg."
c <-  automap::autoKrige(N.P ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "N.P"
c <-  automap::autoKrige(C.P ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "C.P"
c <-  automap::autoKrige(C.N ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "C.N"
c <-  automap::autoKrige(TP..mg.g. ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TP..mg.g."
c <-  automap::autoKrige(TN..mg.g. ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TN..mg.g."
c <-  automap::autoKrige(SOC..mg.g. ~ 1, Soil, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "SOC..mg.g."
library(stringr)
dat.krnw$TagNew = str_pad(dat.krnw$TagNew,7,side = "left", "0")
root_qrl_soil <- left_join(root_qrl,dat.krnw[
  ,c("SOC..mg.g.","TN..mg.g.","TP..mg.g.","C.N","C.P","N.P","AP..mg.kg.","pH" ,"TagNew")],
  by="TagNew")


###To weitao
prd.loc <- as.data.frame(prd.loc)
Soil <- as.data.frame(Soil)
save(prd.loc,Soil,file = "data/Weitao/soildata.RData")
save(dat.krnw,file = "data/Weitao/soildata_vyang.RData")
