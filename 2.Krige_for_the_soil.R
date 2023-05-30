##再用回之前的代码
###插值土壤元素
####导入文启师兄的数据（删掉了单位的一行）
soil <- read.csv("data/soil_origin.csv",header = T)
###画一下
#x11()
#par(mfrow = c(5,7))
#apply(soil[,4:34], 2, hist)
#dev.off()
#正态一下
a <- apply(soil[,c("SOM","TP","EAl","TZn","TFe","TNa","AK","EK","ETB","TCu","AMn","ENa","ECa","AS","TMn","TCa","EC")], 2, log)
a <- cbind(soil[,c("X","gx","gy","pH","AN","AP","TN","TK","Sand","Silt","ACu","AZn","AFe",
                   "EMg","AB","TMg","MWC")],a)
#再画一下
#par(mfrow = c(5,7))
#apply(a[,4:34], 2, hist)
###auto克里金
library(automap)
library(sp)
rm(b)
b <- a
load("E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl.RData")
root_qrl_na <- root_qrl[!is.na(root_qrl$GX), ]
prd.loc <- root_qrl_na[,c("GX","GY")]
colnames(prd.loc) <- c("gx","gy") 
coordinates(b) = ~ gx + gy
coordinates(prd.loc)= ~ gx + gy
rm(dat.krnw)
dat.krnw <- root_qrl_na
c <-  automap::autoKrige(pH ~ 1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "pH"
c <-  automap::autoKrige(EC ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "EC"
c <-  automap::autoKrige(SOM ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "SOM"
c <-  automap::autoKrige(AN ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AN"
c <-  automap::autoKrige(AP ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AP"
c <-  automap::autoKrige(AK ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AK"
c <-  automap::autoKrige(TN ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TN"
c <-  automap::autoKrige(TP ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TP"
c <-  automap::autoKrige(TK ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TK"
c <-  automap::autoKrige(Sand ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "Sand"
c <-  automap::autoKrige(Silt ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "Silt"
c <-  automap::autoKrige(ACu ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "ACu"
c <-  automap::autoKrige(AZn ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AZn"
c <-  automap::autoKrige(AFe ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AFe"
c <-  automap::autoKrige(AMn ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AMn"
c <-  automap::autoKrige(EK ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "EK"
c <-  automap::autoKrige(ENa ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "ENa"
c <-  automap::autoKrige(ECa ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "ECa"
c <-  automap::autoKrige(EMg ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "EMg"
c <-  automap::autoKrige(ETB ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "ETB"
c <-  automap::autoKrige(AB ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AB"
c <-  automap::autoKrige(AS ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "AS"
c <-  automap::autoKrige(EAl ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "EAl"
c <-  automap::autoKrige(TCu ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TCu"
c <-  automap::autoKrige(TZn ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TZn"
c <-  automap::autoKrige(TFe ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TFe"
c <-  automap::autoKrige(TMn ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TMn"
c <-  automap::autoKrige(TNa ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TNa"
c <-  automap::autoKrige(TCa ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TCa"
c <-  automap::autoKrige(TMg ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "TMg"
c <-  automap::autoKrige(MWC ~1, b, prd.loc)
dat.krnw <- cbind(c$krige_output$var1.pred, dat.krnw)
colnames(dat.krnw)[1] <- "MWC"

###搞定了插值，再正态回来
d <- apply(dat.krnw[,c("SOM","TP","EAl","TZn","TFe","TNa","AK","EK","ETB","TCu",
                       "AMn","ENa","ECa","AS","TMn","TCa","EC")],2,exp)
d <- cbind(root_qrl_na, d, dat.krnw[,c("pH","AN","AP","TN","TK","Sand"
                                    ,"Silt","ACu","AZn","AFe","EMg","AB","TMg","MWC")])
##插值搞定！
###pca一下，康康有木有问题
library(ggplot2)
#scale一下数据
e <- apply(d[,c("SOM","TP","EAl","TZn","TFe","TNa","AK","EK","ETB","TCu",
                "AMn","ENa","ECa","AS","TMn","TCa","EC","pH","AN","AP","TN","TK","Sand"
                ,"Silt","ACu","AZn","AFe","EMg","AB","TMg","MWC")],2,scale)
e <- cbind(e,d)
###开始pca
pz <- prcomp(e[,1:31])
summary(pz)
pca <- as.data.frame(pz$x)
#ggbiplot(pz, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE) +
#  scale_color_discrete(name = '') +
#  theme(legend.direction = 'horizontal', legend.position = 'top')

pca$Sp_Name <- e$Family
b <- ggplot(pca, aes(x=PC1, y=PC2,colour=Sp_Name),
            fill = Sp_Name)+geom_point()

b

####保存一下RData
#save(d,file = "E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl_soil.RData")
####整理一下d
###删除一些重复的列，整理成它纯粹的亚子