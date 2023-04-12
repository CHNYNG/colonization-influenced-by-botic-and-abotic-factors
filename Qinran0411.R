
##### set the local folder
setwd("E:/黑石顶测菌根/菌根侵染率/数据整理")

###### input the data
Qrl <- read.csv("菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
N_data <- read.csv("N整理.csv",header = T,fileEncoding = "GBK")
S_data <- read.csv("S整理.csv",header = T,fileEncoding = "GBK")
weigh <- read.csv("weigh.csv",header = T,fileEncoding = "GBK")
Nor_data <- read.csv("整理.csv",header = T,fileEncoding = "GBK")
Qrl <- read.csv("菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
HSD_data <- read.csv("E:/Chu Lab/hsd.species.alive.cy.branch0.csv",header = T,fileEncoding = "GBK")
HSD_env <- read.csv("E:/Chu Lab/HSD.origin/soil_origin.csv",header = T,fileEncoding = "GBK")


###### select the useful data of data morphology
# insert diameters
library(stringr)
diameter <- read.table("E:/黑石顶测菌根/菌根侵染率/数据整理/ylj_for_diameter.txt",
                       header = F,sep="\t",na.strings = "NA",
                       fill = T,fileEncoding = "GBK")
diameter <- diameter[,-c(3:17,19,21,23:109)]
colnames(diameter) <- c("TagNew","Opterator","ProjArea(cm2)","SurfArea(cm2)","AvgDiam(mm)")
diameter$TagNew <- as.character(diameter$TagNew)
str(diameter$TagNew)
diameter$TagNew = str_pad(diameter$TagNew,7,side = "left", "0")

#insert length
length <- read.table("E:/黑石顶测菌根/菌根侵染率/数据整理/ylj_for_length.txt",
header = F,sep="\t",na.strings = "NA",
fill = T,fileEncoding = "GBK")
length <- length[-c(1:5),-c(3:15,17:109)]
colnames(length) <- c("TagNew","Opterator","Length(cm)")
length$TagNew <- as.character(length$TagNew)
str(length$TagNew)
length$TagNew = str_pad(length$TagNew,7,side = "left", "0")


##### merge the root morphology data
library(dplyr)
root_morphology <- left_join(length,diameter,by = "TagNew")
weigh$TagNew <- as.character(weigh$TagNew)
str(weigh$TagNew)
weigh$TagNew = str_pad(weigh$TagNew,7,side = "left", "0")
root_morphology <- left_join(root_morphology,weigh,by = "TagNew")
#root_warning <- subset(root_morphology, is.na(Weight.g))
#View(root_warning)
#write.csv(root_warning,"root_warning.csv",fileEncoding = "GBK")

root_morphology$'SRL(gcm)' <- as.numeric(root_morphology$Weight.g)/
  as.numeric(root_morphology$`Length(cm)`)

###calculate the colonization
gcsys <- aggregate(Qrl$观察视野总数, by = list(type = Qrl$Numbers), sum)
EMsys <- aggregate(Qrl$ECM, by = list(type = Qrl$Numbers), sum)
jssys <- aggregate(Qrl$菌丝出现视野数, by = list(type = Qrl$Numbers), sum)
BZ <- aggregate(Qrl$孢子, by = list(type = Qrl$Numbers), sum)
Hbbz <- aggregate(Qrl$厚壁孢子, by = list(type = Qrl$Numbers), sum)
Pn <- aggregate(Qrl$泡囊, by = list(type = Qrl$Numbers), sum)
Pq <- aggregate(Qrl$盘曲结构, by = list(type = Qrl$Numbers), sum)
Fzxb <- aggregate(Qrl$辅助细胞群数, by = list(type = Qrl$Numbers), sum)
Qrl_ca <- data.frame()
Qrl_ca <- cbind(gcsys$type,gcsys$x,EMsys$x,jssys$x,BZ$x,Hbbz$x,Pn$x,Pq$x,Fzxb$x)
colnames(Qrl_ca) <- c("Numbers","gcsys","EMsys","jssys","BZ","Hbbz","Pn","Pq","Fzxb")
Qrl_ca <- as.data.frame(Qrl_ca)
Qrl_ca[,2:9] <- as.numeric(unlist(Qrl_ca[,2:9]))
Qrl_ca$qr_AM <- Qrl_ca$jssys/Qrl_ca$gcsys
Qrl_ca$qr_EM <- Qrl_ca$EM/Qrl_ca$gcsys
Qrl_ca$qr_BZ <- Qrl_ca$BZ/Qrl_ca$gcsys
Qrl_ca$qr_Hbbz <- Qrl_ca$Hbbz/Qrl_ca$gcsys
Qrl_ca$qr_Pn <- Qrl_ca$Pn/Qrl_ca$gcsys
Qrl_ca$qr_Pq <- Qrl_ca$Pq/Qrl_ca$gcsys
Qrl_ca$qr_fzxb <- Qrl_ca$Fzxb/Qrl_ca$gcsys

#####merge the colonization data to TagNew
total_data <- rbind(Nor_data,S_data,N_data)
colnames(total_data) <- c("TagNew","Numbers")
total_data$TagNew <- as.character(total_data$TagNew)
str(total_data$TagNew)
total_data$TagNew = str_pad(total_data$TagNew,7,side = "left", "0")

root_qrl <- data.frame()
root_qrl <- left_join(Qrl_ca,total_data,by="Numbers")

HSD_data$TagNew <- as.character(HSD_data$TagNew)
str(HSD_data$TagNew)
HSD_data$TagNew = str_pad(HSD_data$TagNew,7,side = "left", "0")
root_qrl <- left_join(root_qrl,HSD_data,by="TagNew")
root_qrl <- left_join(root_qrl,root_morphology,by="TagNew")


#####
#Then the krige
####马珂 黑石顶数据插值
#加载这个包
library(geoR)
###将root_qrl里的NA筛掉，填进root_qrl_env里
root_qrl_env <- root_qrl %>% filter(!is.na(GX) | !is.na(GY))
# 需要插值的坐标
dat <- data.frame(x = root_qrl_env$GX, y = root_qrl_env$GY)


for (ii in 1:31){
###在原始的数据里，选项一列数据
elev.g <- as.geodata(HSD_env, coords.col = 2:3, data.col = ii+3)
###科技黑箱，不动就行
elev.variog <- variog(elev.g,uvec = seq(0, 500, by = 5), max.dist = 500) # trend="2nd"可以不要
elev.sph <- variofit(elev.variog, cov.model = "spherical") # 简单的选择了球面模型，模型选择要拟合，比较麻烦我一般不做


# 最后插值
elev.prd <- krige.conv(elev.g, loc = dat,
                      krige = krige.control(
                        cov.model = "spherical",
                        cov.pars = c(elev.sph$cov.pars[1],
                                     elev.sph$cov.pars[2])))

# 结果提取：
root_qrl_env <- cbind(root_qrl_env, elev.prd$predict)
colnames(root_qrl_env)[ii+55] <- paste0(colnames(HSD_env)[ii+3])
}

#write.csv(root_qrl_env,"BS_root_qrl_env.csv",fileEncoding = "GBK")
