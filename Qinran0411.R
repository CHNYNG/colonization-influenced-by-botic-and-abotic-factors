Qrl <- read.csv("菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
###set the local folder
setwd("E:/黑石顶测菌根/菌根侵染率/数据整理")
Qrl <- read.csv("菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
N_data <- read.csv("N整理.csv",header = T,fileEncoding = "GBK")
S_data <- read.csv("S整理.csv",header = T,fileEncoding = "GBK")
weigh <- read.csv("weigh.csv",header = T,fileEncoding = "GBK")
Nor_data <- read.csv("整理.csv",header = T,fileEncoding = "GBK")
Qrl <- read.csv("菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
HSD_data <- read.csv("E:/Chu Lab/HSD.origin/HSD.csv",header = T,fileEncoding = "GBK")
HSD_env <- read.csv("E:/Chu Lab/HSD.origin/soil_origin.csv",header = T,fileEncoding = "GBK")
###select the useful data of data morphology
library(stringr)
diameter <- read.table("E:/黑石顶测菌根/菌根侵染率/数据整理/ylj_for_diameter.txt",
header = F,sep="\t",na.strings = "NA",
fill = T,fileEncoding = "GBK")
diameter <- diameter[,-c(3:17,19,21,23:109)]
colnames(diameter) <- c("TagNew","Opterator","ProjArea(cm2)","SurfArea(cm2)","AvgDiam(mm)")
diameter$TagNew <- as.character(diameter$TagNew)
str(diameter$TagNew)
diameter$TagNew = str_pad(diameter$TagNew,7,side = "left", "0")
length <- read.table("E:/黑石顶测菌根/菌根侵染率/数据整理/ylj_for_length.txt",
header = F,sep="\t",na.strings = "NA",
fill = T,fileEncoding = "GBK")
length <- length[-c(1:5),-c(3:15,17:109)]
colnames(length) <- c("TagNew","Opterator","Length(cm)")
length$TagNew <- as.character(length$TagNew)
str(length$TagNew)
length$TagNew = str_pad(length$TagNew,7,side = "left", "0")
###merge the root morphology data
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
#####merge the colonization data to TagNew
###calculate the colonization
root_qrl <- left_join(Qrl,Nor_data,by="Numbers")
#####merge the colonization data to TagNew
###calculate the colonization
Nor_data$Numbers <- as.character(Nor_data$Numbers)
root_qrl <- left_join(Qrl,Nor_data,by="Numbers")
View(root_qrl)
root_qrl <- left_join(Qrl,S_data,by="Numbers")
View(root_qrl)
root_qrl <- left_join(Qrl, N_data, by = "Numbers",relationship = "many-to-many")
rm(root_qrl)
#####merge the colonization data to TagNew
###calculate the colonization
Nor_data$Numbers <- as.character(Nor_data$Numbers)
root_qrl <- left_join(Qrl,Nor_data,by="Numbers")
View(root_qrl)
root_qrl <- left_join(root_qrl, S_data, by = "Numbers",relationship = "many-to-many")
root_qrl <- left_join(root_qrl, N_data, by = "Numbers",relationship = "many-to-many")
root_qrl$TagNew <- paste(root_qrl$TagNew,root_qrl$TagNew.y,sep = "")
root_qrl$TagNew <- paste(root_qrl$TagNew,root_qrl$TagNew.x,sep = "")
rm(root_qrl)
root_qrl <- left_join(Qrl,Nor_data,by="Numbers")
root_qrl <- left_join(root_qrl, S_data, by = "Numbers",relationship = "many-to-many")
root_qrl <- left_join(root_qrl, N_data, by = "Numbers",relationship = "many-to-many")
root_qrl$TagNew <- paste(root_qrl$TagNew,root_qrl$TagNew.x,sep = "")
root_qrl$TagNew <- paste(root_qrl$TagNew,root_qrl$TagNew.y,sep = "")
View(root_qrl)
