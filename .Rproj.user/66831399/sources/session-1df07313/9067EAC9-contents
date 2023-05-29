
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
#write.csv(length,"E:/黑石顶测菌根/菌根侵染率/数据整理/length_correct.csv")
#write.csv(diameter,"E:/黑石顶测菌根/菌根侵染率/数据整理/diameter_correct.csv")

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

HSD_data_0 <- subset(HSD_data,Branch==0)
HSD_data_0$TagNew <- as.character(HSD_data_0$TagNew)
str(HSD_data_0$TagNew)
HSD_data_0$TagNew = str_pad(HSD_data_0$TagNew,7,side = "left", "0")
root_qrl <- left_join(root_qrl,HSD_data_0,by="TagNew")
root_qrl <- left_join(root_qrl,root_morphology,by="TagNew")

#整理下root_qrl
root_qrl <- as.data.frame(root_qrl)
root_qrl <- root_qrl %>% select(-X.x, -Species.y, -Species.y.y)

save(root_qrl,file = "E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl.RData")

#### 获取物种的菌根类型
# 创建一个新的数据框"M_Type_deciede"
M_Type_deciede <- data.frame()

# 获取唯一的Latin物种
unique_species <- unique(root_qrl$Latin)

# 循环遍历每个唯一的Latin物种
for (species in unique_species) {
  # 选择具有相同Latin物种的行
  subset_data <- root_qrl[root_qrl$Latin == species, ]
  
  # 计算qr_AM和qr_EM的总和（忽略NA值）
  total_qr_AM <- sum(subset_data$qr_AM, na.rm = TRUE)
  total_qr_EM <- sum(subset_data$qr_EM, na.rm = TRUE)
  
  # 创建一个新的行，将Latin、qr_AM和qr_EM添加到"M_Type_deciede"中
  new_row <- data.frame(Latin = species, qr_AM = total_qr_AM, qr_EM = total_qr_EM)
  
  # 将新行添加到"M_Type_deciede"
  M_Type_deciede <- rbind(M_Type_deciede, new_row)
}

# 查看"M_Type_deciede"数据框
M_Type_deciede

###给它们赋予菌根类型
library(dplyr)

# 创建 M_Type 列并初始化为 "NA"
M_Type_deciede$M_Type <- "NA"

# 根据条件更新 M_Type 列的值
M_Type_deciede <- M_Type_deciede %>%
  mutate(M_Type = ifelse(qr_AM > 0, "AM", M_Type),
         M_Type = ifelse(qr_EM > 0, "EM", M_Type),
         M_Type = ifelse(qr_AM > 0 & qr_EM > 0, "AM-EM", M_Type))

# 试试看如果按照个体
library(dplyr)

root_qrl <- root_qrl %>%
  mutate(M_Type = case_when(
    !is.na(qr_AM) & qr_AM > 0 & !is.na(qr_EM) & qr_EM > 0 ~ "AM-EM",
    !is.na(qr_AM) & qr_AM > 0 ~ "AM",
    !is.na(qr_EM) & qr_EM > 0 ~ "EM",
    TRUE ~ "NA"  # 二者均为 NA 或其他情况
  ))

###没有了
