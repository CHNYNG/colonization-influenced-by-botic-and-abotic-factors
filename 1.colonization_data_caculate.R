######
##这部分用来整理数据##
###整理了所选所有个体的AM侵染率、EcM侵染率，个体信息
#（包含TagNew、个体DBH和物种中文名、拉丁名和科属种）
###根系信息（包含根系表面积、平均根直径、比根长

###### input the data
Qrl <- read.csv("data/菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
N_data <- read.csv("data/N整理.csv",header = T,fileEncoding = "GBK")
S_data <- read.csv("data/S整理.csv",header = T,fileEncoding = "GBK")
weigh <- read.csv("data/weigh.csv",header = T,fileEncoding = "GBK")
Nor_data <- read.csv("data/整理.csv",header = T,fileEncoding = "GBK")
Qrl <- read.csv("data/菌根侵染率_观察记录_整合.csv",header = T,fileEncoding = "GBK")
HSD_data <- read.csv("data/HSD.csv",header = T,fileEncoding = "GBK")


###### select the useful data of data morphology
# insert diameters
library(stringr)
library(dplyr)
diameter <- read.table("data/ylj_for_diameter.txt",
                       header = F,sep="\t",na.strings = "NA",
                       fill = T,fileEncoding = "GBK")
diameter <- diameter[,-c(3:17,19,21,23:109)]
colnames(diameter) <- c("TagNew","Opterator","ProjArea(cm2)","SurfArea(cm2)","AvgDiam(mm)")
diameter$TagNew <- as.character(diameter$TagNew)
str(diameter$TagNew)
diameter$TagNew = str_pad(diameter$TagNew,7,side = "left", "0")
HSD_data$TagNew = str_pad(HSD_data$TagNew,7,side = "left", "0")
HSD_data <- HSD_data %>%
  filter(!(Latin %in% c("A-VINE", "unknown", "T_VINE")))

#insert length
length <- read.table("data/ylj_for_length.txt",
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

#这里做完了，可以把已经合并的数据集删掉了
rm(gcsys,EMsys,jssys,BZ,Hbbz,Pn,Pq,Fzxb)


#####merge the colonization data to TagNew
total_data <- rbind(Nor_data,S_data,N_data)
colnames(total_data) <- c("TagNew","Numbers")
total_data$TagNew <- as.character(total_data$TagNew)
str(total_data$TagNew)
total_data$TagNew = str_pad(total_data$TagNew,7,side = "left", "0")

root_qrl <- data.frame()
root_qrl <- left_join(Qrl_ca,total_data,by="Numbers")



root_qrl <- left_join(root_qrl,HSD_data,by="TagNew",relationship = "many-to-many")
root_qrl <- left_join(root_qrl,root_morphology,by="TagNew")

#整理下root_qrl
root_qrl <- as.data.frame(root_qrl)
root_qrl <- root_qrl %>% select(-X, -Species.y)

#### 这里处理完了，把已经有的数据删一删
rm(diameter,weigh,total_data,S_data,N_data,length,Nor_data,Qrl,Qrl_ca,root_morphology)
#save(root_qrl,file = "E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl.RData")

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
###这里用完了species，把它删掉吧
rm(species,total_qr_AM,total_qr_EM,new_row,subset_data)

###给它们赋予菌根类型
library(dplyr)

# 创建 M_Type 列并初始化为 "NA"
M_Type_deciede$M_Type <- "NA"

# 根据条件更新 M_Type 列的值
M_Type_deciede <- M_Type_deciede %>%
  mutate(M_Type = ifelse(qr_AM > 0, "AM", M_Type),
         M_Type = ifelse(qr_EM > 0, "EM", M_Type),
         M_Type = ifelse(qr_AM > 0 & qr_EM > 0, "AM-EM", M_Type))

###把物种的数据和单个个体的数据合并
library(dplyr)
M_Type_deciede <- M_Type_deciede[, c("Latin", "M_Type")]

root_qrl<- left_join(root_qrl, M_Type_deciede, by = "Latin")

# 试试看如果按照个体
library(dplyr)

root_qrl <- root_qrl %>%
  mutate(Type = case_when(
    !is.na(qr_AM) & qr_AM > 0 & !is.na(qr_EM) & qr_EM > 0 ~ "AM-EM",
    !is.na(qr_AM) & qr_AM > 0 ~ "AM",
    !is.na(qr_EM) & qr_EM > 0 ~ "EM",
    TRUE ~ "NA"  # 二者均为 NA 或其他情况
  ))



###这里用species来merge一下数据库里的菌根类型
##### 给数据库加载物种信息
#library(devtools)
#install_github("helixcn/plantlist", build_vignettes = TRUE)
#加载plantlist包
library(plantlist)
#先筛选出所有的物种
# 获取不重复的拉丁名称列
unique_species <- unique(HSD_data$Latin)
specieslist <- TPL(unique_species)
colnames(specieslist) <- c("Latin","Genus","Family","Family_number","Order","Group")
#把正名后的物种加入到大数据里
root_qrl <- left_join(root_qrl,specieslist[,c("Latin","Genus","Family","Order")],
                      by = "Latin")
#给物种信息合并菌根类型
###先不加，嘤嘤嘤-
