######
##这部分用来整理数据##
###整理了所选所有个体的AM侵染率、EcM侵染率，个体信息
#（包含TagNew、个体DBH和物种中文名、拉丁名和科属种）
###根系信息（包含根系表面积、平均根直径、比根长

###### input the data
Qrl <- read.csv("data/菌根侵染率_整合.csv",header = T,fileEncoding = "GBK")
N_data <- read.csv("data/N整理.csv",header = T,fileEncoding = "GBK",stringsAsFactors = FALSE)
S_data <- read.csv("data/S整理.csv",header = T,fileEncoding = "GBK",stringsAsFactors = FALSE)
weigh <- read.csv("data/weigh.csv",header = T,fileEncoding = "GBK")
Nor_data <- read.csv("data/整理.csv",header = T,fileEncoding = "GBK",stringsAsFactors = FALSE)
HSD_data <- read.csv("data/HSD.csv",header = T,fileEncoding = "GBK")


###### select the useful data of data morphology
# insert diameters
library(stringr)
library(dplyr)
diameter <- read.table("data/ylj_for_diameter.txt",
                       header = F,sep="\t",na.strings = "NA",
                       fill = T,fileEncoding = "GBK")
diameter <- diameter%>%
  select(-c(2:17, 19, 21, 23:103)) %>%
  rename(TagNew = V1, ProjArea = V18, SurfArea = V20, AvgDiam = V22) %>%
  mutate(TagNew = str_pad(as.character(TagNew), width = 7, side = "left", pad = "0"))

HSD_data <- HSD_data %>%
  mutate(TagNew = str_pad(TagNew, width = 7, side = "left", pad = "0")) %>%
  filter(!(Latin %in% c("A-VINE", "unknown", "T_VINE")))

#insert length
length <- read.table("data/ylj_for_length.txt",
header = F,sep="\t",na.strings = "NA",
fill = T,fileEncoding = "GBK")

length <- length[-c(1:5), -c(3:15, 17:109)] %>%
  rename(TagNew = V1, Opterator = V2, Length = V16) %>%
  mutate(TagNew = str_pad(as.character(TagNew), width = 7, side = "left", pad = "0"))
#write.csv(length,"E:/黑石顶测菌根/菌根侵染率/数据整理/length_correct.csv")
#write.csv(diameter,"E:/黑石顶测菌根/菌根侵染率/数据整理/diameter_correct.csv")

##### merge the root morphology data
weigh <- weigh %>%
  mutate(TagNew = str_pad(as.character(TagNew), width = 7, side = "left", pad = "0"))

root_morphology <- length %>%
  left_join(diameter, by = "TagNew") %>%
  left_join(weigh, by = "TagNew")
#root_warning <- subset(root_morphology, is.na(Weight.g))
#View(root_warning)
#write.csv(root_warning,"root_warning.csv",fileEncoding = "GBK")

root_morphology <- root_morphology %>%
  mutate(SRL = as.numeric(Length) / as.numeric(Weight.g),
         SRA = as.numeric(SurfArea) / as.numeric(Weight.g))


Qrl <- Qrl %>%
  mutate(AM_range = 菌丝出现视野数 / 观察视野总数,
         EM_range = ECM / 观察视野总数,
         am = case_when(
           AM_range <= 0.25 ~ 0.125,
           AM_range > 0.25 & AM_range <= 0.5 ~ 0.375,
           AM_range > 0.5 & AM_range <= 0.75 ~ 0.625,
           AM_range > 0.75 ~ 0.875
         ),
         em = case_when(
           EM_range <= 0.25 ~ 0.125,
           EM_range > 0.25 & EM_range <= 0.5 ~ 0.375,
           EM_range > 0.5 & EM_range <= 0.75 ~ 0.625,
           EM_range > 0.75 ~ 0.875
         ))


###calculate the colonization
Qrl_ca <- Qrl %>%
  group_by(Numbers) %>%
  summarize(
    gcsys = sum(观察视野总数),
    EMsys = sum(ECM),
    jssys = sum(菌丝出现视野数),
    BZ = sum(孢子),
    Hbbz = sum(厚壁孢子),
    Pn = sum(泡囊),
    Pq = sum(盘曲结构),
    Fzxb = sum(辅助细胞群数),
    am = mean(am),
    em = mean(em)
  ) %>%
  mutate(
    qr_AM = jssys / gcsys,
    qr_EM = EMsys / gcsys,
    qr_BZ = BZ / gcsys,
    qr_Hbbz = Hbbz / gcsys,
    qr_Pn = Pn / gcsys,
    qr_Pq = Pq / gcsys,
    qr_fzxb = Fzxb / gcsys
  ) %>%
  select(qr_AM, qr_EM, qr_BZ, qr_Hbbz, qr_Pn, qr_Pq, qr_fzxb, Numbers,am, em)


#####merge the colonization data to TagNew
Nor_data <- Nor_data %>%
  mutate(Numbers = as.character(Numbers))

S_data <- S_data %>%
  mutate(Numbers = as.character(Numbers))

N_data <- N_data %>%
  mutate(Numbers = as.character(Numbers))

total_data <- bind_rows(Nor_data, S_data, N_data) %>%
  mutate(TagNew = str_pad(as.character(TagNew), width = 7, side = "left", pad = "0"))

root_qrl <- data.frame()

root_qrl <- Qrl_ca %>%
  left_join(total_data, by = "Numbers") %>%
  left_join(HSD_data, by = "TagNew", relationship = "many-to-many") %>%
  left_join(root_morphology, by = "TagNew")

#整理下root_qrl
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

# 创建 M_Type 列并初始化为 "NA"
M_Type_deciede$M_Type <- "NA"

# 根据条件更新 M_Type 列的值
M_Type_deciede <- M_Type_deciede %>%
  mutate(M_Type = ifelse(qr_AM > 0, "AM", M_Type),
         M_Type = ifelse(qr_EM > 0, "EM", M_Type),
         M_Type = ifelse(qr_AM > 0 & qr_EM > 0, "AM-EM", M_Type))

###把物种的数据和单个个体的数据合并
M_Type_deciede <- M_Type_deciede[, c("Latin", "M_Type")]

root_qrl<- left_join(root_qrl, M_Type_deciede, by = "Latin")

# 试试看如果按照个体
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
#这几个包服务器都加载不了，所以用电脑把数据抠出来，在存进来导入
#install_github("helixcn/plantlist", build_vignettes = TRUE)
#加载plantlist包
library(plantlist)
#先筛选出所有的物种
# 获取不重复的拉丁名称列
unique_species <- unique(HSD_data$Latin)

specieslist <- TPL(unique_species)
colnames(specieslist) <- c("Latin","Genus","Family","Family_number","Order","Group")

#把正名后的物种加入到大数据里

#读进来物种数据
#这个unique_species之后还会用到
#unique_species <- unique(HSD_data$Latin)
#load(file =  "data/specieslist.RData")
root_qrl <- left_join(root_qrl,specieslist[,c("Latin","Genus","Family","Order")],
                      by = "Latin")
###
#先把branch=0筛出来？不需要那么多
root_qrl <- root_qrl %>%
  filter(Branch ==0) %>%
  select(-Numbers) %>%
  distinct(TagNew, .keep_all = TRUE) %>%
  mutate(Latin = gsub(" ", "_", Latin))


#给物种信息合并菌根类型
###先不加，嘤嘤嘤-


