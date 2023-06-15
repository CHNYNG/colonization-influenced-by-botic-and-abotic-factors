###一个题外工作——to 王汝
####确定一下需要的数据

####一个是可以给她菌根侵染率&根系数据
####另一个是可以给她新测的土壤数据


#导入几个数据
soil_to_liu <- read.csv("data/Ru Wang/HSD_Soil_To Liu.csv",header = T,fileEncoding = "GBK")
soil_from_liu <- read.csv("data/Ru Wang/HSD_Soil_Nutrition_Liu.csv",header = T,fileEncoding = "GBK")
wang <- read.csv("data/Ru Wang/group.csv",header = T,fileEncoding = "GBK")
towang <- d[,c("qr_AM","qr_EM","qr_BZ","qr_Pn","qr_Pq","qr_fzxb","TagNew",
               "Length(cm)","ProjArea(cm2)","AvgDiam(mm)","Weight.g")]

#首先把她给我的数据的TagNew提出来
library(dplyr)
wang <- left_join(wang,soil_to_liu,by = "SampleID")
library(stringr)
wang$TagNew = str_pad(wang$TagNew,7,side = "left", "0")
#菌根数据
wang <- left_join(wang,towang,by = "TagNew")
#土壤数据
wang <- left_join(wang,soil_from_liu,by = "SampleID")
#搞定！保存一下！耶耶耶！
write.csv(wang,"data/Ru Wang/data_with_roottraits_soilCNP.csv",fileEncoding = "GBK")
