#####
#回归#
#####

#把所有的数据集合起来


library(dplyr)
reg <- left_join(root_qrl, soil_pred[,-which(names(soil_pred) %in% c("GX","GY"))], by = "TagNew",relationship = "many-to-many")
reg <- left_join(reg, sp_loc[,-which(names(sp_loc) %in% c("GX","GY"))], by = "TagNew", relationship = "many-to-many")
#整理一下reg

reg <- reg %>%
  #筛一下reg，留下Branch==0的数值
  filter(Branch == 0) %>%
  #删除重复项
  distinct() %>%
  #删除用不到的列
  select(-Numbers, -gcsys, -EMsys, -jssys, -BZ, -Pn, -Pq,
         -Hbbz, -Fzxb, -qr_Hbbz, -qr_fzxb, -DATE1,
         -DATE2, -Tag, -Branch, -Species.x, -Opterator.x,
         -Opterator.y, -`Length(cm)`, -`根系扫描日期`,
         -No_Dat, -`菌根侵染编号`, -`称重日期`, -`错写为`,
         `ProjArea(cm2)`, `SurfArea(cm2)`,-`Weight.g`)%>%
#重命名一下数据框
  rename(SRL = "SRL(gcm)",
         AD = "AvgDiam(mm)")%>%
#把概率NA的列转为0
  mutate(qr_AM = coalesce(qr_AM, 0),
         qr_EM = coalesce(qr_EM, 0),
         qr_BZ = coalesce(qr_BZ, 0),
         qr_Pn = coalesce(qr_Pn, 0),
         qr_Pq = coalesce(qr_Pq, 0),
         SRL = coalesce(SRL, 0),
         AD = coalesce(AD, 0),
         SRA = coalesce(SRA, 0))%>%
  #删除H和DBH1为NA的行
  filter(!is.na(H))
###整理一下，做个重新的排序
new_order <- c("TagNew", "Qudrat", "Latin", "Genus", "Family", "Order", "M_Type", "Type", "GF", "GX", "GY",
               "Status1" , "Status2", "DBH1", "DBH2", "H",
               "qr_AM", "qr_EM", "qr_BZ", "qr_Pn", "qr_Pq" ,  "AD","SRL", "SRA",
               "soc", "tn", "tp", "ap", "ph",
               "pd20_unweigh", "mpd20_unweigh", "mpd20_weigh", "mntd20_unweigh", "mntd20_weigh",
               "shannon_div_20", "invsimpson_div_20", "simpson_div_20" )

reg <- reg %>%
  select(all_of(new_order))

#做回归之前先检查下自相关
correlation_matrix <- cor(reg[,14:37])
# 查看相关性矩阵
print(correlation_matrix)

#按照侵染率、根系、土壤、植物、谱系参数、多样性指数
#分别计算变量间的相关性
Scatterplot <- ggpairs(reg[,14:37], title="correlogram with ggpairs()") 

#β回归
library(tidyverse)
beita_reg <- lm(qr_AM ~ shannon_div_20 + SRL + AD + mpd20_weigh + soc + tn + tp + ap + ph+ DBH2, data = reg)
step_beita_reg <- step(beita_reg)
summary(beita_reg)
summary(step_beita_reg)
