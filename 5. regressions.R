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
  select(-gcsys, -EMsys, jssys, -BZ, -Pn, -Pq,
         -Hbbz, -Fzxb, -qr_Hbbz, -qr_fzxb,
         )
#重命名一下数据框
  rename(SRL = "SRL(gcm)",
         AD = "AvgDiam(mm)")
#做回归之前先检查下自相关

#β回归
library(tidyverse)
beita_reg <- lm(qr_AM ~ shannon_div_20 + SRL + AD + mpd20_weigh + soc + tn + tp + ap + ph+ DBH2, data = reg)
step_beita_reg <- step(beita_reg)
summary(beita_reg)
summary(step_beita_reg)
