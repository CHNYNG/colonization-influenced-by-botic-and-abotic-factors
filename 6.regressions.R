#####
#回归#
#####

#把所有的数据集合起来


library(dplyr)
reg <- data.frame()
reg <- root_qrl %>%
  left_join(soil_pred[, -which(names(soil_pred) %in% c("GX", "GY"))],
            by = "TagNew", suffix = c("", "_soil"), relationship = "many-to-many") %>%
  left_join(sp_loc[, -which(names(sp_loc) %in% c("GX", "GY"))],
            by = "TagNew", suffix = c("", "_sp_loc"), relationship = "many-to-many") %>%
  left_join(hsd_neighbor, 
            by = "TagNew", suffix = c("", "_hsd_neighbor"), relationship = "many-to-many")
#整理一下reg

reg <- reg %>%
  #筛一下reg，留下Branch==0的数值
  filter(Branch == 0) %>%
  #删除重复项
  distinct() %>%
  #删除用不到的列
  select(-Numbers,-qr_Hbbz, -qr_fzxb, -DATE1, -DATE2,
         -Tag, -Branch, -Species.x, -Opterator,
         -Length, -ProjArea, -SurfArea, -`根系扫描日期`,
         -No_Dat, -`菌根侵染编号`, -`称重日期`, -`错写为`,
         -`Weight.g`, -H)%>%
#重命名一下数据框
  rename(AD = AvgDiam)%>%
#把概率NA的列转为0
  mutate(qr_AM = coalesce(qr_AM, 0),
         qr_EM = coalesce(qr_EM, 0),
         qr_BZ = coalesce(qr_BZ, 0),
         qr_Pn = coalesce(qr_Pn, 0),
         qr_Pq = coalesce(qr_Pq, 0),
         SRL = coalesce(SRL, 0),
         AD = coalesce(AD, 0),
         SRA = coalesce(SRA, 0),
         DBH1 =coalesce(DBH1, 0))#%>%
  #删除H和DBH1为NA的行
  #filter(!is.na(H))



#####
#weitao:做回归之前先检查下自相关
#检查一下可能相关的变量的自相关关系
library(GGally)
reg_myco <- reg[,c("qr_AM", "qr_EM", "qr_BZ", "qr_Pn","qr_Pq" )]
p_reg_myco <- ggpairs(reg_myco, title="correlogram of mycorrhiza") 
#qr_AM与qr_BZ、qr_Pn不显著
#qr_EM与qr_BZ、qr_Pn、qr_Pq均不显著
#qr_Pn与qr_Pq不显著
ggsave("pic/p_reg_myco.pdf", plot = p_reg_myco, width = 11.69, height = 8.27, units = "in", device = "pdf")

#关于生长的变量需要再算一个growth_rate
reg <- reg %>%
  mutate(growth_rate= log(DBH2)/log(DBH1),
  growth_rate = ifelse(is.infinite(growth_rate), 0,
                       growth_rate))

reg_gr <- reg[,c("DBH1", "DBH2", "growth_rate")]
p_reg_gr<- ggpairs(reg_gr, title="correlogram of growth")
#三者相关性显著
ggsave("pic/p_reg_gr.pdf", plot = p_reg_gr, width = 11.69, height = 8.27, units = "in", device = "pdf")


#关于根系性状的变量的相关
reg_root <- reg[,c("AD", "SRL", "SRA")]
p_reg_root <- ggpairs(reg_root, title="correlogram of root")
#三者相关性显著
ggsave("pic/p_reg_root.pdf", plot = p_reg_root, width = 11.69, height = 8.27, units = "in", device = "pdf")


#关于土壤元素
reg_soil <- reg[,c("soc", "tn", "tp", "ap", "ph")]
p_reg_soil <- ggpairs(reg_soil, title="correlogram of soil")
#其中ap和soc的相关性不显著
ggsave("pic/p_reg_soil.pdf", plot = p_reg_soil, width = 11.69, height = 8.27, units = "in", device = "pdf")

#关于谱系（我好懒，想把它们搞一起
reg_phylo <- reg[,c("pd20_unweigh", "mpd20_unweigh", "mpd20_weigh", "mntd20_unweigh", "mntd20_weigh",
                    "pd10_unweigh", "mpd10_unweigh", "mpd10_weigh", "mntd10_unweigh", "mntd10_weigh",
                    "pd50_unweigh", "mpd50_unweigh", "mpd50_weigh", "mntd50_unweigh", "mntd50_weigh")]

p_reg_phylo <- ggpairs(reg_phylo, title="correlogram of phylo")
ggsave("pic/p_reg_phylo.pdf", plot = p_reg_phylo, width = 11.69, height = 8.27, units = "in", device = "pdf")

#关于邻体多样性
reg_div <- reg[,c("shannon_div_20", "invsimpson_div_20", "simpson_div_20",
                  "shannon_div_10", "invsimpson_div_10", "simpson_div_10",
                  "shannon_div_50", "invsimpson_div_50", "simpson_div_50",)]

p_reg_div <- ggpairs(reg_div, title="correlogram of divsity")
ggsave("pic/p_reg_div.pdf", plot = p_reg_div, width = 11.69, height = 8.27, units = "in", device = "pdf")

#关于邻体影响力
reg_nei <- reg[,c("BD_20", "CBD_20", "HBD_20",
                  "BD_10", "CBD_10", "HBD_10",
                  "BD_50", "CBD_50", "HBD_50")]

p_reg_nei <- ggpairs(reg_nei, title="correlogram of neighbors")
ggsave("pic/p_reg_nei.pdf", plot = p_reg_nei, width = 11.69, height = 8.27, units = "in", device = "pdf")

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
