#####
#回归#
#####

#把所有的数据集合起来

######
#加载前步骤
#####
#source("5.neighbors.R")
library(dplyr)

#把所有的数据集合起来

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
  select(-qr_Hbbz, -qr_fzxb, -DATE1, -DATE2,
         -Tag, -Branch, -Species.x, -Opterator,
         -Length, -ProjArea, -SurfArea, -`根系扫描日期`,
         -No_Dat, -`菌根侵染编号`, -`称重日期`, -`错写为`,
         -`Weight.g`, -H)%>%
#重命名一下数据框
  rename(AD = AvgDiam)%>%
#把概率NA的列转为0
  mutate(qr_AM = coalesce(qr_AM, 0),
         GX = as.numeric(as.character(GX)),
         qr_AM = coalesce(qr_AM, 0),
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
#加一列growth_rate
reg <- reg %>%
  mutate(growth_rate = log(DBH2)/log(DBH1),
         growth_rate = ifelse(is.infinite(growth_rate) | is.nan(growth_rate), 0, growth_rate))

###
#把10/20/50分开数据集,把AM=0和EM=0分别数据集，也就是有
reg_AM_10 <- reg %>%
  select(-qr_EM, -pd50_unweigh, -mpd50_unweigh, -mpd50_weigh, -mntd50_unweigh,
         -mntd50_weigh, -pd20_unweigh, -mpd20_unweigh, -mpd20_weigh,
         -mntd20_unweigh, -mntd20_weigh, -shannon_div_20, -invsimpson_div_20,
         -simpson_div_20, -shannon_div_50, -invsimpson_div_50,
         -simpson_div_50, -BD_20, -BD_50, -CBD_20, -CBD_50,-HBD_20, -HBD_50) %>%
  filter(qr_AM > 0)
  
reg_AM_20<- reg %>%
  select(-qr_EM, -pd50_unweigh, -mpd50_unweigh, -mpd50_weigh, -mntd50_unweigh,
         -mntd50_weigh, -pd10_unweigh, -mpd10_unweigh, -mpd10_weigh,
         -mntd10_unweigh, -mntd10_weigh,-shannon_div_10,-invsimpson_div_10,
         -simpson_div_10,-shannon_div_50,-invsimpson_div_50,
         -simpson_div_50, -BD_10, -BD_50, -CBD_10, -CBD_50,-HBD_10, -HBD_50) %>%
  filter(qr_AM > 0)

reg_AM_50 <- reg %>%
  select(-qr_EM, -pd20_unweigh, -mpd20_unweigh, -mpd20_weigh, -mntd20_unweigh,
         -mntd20_weigh, -pd10_unweigh, -mpd10_unweigh, -mpd10_weigh,
         -mntd10_unweigh, -mntd10_weigh,-shannon_div_10,-invsimpson_div_10,
         -simpson_div_10,-shannon_div_20,-invsimpson_div_20,
         -simpson_div_20, -BD_10, -BD_20, -CBD_10, -CBD_20,-HBD_10, -HBD_20) %>%
  filter(qr_AM > 0)

reg_EM_10 <- reg %>%
  select(-qr_AM,-qr_BZ,-qr_Pn, -qr_Pq,
         -pd50_unweigh, -mpd50_unweigh, -mpd50_weigh, -mntd50_unweigh,
         -mntd50_weigh, -pd20_unweigh, -mpd20_unweigh, -mpd20_weigh,
         -mntd20_unweigh, -mntd20_weigh, -shannon_div_20, -invsimpson_div_20,
         -simpson_div_20, -shannon_div_50, -invsimpson_div_50,
         -simpson_div_50, -BD_20, -BD_50, -CBD_20, -CBD_50,-HBD_20, -HBD_50) %>%
  filter(qr_EM > 0)

reg_EM_20 <- reg %>%
  select(-qr_AM,-qr_BZ,-qr_Pn, -qr_Pq,
         -pd50_unweigh, -mpd50_unweigh, -mpd50_weigh, -mntd50_unweigh,
         -mntd50_weigh, -pd10_unweigh, -mpd10_unweigh, -mpd10_weigh,
         -mntd10_unweigh, -mntd10_weigh, -shannon_div_10, -invsimpson_div_10,
         -simpson_div_10, -shannon_div_50, -invsimpson_div_50,
         -simpson_div_50, -BD_10, -BD_50, -CBD_10, -CBD_50,-HBD_10, -HBD_50) %>%
  filter(qr_EM > 0)

reg_EM_50 <- reg %>%
  select(-qr_AM,-qr_BZ,-qr_Pn, -qr_Pq,
         -pd20_unweigh, -mpd20_unweigh, -mpd20_weigh, -mntd20_unweigh,
         -mntd20_weigh, -pd10_unweigh, -mpd10_unweigh, -mpd10_weigh,
         -mntd10_unweigh, -mntd10_weigh, -shannon_div_10, -invsimpson_div_10,
         -simpson_div_10, -shannon_div_20, -invsimpson_div_20,
         -simpson_div_20, -BD_10, -BD_20, -CBD_10, -CBD_20,-HBD_10, -HBD_20) %>%
  filter(qr_EM > 0)



save(reg, reg_AM_10, reg_AM_20, reg_AM_50, 
     reg_EM_10, reg_EM_20, reg_EM_50,
     root_qrl, soil_pred, sp_loc, hsd_neighbor, file = "data/data_for_reg.RData")
#####
#加载之前数据
#####

#load("data/data_for_reg.RData")


#####
#weitao:做回归之前先检查下自相关

#####
#未标准化的数据
#####
library(ggplot2)
library(dplyr)
library(tidyverse)


# 选择数值型变量列
numeric_vars <- reg_AM_20[, c("qr_AM","qr_BZ", "qr_Pn", "qr_Pq",
                              #"qr_EM",
                              "GX", "GY","DBH1", "DBH2",
                              "AD", "SRL", "SRA",
                              "soc", "tn", "tp", "ap", "ph",
                              #"pd50_unweigh", "mpd50_unweigh", "mpd50_weigh", "mntd50_unweigh", "mntd50_weigh",
                              "pd20_unweigh", "mpd20_unweigh", "mpd20_weigh", "mntd20_unweigh", "mntd20_weigh",
                              #"pd10_unweigh", "mpd10_unweigh", "mpd10_weigh", "mntd10_unweigh", "mntd10_weigh",
                              "shannon_div_20", "invsimpson_div_20", "simpson_div_20",
                              #"shannon_div_10", "invsimpson_div_10", "simpson_div_10",
                              #"shannon_div_50", "invsimpson_div_50", "simpson_div_50",
                              "BD_20", "CBD_20", "HBD_20",
                              #"BD_10", "CBD_10", "HBD_10",
                              #"BD_50", "CBD_50", "HBD_50",
                              "growth_rate")]
library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#fc9272", mid = "#f7fcb9", high = "#31a354", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (EM_10)") +
  # 添加相关系数标签
 # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = ifelse(cor_data$value > 0.9 & cor_data$value < 1, paste0(round(cor_data$value, 2), "*"), round(cor_data$value, 2))), 
            vjust = 1, size = 3, color = ifelse(cor_data$value > 0.9 & cor_data$value < 1, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_EM_10.pdf", p, width = 24, height = 12)

library(randomForest)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(qr_AM ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$qr_AM,
  ntreeTry = 1000,
  stepFactor = 1.5,
  improve = 0.01,
  trace = TRUE,
  plot = TRUE
)

best_m <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
print(mtry) # best_m with lowest oob error
print(best_m)

# 使用最佳的mtry值重新训练模型
env_rf_mtry <-
  randomForest(
    qr_AM ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)
var_ip

numeric_vars <- as.data.frame(numeric_vars)

numeric_vars$qr_AM <- ifelse(numeric_vars$qr_AM <= 0, 0.01, ifelse(numeric_vars$qr_AM >= 1, 0.99, numeric_vars$qr_AM))
numeric_vars$qr_EM <- ifelse(numeric_vars$qr_EM >= 1, 0.99, numeric_vars$qr_EM)

library(betareg)
#AM_20
reg$qr_AM <- ifelse(reg$qr_AM <= 0, 0.01, ifelse(reg$qr_AM >= 1, 0.99, reg$qr_AM))
beita_model_forest <- betareg(qr_AM ~
                                 mntd20_unweigh * mntd20_weigh * pd20_unweigh + ap * tp * soc * tn + SRL * SRA * AD + CBD_20, data = reg_AM_20, link = "logit")
summary(beita_model_forest)
library(glmmTMB)
glm_20 <- glmmTMB(qr_AM ~ mntd20_unweigh * mntd20_weigh * pd20_unweigh + ap * tp * soc * tn + SRL * SRA * AD + CBD_20 + (1|Order), reg, family=beta_family)
glm_20 <- glmmTMB(qr_AM ~ mntd20_unweigh + mntd20_weigh + pd20_unweigh + ap + tp + soc + tn + SRL + SRA + AD + CBD_20 + (1|Order), reg, family=beta_family)
#EM_10
beita_model_forest <- betareg(qr_EM ~
                         soc + mpd10_unweigh + tn + pd10_unweigh + DBH1 + shannon_div_10 + invsimpson_div_10 + GX, data = numeric_vars, link = "logit")

summary(beita_model_forest)
vif_values <- vif(beita_model_forest)
print(vif_values)


beita_model_cor <- betareg(qr_AM ~ 
                              AD + SRL + soc+ tp + ap + ph + pd50_unweigh +mpd50_unweigh + mpd50_weigh + mntd50_unweigh + mntd50_weigh + shannon_div_50 + BD_50 + CBD_50 , data = numeric_vars, link = "probit")
beita_model_cor <- betareg(qr_EM ~
                             DBH2 + AD + SRL + soc+ tp + ap + ph + pd10_unweigh +mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 + growth_rate, data = numeric_vars, link = "logit")

summary(beita_model_cor)

vif_values <- vif(beita_model_cor)
print(vif_values)



library(lme4)

reg_AM_10$qr_AM <- ifelse(reg_AM_10$qr_AM <= 0, 0.01, ifelse(reg_AM_10$qr_AM >= 1, 0.99, reg_AM_10$qr_AM))
reg_AM_10$Order <- as.factor(reg_AM_10$Order)
model <- betareg(qr_AM ~ AD + SRL + soc + tp + ap + ph + pd10_unweigh + mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 ,
                 data = reg_AM_10)
glmer(qr_AM ~ AD + SRL + soc+ tp + ap + ph + pd10_unweigh + mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 + Order , data = reg_AM_10, family = betabinomial(link = "logit"))

# 假设 "Order" 是一个因子型变量
reg_AM_10$Order <- factor(reg_AM_10$Order)

# 创建虚拟变量
dummy_vars <- model.matrix(~ Order - 1, data = reg_AM_10)

# 合并虚拟变量到数据框
reg_AM_10 <- cbind(reg_AM_10, dummy_vars)

# 构建Beta回归模型，包括虚拟变量
model <- betareg(qr_AM ~ AD * SRL * soc * tp * ap * ph * pd10_unweigh * mpd10_unweigh * mpd10_weigh * mntd10_unweigh * mntd10_weigh * shannon_div_10 * BD_10 * CBD_10 * Order, 
                 link = "cloglog", data = reg_AM_10)

# 查看模型摘要
summary(model)
