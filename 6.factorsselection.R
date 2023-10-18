#把所有的数据集合起来

######
#加载前步骤
#####
#source("5.neighbors.R")
library(dplyr)

#把所有的数据集合起来
load("data/env_weitao.RData")

reg <- data.frame()
reg <- root_qrl %>%
  left_join(soil_pred[, -which(names(soil_pred) %in% c("GX", "GY"))],
            by = "TagNew", suffix = c("", "_soil"), relationship = "many-to-many") %>%
  left_join(sp_loc[, -which(names(sp_loc) %in% c("GX", "GY"))],
            by = "TagNew", suffix = c("", "_sp_loc"), relationship = "many-to-many") %>%
  left_join(hsd_neighbor, 
            by = "TagNew", suffix = c("", "_hsd_neighbor"), relationship = "many-to-many")
reg <- cbind(reg,pd_ind_all,env)
#整理一下reg

reg <- reg %>%
  #筛一下reg，留下Branch==0的数值
  # filter(Branch == 0) %>%
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
  mutate(#qr_AM = coalesce(qr_AM, 0),
    GX = as.numeric(as.character(GX)),
    #qr_AM = coalesce(qr_AM, 0),
    #qr_EM = coalesce(qr_EM, 0),
    #qr_BZ = coalesce(qr_BZ, 0),
    #qr_Pn = coalesce(qr_Pn, 0),
    #qr_Pq = coalesce(qr_Pq, 0),
    #SRL = coalesce(SRL, 0),
    #AD = coalesce(AD, 0),
    #SRA = coalesce(SRA, 0),
    #DBH1 =coalesce(DBH1, 0)
  )#%>%
#删除H和DBH1为NA的行
#filter(!is.na(H))
#加一列growth_rate
#删除可能有tag错误的行（这样的行qr_AM、EM和AD、SRL、SRA都为0，应该有10行左右
#reg <- reg[!(reg$qr_AM == 0 & reg$qr_EM == 0 & reg$AD == 0 & reg$SRL == 0 & reg$SRA == 0), ]

#reg <- reg %>%
#  mutate(growth_rate = log(DBH2)/log(DBH1),
#         growth_rate = ifelse(is.infinite(growth_rate) | is.nan(growth_rate), 0, growth_rate))

pca_data <- reg[, c("soc","tn","tp","ap","ph", "moisture")]
pca_result <- prcomp(pca_data, scale = TRUE)
pca_result <- pca_result$x
pca_results <- data.frame(PC1 = pca_result[, 1], PC2 = pca_result[, 2], reg[,c("qr_AM","M_Type")])
summary(pca_result)
reg <- reg %>%
  mutate(soil_pc1 = pca_results[,"PC1"],
         soil_pc2 = pca_results[,"PC2"])

#####
#数据分类
#####

#####
#给多度分类
#####


hsd_ba <- hsd_alive_singlebr %>%
  group_by(scientific.name) %>%
  summarise(abundance = n(),
            dbh_multi = sum(dbh_multi)) %>%
  rename(Latin = scientific.name) %>%
  mutate(rel_abun = abundance / sum(abundance),
         rel_dbh_multi = dbh_multi / sum(dbh_multi)) %>%
  arrange(desc(abundance)) %>%
  mutate(cumulative_abun = cumsum(rel_abun),
         cumulative_dbh_multi = cumsum(rel_dbh_multi)) %>%
  mutate(
    sptype1 = ifelse(
      cumulative_abun <= 0.2,
      "dominant",
      ifelse(cumulative_abun > 0.2 &
               cumulative_abun <= 0.9,
             "common", "rare")
    ),
    sptype2 = ifelse(
      cumulative_dbh_multi <= 0.1,
      "dominant",
      ifelse(cumulative_dbh_multi > 0.1 &
               cumulative_dbh_multi <= 0.9,
             "common", "rare")
    )
  )
###这里只有显脉新木姜子和黄果厚壳桂两个优势种，检查完数据也就这两个！！！！师兄的数据有问题！！！！
####(师兄的数据做出来也只有3个，嘤嘤嘤，不是我的错)

reg <- reg %>%
  left_join(hsd_ba, by = "Latin") %>%
  select(-abundance, -dbh_multi, -rel_abun, rel_dbh_multi, -cumulative_abun, 
         -cumulative_dbh_multi)


#####
#给dbh分类
#####
#reg <- reg[order(reg$Latin,reg$DBH2),]

#reg <- reg %>%
#  group_by(Latin) %>%
#  mutate(DBH_order = row_number(),
#         DBH_order = as.character(DBH_order))
#reg <- reg %>%
#  ungroup()

hsd_alive_singlebr_dbh <- hsd_alive_singlebr %>%
  rename(Latin=scientific.name)%>%
  group_by(Latin) %>%
  mutate(dbh_multi_rank = rank(dbh_multi)) %>%
  ungroup()

hsd_alive_singlebr_dbh <- hsd_alive_singlebr_dbh %>%
  group_by(Latin) %>%
  mutate(total_count = n()) %>%
  ungroup() %>%
  mutate(seedling_cutoff = total_count * 0.25) %>%
  group_by(Latin) %>%
  mutate(history = ifelse(dbh_multi_rank <= seedling_cutoff, "seedling", "adult")) %>%
  ungroup() %>%
  select(-total_count, -seedling_cutoff, -dbh_multi_rank, -Latin)

reg <- reg %>%
  left_join(hsd_alive_singlebr_dbh %>% select(TagNew, history), by = "TagNew")%>%
  distinct()

#####
#给土壤聚类
#####

#####
#第一种聚类
#####

# 将数据分为k个簇（k是您选择的簇数）
library(cluster)
# 选择要聚类的变量，可以根据需要进行选择
soil_cluster <- soil_pred[, c("soc", "tn", "tp", "ap", "ph", "moisture")]

# 执行K均值聚类
k <- 3  # 设置簇的数量
kmeans_result <- kmeans(soil_cluster, centers = k)

# 查看每个样本所属的簇
#cluster_assignments <- kmeans_result$cluster

reg <- reg %>%
  mutate(soil_cluster = kmeans_result$cluster)


#####
#将数据进行scale
#####

scale_to_01 <- function(x) {
  # 使用na.rm参数忽略NA值
(x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
scale_to_01(c(1,1,1,1,2))
scale_data <- reg %>%
  select(qr_AM, qr_BZ, qr_Pn, qr_Pq,qr_EM,
         GX, GY,DBH1, DBH2, rel_dbh_multi, AD, SRL, SRA,
         soc, tn, tp, ap, ph, moisture,
         elevation, aspect, slope, convexity,
         pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
         pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
         pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
         totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
         totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
         totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
         shannon_div_20, invsimpson_div_20, simpson_div_20,
         shannon_div_10, invsimpson_div_10, simpson_div_10,
         shannon_div_50, invsimpson_div_50, simpson_div_50,
         BD_20, CBD_20, HBD_20,
         BD_10, CBD_10, HBD_10,
         BD_50, CBD_50, HBD_50)
scale_data <- as.matrix(scale_data)
scale_data <- apply(scale_data,2,scale_to_01)
scale_data <- as.data.frame(scale_data)
reg_sc <- reg %>%
  select( -qr_AM, -qr_BZ, -qr_Pn, -qr_Pq, -qr_EM, 
          -GX, -GY, -DBH1, -DBH2, -rel_dbh_multi, -AD, -SRL, -SRA,
          -soc, -tn, -tp, -ap, -ph, -moisture, 
          -elevation, -aspect, -slope, -convexity,
          -pd50_unweigh, -mpd50_unweigh, -mpd50_weigh, -mntd50_unweigh, -mntd50_weigh,
          -pd20_unweigh, -mpd20_unweigh, -mpd20_weigh, -mntd20_unweigh, -mntd20_weigh,
          -pd10_unweigh, -mpd10_unweigh, -mpd10_weigh, -mntd10_unweigh, -mntd10_weigh,
          -totpd_20, -avepd_20, -minpd_20, -apd_20, -ntpd_20,
          -totpd_10, -avepd_10, -minpd_10, -apd_10, -ntpd_10,
          -totpd_50, -avepd_50, -minpd_50, -apd_50, -ntpd_50,
          -shannon_div_20, -invsimpson_div_20, -simpson_div_20,
          -shannon_div_10, -invsimpson_div_10, -simpson_div_10,
          -shannon_div_50, -invsimpson_div_50, -simpson_div_50,
          -BD_20, -CBD_20, -HBD_20,
          -BD_10, -CBD_10, -HBD_10,
          -BD_50, -CBD_50, -HBD_50)%>%
  bind_cols(scale_data)

rm(scale_data)

save(reg, file = "data/data_for_reg.RData")

#####
#选择回归用的变量
#####

#20
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,
  am,
  #qr_BZ, qr_Pn, qr_Pq,
  #qr_EM,
  #GX, GY,DBH1,
  DBH2, rel_dbh_multi, RDi,
  AD, SRL, SRA,
  soc, tn, tp, ap, ph, moisture, 
  elevation, aspect, slope, convexity,
  soil_pc1, soil_pc2,
  #pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
  pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
  #pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
  totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
  #totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
  #totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
  shannon_div_20, invsimpson_div_20, simpson_div_20,
  #shannon_div_10, invsimpson_div_10, simpson_div_10,
  #shannon_div_50, invsimpson_div_50, simpson_div_50,
  BD_20, CBD_20, HBD_20
  #BD_10, CBD_10, HBD_10,
  #BD_50, CBD_50, HBD_50,
)

numeric_vars <- na.omit(numeric_vars)


library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#f7fcb9", mid = "#cbc9e2", high = "#f7fcb9", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (am_20)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_AM_20.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(am ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$am,
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
    am ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)

print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)
#apd_20 + avepd_20 + minpd_20 + soc + pd20_unweigh + tp + RDi + tn + moisture + ap

glm_20_M_1 <- glmmTMB(am ~  apd_20 + minpd_20 + soc + pd20_unweigh + tp + tn + moisture + ap + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20_M_1b <- update(glm_20_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_1b)
glm_20_M_2 <- glmmTMB(am ~  apd_20 + minpd_20 + soc + pd20_unweigh + tp + tn + moisture + ap + (1|resource), reg_sc, family=beta_family)
glm_20_M_2b <- update(glm_20_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_2b)
glm_20_M_3 <- glmmTMB(am ~  apd_20 + minpd_20 + soc + pd20_unweigh + tp + tn + moisture + ap + (1|sptype2), reg_sc, family=beta_family)
glm_20_M_3b <- update(glm_20_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_3b)
#mntd20_unweigh + aspect + AD + avepd_20 + DBH2 + pd20_unweigh + apd_20 + minpd_20 + convexity
glm_20_P_1 <- glmmTMB(am ~  mntd20_unweigh + aspect + AD + avepd_20 + DBH2 + pd20_unweigh + apd_20 + minpd_20 + convexity + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20_P_1b <- update(glm_20_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_1b)
glm_20_P_2 <- glmmTMB(am ~  mntd20_unweigh + aspect + AD + avepd_20 + DBH2 + pd20_unweigh + apd_20 + minpd_20 + convexity + (1|resource), reg_sc, family=beta_family)
glm_20_P_2b <- update(glm_20_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_2b)
glm_20_P_3 <- glmmTMB(am ~  mntd20_unweigh + aspect + AD + avepd_20 + DBH2 + pd20_unweigh + apd_20 + minpd_20 + convexity + (1|sptype2), reg_sc, family=beta_family)
glm_20_P_3b <- update(glm_20_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_3b)

#10
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,
    am,
    #qr_BZ, qr_Pn, qr_Pq,
    #qr_EM,
    #GX, GY,DBH1,
    DBH2, rel_dbh_multi, RDi,
    AD, SRL, SRA,
    soc, tn, tp, ap, ph, moisture, 
    elevation, aspect, slope, convexity,
    soil_pc1, soil_pc2,
    #pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
    #pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
    pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
    #totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
    totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
    #totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
    #shannon_div_20, invsimpson_div_20, simpson_div_20,
    shannon_div_10, invsimpson_div_10, simpson_div_10,
    #shannon_div_50, invsimpson_div_50, simpson_div_50,
    #BD_20, CBD_20, HBD_20
    BD_10, CBD_10, HBD_10,
    #BD_50, CBD_50, HBD_50,
  )

numeric_vars <- na.omit(numeric_vars)


library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#cbc9e2", mid = "#f7fcb9", high = "#cbc9e2", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (am_10)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_AM_10.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(am ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$am,
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
    am ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)

print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
  head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)

#apd_10 + RDi + totpd_10 + moisture + minpd_10 + shannon_div_10 + soc
glm_10_M_1 <- glmmTMB(am ~  apd_10 + RDi + totpd_10 + moisture + minpd_10 + shannon_div_10 + soc   + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10_M_1b <- update(glm_10_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_1b)
glm_10_M_2 <- glmmTMB(am ~  apd_10 + RDi + totpd_10 + moisture + minpd_10 + shannon_div_10 + soc + (1|resource), reg_sc, family=beta_family)
glm_10_M_2b <- update(glm_10_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_2b)
glm_10_M_3 <- glmmTMB(am ~  apd_10 + RDi + totpd_10 + moisture + minpd_10 + shannon_div_10 + soc + (1|sptype2), reg_sc, family=beta_family)
glm_10_M_3b <- update(glm_10_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_3b)

#mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10
glm_10_P_1 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10_P_1b <- update(glm_10_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_1b)
glm_10_P_2 <- glmmTMB(am ~ mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|resource), reg_sc, family=beta_family)
glm_10_P_2b <- update(glm_10_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_2b)
glm_10_P_3 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|sptype2), reg_sc, family=beta_family)
glm_10_P_3b <- update(glm_10_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_3b)


#50
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,
    am,
    #qr_BZ, qr_Pn, qr_Pq,
    #qr_EM,
    #GX, GY,DBH1,
    DBH2, rel_dbh_multi, RDi,
    AD, SRL, SRA,
    soc, tn, tp, ap, ph, moisture, 
    elevation, aspect, slope, convexity,
    soil_pc1, soil_pc2,
    pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
    #pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
    #pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
    #totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
    #totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
    totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
    #shannon_div_20, invsimpson_div_20, simpson_div_20,
    #shannon_div_10, invsimpson_div_10, simpson_div_10,
    shannon_div_50, invsimpson_div_50, simpson_div_50,
    #BD_20, CBD_20, HBD_20
    #BD_10, CBD_10, HBD_10,
    BD_50, CBD_50, HBD_50,
  )

numeric_vars <- na.omit(numeric_vars)


library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#cbc9e2", mid = "#f7fcb9", high = "#cbc9e2", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (am_50)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_AM_50.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(am ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$am,
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
    am ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)
print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
  head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)
#apd_50 + tp + minpd_50 + totpd_50 + AD + ap

glm_50_M_1 <- glmmTMB(am ~  apd_50 + tp + minpd_50 + totpd_50 + AD + ap  + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50_M_1b <- update(glm_50_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_1b)
glm_50_M_2 <- glmmTMB(am ~  apd_50 + tp + minpd_50 + totpd_50 + AD + ap + (1|resource), reg_sc, family=beta_family)
glm_50_M_2b <- update(glm_50_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_2b)
glm_50_M_3 <- glmmTMB(am ~  apd_50 + tp + minpd_50 + totpd_50 + AD + ap + (1|sptype2), reg_sc, family=beta_family)
glm_50_M_3b <- update(glm_50_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_3b)
#apd_50 + AD + aspect + minpd_50 + pd50_unweigh + DBH2 + mntd50_unweigh + tp  
glm_50_P_1 <- glmmTMB(am ~  apd_50 + AD + aspect + minpd_50 + pd50_unweigh + DBH2 + mntd50_unweigh + tp   + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50_P_1b <- update(glm_50_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_1b)
glm_50_P_2 <- glmmTMB(am ~  apd_50 + AD + aspect + minpd_50 + pd50_unweigh + DBH2 + mntd50_unweigh + tp   + (1|resource), reg_sc, family=beta_family)
glm_50_P_2b <- update(glm_50_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_2b)
glm_50_P_3 <- glmmTMB(am ~  apd_50 + AD + aspect + minpd_50 + pd50_unweigh + DBH2 + mntd50_unweigh + tp   + (1|sptype2), reg_sc, family=beta_family)
glm_50_P_3b <- update(glm_50_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_3b)

#####
#关于em
#####
#20
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,
    #am,
    em,
    #qr_BZ, qr_Pn, qr_Pq,
    #qr_EM,
    #GX, GY,DBH1,
    DBH2, rel_dbh_multi, RDi,
    AD, SRL, SRA,
    soc, tn, tp, ap, ph, moisture, 
    elevation, aspect, slope, convexity,
    soil_pc1, soil_pc2,
    #pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
    pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
    #pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
    totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
    #totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
    #totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
    shannon_div_20, invsimpson_div_20, simpson_div_20,
    #shannon_div_10, invsimpson_div_10, simpson_div_10,
    #shannon_div_50, invsimpson_div_50, simpson_div_50,
    BD_20, CBD_20, HBD_20
    #BD_10, CBD_10, HBD_10,
    #BD_50, CBD_50, HBD_50,
  )

numeric_vars <- na.omit(numeric_vars)

library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#cbc9e2", mid = "#f7fcb9", high = "#cbc9e2", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (em_20)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_EM_10.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(em ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$em,
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
    em ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)
print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
  head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)

glm_20em_M_1 <- glmmTMB(em ~  tn + RDi + apd_20 + soil_pc2 + elevation + slope + soc + pd20_unweigh + SRL + CBD_20 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20em_M_1b <- update(glm_20em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_20em_M_1b)
glm_20em_M_2 <- glmmTMB(em ~   tn + RDi + apd_20 + soil_pc2 + elevation + slope + soc + pd20_unweigh + SRL + CBD_20
 + (1|resource), reg_sc, family=beta_family)
glm_20em_M_2b <- update(glm_20em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_2b)
glm_20em_M_3 <- glmmTMB(em ~   tn + RDi + apd_20 + soil_pc2 + elevation + slope + soc + pd20_unweigh + SRL + CBD_20 + (1|sptype2), reg_sc, family=beta_family)
glm_20em_M_3b <- update(glm_20em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_3b)
glm_20em_M_4 <- glmmTMB(em ~ tn + RDi + apd_20 + soil_pc2 + elevation + slope + soc + pd20_unweigh + SRL + CBD_20 , reg_sc, family=beta_family)
glm_20em_M_4b <- update(glm_20em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_4b)
glm_20em_P_1 <- glmmTMB(em ~  mntd20_weigh + tn + soc + slope + mpd20_unweigh + aspect + ap + apd_20 + convexity + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20em_P_1b <- update(glm_20em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_1b)
glm_20em_P_2 <- glmmTMB(em ~  mntd20_weigh + tn + soc + slope + mpd20_unweigh + aspect + ap + apd_20 + convexity + (1|resource), reg_sc, family=beta_family)
glm_20em_P_2b <- update(glm_20em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_2b)
glm_20em_P_3 <- glmmTMB(em ~  mntd20_weigh + tn + soc + slope + mpd20_unweigh + aspect + ap + apd_20 + convexity + (1|sptype2), reg_sc, family=beta_family)
glm_20em_P_3b <- update(glm_20em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_3b)
glm_20em_P_4 <- glmmTMB(em ~  mntd20_weigh + tn + soc + slope + mpd20_unweigh + aspect + ap + apd_20 + convexity , reg_sc, family=beta_family)
glm_20em_P_4b <- update(glm_20em_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_4)

#10
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,am,
    em,
    #qr_BZ, qr_Pn, qr_Pq,
    #qr_EM,
    #GX, GY,DBH1,
    DBH2, rel_dbh_multi, RDi,
    AD, SRL, SRA,
    soc, tn, tp, ap, ph, moisture, 
    elevation, aspect, slope, convexity,
    soil_pc1, soil_pc2,
    #pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
    #pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
    pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
    #totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
    totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
    #totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
    #shannon_div_20, invsimpson_div_20, simpson_div_20,
    shannon_div_10, invsimpson_div_10, simpson_div_10,
    #shannon_div_50, invsimpson_div_50, simpson_div_50,
    #BD_20, CBD_20, HBD_20
    BD_10, CBD_10, HBD_10,
    #BD_50, CBD_50, HBD_50,
  )

numeric_vars <- na.omit(numeric_vars)


library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#cbc9e2", mid = "#f7fcb9", high = "#cbc9e2", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (em_10)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_EM_10.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(em ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$em,
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
    em ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)
print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
  head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)

glm_10em_M_1 <- glmmTMB(em ~  tn + mpd10_weigh + RDi + soil_pc2 + moisture + avepd_10 + SRL + pd10_unweigh  + rel_dbh_multi + totpd_10 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10em_M_1b <- update(glm_10em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_10em_M_1b)
glm_10em_M_2 <- glmmTMB(em ~   tn + mpd10_weigh + RDi + soil_pc2 + moisture + avepd_10 + SRL + pd10_unweigh  + rel_dbh_multi + totpd_10 + (1|resource), reg_sc, family=beta_family)
glm_10em_M_2b <- update(glm_10em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_2b)
glm_10em_M_3 <- glmmTMB(em ~   tn + mpd10_weigh + RDi + soil_pc2 + moisture + avepd_10 + SRL + pd10_unweigh  + rel_dbh_multi + totpd_10 + (1|sptype2), reg_sc, family=beta_family)
glm_10em_M_3b <- update(glm_10em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_3b)
glm_10em_M_4 <- glmmTMB(em ~   tn + mpd10_weigh + RDi + soil_pc2 + moisture + avepd_10 + SRL + pd10_unweigh  + rel_dbh_multi + totpd_10 , reg_sc, family=beta_family)
glm_10em_M_4b <- update(glm_10em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_4b)
glm_10em_P_1 <- glmmTMB(em ~  mntd10_weigh + tn + slope + soc + mpd10_unweigh + ap + shannon_div_10 + soil_pc1 + elevation + AD + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10em_P_1b <- update(glm_10em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_1b)
glm_10em_P_2 <- glmmTMB(em ~  mntd10_weigh + tn + slope + soc + mpd10_unweigh + ap + shannon_div_10 + soil_pc1 + elevation + AD + (1|resource), reg_sc, family=beta_family)
glm_10em_P_2b <- update(glm_10em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_2b)
glm_10em_P_3 <- glmmTMB(em ~  mntd10_weigh + tn + slope + soc + mpd10_unweigh + ap + shannon_div_10 + soil_pc1 + elevation + AD + (1|sptype2), reg_sc, family=beta_family)
glm_10em_P_3b <- update(glm_10em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_3b)

#50
# 选择数值型变量列
numeric_vars <- reg_sc %>%
  select(#qr_AM,am,
    em,
    #qr_BZ, qr_Pn, qr_Pq,
    #qr_EM,
    #GX, GY,DBH1,
    DBH2, rel_dbh_multi, RDi,
    AD, SRL, SRA,
    soc, tn, tp, ap, ph, moisture, 
    elevation, aspect, slope, convexity,
    soil_pc1, soil_pc2,
    pd50_unweigh, mpd50_unweigh, mpd50_weigh, mntd50_unweigh, mntd50_weigh,
    #pd20_unweigh, mpd20_unweigh, mpd20_weigh, mntd20_unweigh, mntd20_weigh,
    #pd10_unweigh, mpd10_unweigh, mpd10_weigh, mntd10_unweigh, mntd10_weigh,
    #totpd_20, avepd_20, minpd_20, apd_20, ntpd_20, 
    #totpd_10, avepd_10, minpd_10, apd_10, ntpd_10, 
    totpd_50, avepd_50, minpd_50, apd_50, ntpd_50, 
    #shannon_div_20, invsimpson_div_20, simpson_div_20,
    #shannon_div_10, invsimpson_div_10, simpson_div_10,
    shannon_div_50, invsimpson_div_50, simpson_div_50,
    #BD_20, CBD_20, HBD_20
    #BD_10, CBD_10, HBD_10,
    BD_50, CBD_50, HBD_50,
  )

numeric_vars <- na.omit(numeric_vars)


library(ggplot2)
library(reshape2)

# 计算相关性矩阵

cor_matrix <- cor(numeric_vars, use = "complete.obs", method = "spearman")

# 将相关性矩阵转换为长格式
cor_data <- melt(cor_matrix)

# 创建相关性图
p <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#cbc9e2", mid = "#f7fcb9", high = "#cbc9e2", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot (em_50)") +
  # 添加相关系数标签
  # geom_text(aes(label = round(value, 2)), vjust = 1, size = 3, color = "black")+
  
  geom_text(aes(label = round(cor_data$value, 2)), 
            vjust = 1, size = 3, color = ifelse(abs(cor_data$value) > 0.75, "red", "black"), show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # 调整横坐标字体大小
        axis.text.y = element_text(size = 12),  # 调整纵坐标字体大小
        axis.title.x = element_blank(),  # 去掉横坐标标题
        axis.title.y = element_blank())  # 去掉纵坐标标题
print(p)
ggsave("pic/correlation_plot_EM_50.pdf", p, width = 24, height = 12)

library(randomForest)
library(tibble)
####使用随机森林选择变量
# 使用你的数据集和目标变量
set.seed(0725)
env_rf <- randomForest(em ~ ., 
                       data = numeric_vars, ntree = 1000)
print(env_rf)

# 进行参数调优
mtry <- tuneRF(
  numeric_vars,
  numeric_vars$em,
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
    em ~ .,
    data = numeric_vars,
    mtry = best_m,
    importance = TRUE,
    ntree = 1000
  )
print(env_rf_mtry)

# 可以继续使用你的代码进行特征选择
varImpPlot(env_rf_mtry)
var_ip_MSE <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(percentage_increase_in_MSE)) %>% 
  head(10)
print(var_ip_MSE$variables)

var_ip_purity <- importance(env_rf_mtry) %>% 
  as.data.frame() %>% 
  rownames_to_column("variables") %>% 
  rename(
    'percentage_increase_in_MSE' = `%IncMSE`
  ) %>% 
  arrange(desc(IncNodePurity)) %>% 
  head(10)
print(var_ip_purity$variables)

numeric_vars <- as.data.frame(numeric_vars)

glm_50em_M_1 <- glmmTMB(em ~  mntd50_unweigh + RDi + SRL + mntd50_weigh + DBH2 + invsimpson_div_50 + BD_50 + totpd_50 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50em_M_1b <- update(glm_50em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_50em_M_1b)
glm_50em_M_2 <- glmmTMB(em ~  mntd50_unweigh + RDi + SRL + mntd50_weigh + DBH2 + invsimpson_div_50 + BD_50 + totpd_50 , reg_sc, family=beta_family)
glm_50em_M_2b <- update(glm_50em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_2b)
glm_50em_M_3 <- glmmTMB(em ~  mntd50_unweigh + RDi + SRL + mntd50_weigh + DBH2 + invsimpson_div_50 + BD_50 + totpd_50 + (1|sptype2), reg_sc, family=beta_family)
glm_50em_M_3b <- update(glm_50em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_3b)
glm_50em_P_1 <- glmmTMB(em ~  mntd50_weigh + slope + mntd50_unweigh + tn + ap + soil_pc1 +  convexity +  aspect +  SRL + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50em_P_1b <- update(glm_50em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_1b)
glm_50em_P_2 <- glmmTMB(em ~  mntd50_weigh + slope + mntd50_unweigh + tn + ap + soil_pc1 +  convexity +  aspect +  SRL + (1|resource), reg_sc, family=beta_family)
glm_50em_P_2b <- update(glm_50em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_2b)
glm_50em_P_3 <- glmmTMB(em ~  mntd50_weigh + slope + mntd50_unweigh + tn + ap + soil_pc1 +  convexity +  aspect +  SRL + (1|sptype2), reg_sc, family=beta_family)
glm_50em_P_3b <- update(glm_50em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_3b)


