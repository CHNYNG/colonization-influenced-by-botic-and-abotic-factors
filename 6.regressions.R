#####
#回归#
#####

#把所有的数据集合起来

library(dplyr)
=======
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
  select(-Numbers,-qr_Hbbz, -qr_fzxb, -DATE1, -DATE2,
         -Tag, -Branch, -Species.x, -Opterator,
         -Length, -ProjArea, -SurfArea, -`根系扫描日期`,
         -No_Dat, -`菌根侵染编号`, -`称重日期`, -`错写为`,
         -`Weight.g`, -H)%>%
#重命名一下数据框
  rename(AD = AvgDiam)%>%
#把概率NA的列转为0
  mutate(qr_AM = coalesce(qr_AM, 0),
  mutate(GX = as.numeric(as.character(GX)),
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
=======
reg <- reg %>%
  mutate(growth_rate = log(DBH2)/log(DBH1),
         growth_rate = ifelse(is.infinite(growth_rate) | is.nan(growth_rate), 0, growth_rate))


save(reg, root_qrl, soil_pred, sp_loc, hsd_neighbor, file = "data/data_for_reg.RData")
#####
#加载之前数据
#####

load("data/data_for_reg.RData")

#####
#weitao:做回归之前先检查下自相关

#####
#未标准化的数据
#####
library(ggplot2)

# 选择数值型变量列
numeric_vars <- reg[, c("qr_AM", "qr_EM", "qr_BZ", "qr_Pn", "qr_Pq", "GX", "GY", "DBH1", "DBH2", "AD", "SRL", "SRA", "soc", "tn", "tp", "ap", "ph", "pd50_unweigh", "mpd50_unweigh", "mpd50_weigh", "mntd50_unweigh", "mntd50_weigh", "pd20_unweigh", "mpd20_unweigh", "mpd20_weigh", "mntd20_unweigh", "mntd20_weigh", "pd10_unweigh", "mpd10_unweigh", "mpd10_weigh", "mntd10_unweigh", "mntd10_weigh", "shannon_div_20", "invsimpson_div_20", "simpson_div_20", "shannon_div_10", "invsimpson_div_10", "simpson_div_10", "shannon_div_50", "invsimpson_div_50", "simpson_div_50", "BD_20", "CBD_20", "HBD_20", "BD_10", "CBD_10", "HBD_10", "BD_50", "CBD_50", "HBD_50", "growth_rate")]

# 提取数值列
numeric_vars <- as.matrix(numeric_vars)

# 计算相关性矩阵
cor_matrix <- cor(numeric_vars, use = "complete.obs")

# 初始化显著性矩阵
p_values <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))

# 计算显著性矩阵
for (i in 1:(ncol(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    cor_test <- cor.test(numeric_vars[, i], numeric_vars[, j], method = "pearson")
    p_values[i, j] <- cor_test$p.value
    p_values[j, i] <- cor_test$p.value
  }
}

# 将显著性矩阵命名为 p_values
colnames(p_values) <- colnames(cor_matrix)
rownames(p_values) <- colnames(cor_matrix)

library(ggplot2)

# 创建相关关系图数据框
high_corr_pairs <- which(cor_matrix > 0.5 & cor_matrix < 1, arr.ind = TRUE)
high_corr_vars <- rownames(cor_matrix)[high_corr_pairs[, 1]]
high_corr_pairs <- data.frame(
  Variable1 = high_corr_vars[high_corr_pairs[, 1]],
  Variable2 = colnames(cor_matrix)[high_corr_pairs[, 2]],
  Correlation = cor_matrix[high_corr_pairs],
  PValue = p_values[high_corr_pairs]
)

# 创建相关关系图
p <- ggplot(high_corr_pairs, aes(x = Variable1, y = Variable2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = ifelse(PValue < 0.05, "*", "")), color = "black", size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot with Significance Annotations (origin)")

print(p)

###为便于查看，把它制作成表格
# 创建一个空的数据框来存储表格内容
table_data <- data.frame(Variable1 = character(0), Variable2 = character(0), Correlation = numeric(0), PValue = numeric(0))

# 提取高相关性变量对及其相关性和显著性值
for (i in 1:nrow(high_corr_pairs)) {
  row <- high_corr_pairs[i, ]
  cor_val <- row$Correlation
  p_val <- row$PValue
  variable1 <- row$Variable1
  variable2 <- row$Variable2
  new_row <- data.frame(Variable1 = variable1, Variable2 = variable2, Correlation = cor_val, PValue = p_val)
  table_data <- rbind(table_data, new_row)
}

# 打印提取的表格内容
print(table_data)

#####
#这个相关关系可能不适合，换kendall的方法
#####

cor_matrix_kendall <- cor(numeric_vars, method = "kendall", use = "complete.obs")
p_values_kendall <- matrix(NA, nrow = ncol(cor_matrix_kendall), ncol = ncol(cor_matrix_kendall))
for (i in 1:(ncol(cor_matrix_kendall) - 1)) {
for (j in (i + 1):ncol(cor_matrix_kendall)) {
cor_test <- cor.test(numeric_vars[, i], numeric_vars[, j], method = "kendall")
p_values_kendall[i, j] <- cor_test$p.value
p_values_kendall[j, i] <- cor_test$p.value
}
}
# 将显著性矩阵命名为 p_values_kendall
colnames(p_values_kendall) <- colnames(cor_matrix_kendall)
rownames(p_values_kendall) <- colnames(cor_matrix_kendall)
library(ggplot2)


# 创建相关关系图数据框
high_corr_pairs <- which(cor_matrix_kendall > 0.5 & cor_matrix_kendall < 1, arr.ind = TRUE)
high_corr_vars <- rownames(cor_matrix_kendall)[high_corr_pairs[, 1]]
high_corr_pairs <- data.frame(
Variable1 = high_corr_vars[high_corr_pairs[, 1]],
Variable2 = colnames(cor_matrix_kendall)[high_corr_pairs[, 2]],
Correlation = cor_matrix_kendall[high_corr_pairs],
PValue = p_values_kendall[high_corr_pairs]
)
# 创建相关关系图
p <- ggplot(high_corr_pairs, aes(x = Variable1, y = Variable2, fill = Correlation)) +
geom_tile() +
geom_text(aes(label = ifelse(PValue < 0.05, "*", "")), color = "black", size = 6) +
scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "Correlation Plot with Significance Annotations (kendall)")
print(p)

#####
#标准化一下
#####

numeric_vars_zscore <- apply(numeric_vars, 2, function(x) (x - mean(x)) / sd(x))

# 提取数值列
numeric_vars_zscore <- as.matrix(numeric_vars_zscore)

# 计算相关性矩阵
cor_matrix <- cor(numeric_vars_zscore, method = "kendall", use = "complete.obs")

# 初始化显著性矩阵
p_values <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))

# 计算显著性矩阵
for (i in 1:(ncol(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    cor_test <- cor.test(numeric_vars_zscore[, i], numeric_vars_zscore[, j], method = "kendall")
    p_values[i, j] <- cor_test$p.value
    p_values[j, i] <- cor_test$p.value
  }
}

# 将显著性矩阵命名为 p_values
colnames(p_values) <- colnames(cor_matrix)
rownames(p_values) <- colnames(cor_matrix)

library(ggplot2)

# 创建相关关系图数据框
high_corr_pairs <- which(cor_matrix > 0.5 & cor_matrix < 1, arr.ind = TRUE)

high_corr_vars <- rownames(cor_matrix)[high_corr_pairs[, 1]]
high_corr_pairs <- data.frame(
  Variable1 = high_corr_vars[high_corr_pairs[, 1]],
  Variable2 = colnames(cor_matrix)[high_corr_pairs[, 2]],
  Correlation = cor_matrix[high_corr_pairs],
  PValue = p_values[high_corr_pairs]
)

# 创建相关关系图
p <- ggplot(high_corr_pairs, aes(x = Variable1, y = Variable2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = ifelse(PValue < 0.05, "*", "")), color = "black", size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Plot with Significance Annotations(zscore)")

print(p)

###为便于查看，把它制作成表格
# 创建一个空的数据框来存储表格内容
table_data <- data.frame(Variable1 = character(0), Variable2 = character(0), Correlation = numeric(0), PValue = numeric(0))

# 提取高相关性变量对及其相关性和显著性值
for (i in 1:nrow(high_corr_pairs)) {
  row <- high_corr_pairs[i, ]
  cor_val <- row$Correlation
  p_val <- row$PValue
  variable1 <- row$Variable1
  variable2 <- row$Variable2
  new_row <- data.frame(Variable1 = variable1, Variable2 = variable2, Correlation = cor_val, PValue = p_val)
  table_data <- rbind(table_data, new_row)
}

# 打印提取的表格内容
print(table_data)
#β回归
library(tidyverse)
#关于AM的部分
>>>>>>> 26443777ed9086be7130a97033f3852e9444c867
beita_reg <- lm(qr_AM ~ shannon_div_20 + SRL + AD + mpd20_weigh + soc + tn + tp + ap + ph+ DBH2, data = reg)
step_beita_reg <- step(beita_reg)
summary(beita_reg)
summary(step_beita_reg)
