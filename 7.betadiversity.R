#beta diversity 思路

#### 首先是需要所有的β多样性先做20

######整理数据为β多样性矩阵(需要的话)
hsd_vegan_beita_20 <- ifelse(hsd_vegan_alpha_20 >= 0.5, 1, 0)
hsd_vegan_beita_20 <- as.matrix(hsd_vegan_beita_20)



#计算β多样性指数
#####
# 计算 Jaccard 相似性指数作为β多样性的距离
#####
beta_div_jaccard_20 <- vegdist(hsd_vegan_alpha_20, method = "jaccard")
beta_div_jaccard_20_d <- as.data.frame(as.matrix(beta_div_jaccard_20))
# 去掉列名中的 "X" 前缀
colnames(beta_div_jaccard_20_d) <- gsub("^X", "", colnames(beta_div_jaccard_20_d))
# 将行名和列名添加为新列
beta_div_jaccard_20_d$RowName <- rownames(beta_div_jaccard_20_d)
colnames(beta_div_jaccard_20_d) <- make.names(colnames(beta_div_jaccard_20_d))  # 使列名合法

# 将数据框重塑为长格式
library(reshape2)
beta_div_jaccard_20 <- melt(beta_div_jaccard_20_d, id.vars = "RowName")

# 修改列名
colnames(beta_div_jaccard_20) <- c("RowName", "ColName", "jaccard_20")
beta_div_jaccard_20$ColName <- gsub("^X", "", beta_div_jaccard_20$ColName)
# 删除Value为0的行
beta_div_jaccard_20 <- beta_div_jaccard_20 %>%
  filter(jaccard_20 != 0)

beta_div_jaccard_20 <- beta_div_jaccard_20 %>%
  mutate(RowName = pmin(RowName, ColName), ColName = pmax(RowName, ColName)) %>%
  distinct(RowName, ColName, .keep_all = TRUE)


#####
# 计算 Bray-Curtis 相似性指数作为β多样性的距离
#####
beta_div_bray_20 <- vegdist(hsd_vegan_alpha_20, method = "bray")
beta_div_bray_20_d <- as.data.frame(as.matrix(beta_div_bray_20))
# 去掉列名中的 "X" 前缀
colnames(beta_div_bray_20_d) <- gsub("^X", "", colnames(beta_div_bray_20_d))
# 将行名和列名添加为新列
beta_div_bray_20_d$RowName <- rownames(beta_div_bray_20_d)
colnames(beta_div_bray_20_d) <- make.names(colnames(beta_div_bray_20_d))  # 使列名合法

# 将数据框重塑为长格式
library(reshape2)
beta_div_bray_20 <- melt(beta_div_bray_20_d, id.vars = "RowName")

# 修改列名
colnames(beta_div_bray_20) <- c("RowName", "ColName", "bray_20")
beta_div_bray_20$ColName <- gsub("^X", "", beta_div_bray_20$ColName)
# 删除Value为0的行
beta_div_bray_20 <- beta_div_bray_20 %>%
  filter(bray_20 != 0)

beta_div_bray_20 <- beta_div_bray_20 %>%
  mutate(RowName = pmin(RowName, ColName), ColName = pmax(RowName, ColName)) %>%
  distinct(RowName, ColName, .keep_all = TRUE)
###嗷嗷嗷嗷嗷嗷！！！终于好了！！！！我好菜…


#####
# 计算 Sørensen-Dice 相似性指数作为β多样性的距离
#####
beta_div_sorensen_20 <- vegdist(hsd_vegan_alpha_20, method = "kulczynski")
beta_div_sorensen_20_d <- as.data.frame(as.matrix(beta_div_sorensen_20))
# 去掉列名中的 "X" 前缀
colnames(beta_div_sorensen_20_d) <- gsub("^X", "", colnames(beta_div_sorensen_20_d))
# 将行名和列名添加为新列
beta_div_sorensen_20_d$RowName <- rownames(beta_div_sorensen_20_d)
colnames(beta_div_sorensen_20_d) <- make.names(colnames(beta_div_sorensen_20_d))  # 使列名合法

# 将数据框重塑为长格式
library(reshape2)
beta_div_sorensen_20 <- melt(beta_div_sorensen_20_d, id.vars = "RowName")

# 修改列名
colnames(beta_div_sorensen_20) <- c("RowName", "ColName", "sorensen_20")
beta_div_sorensen_20$ColName <- gsub("^X", "", beta_div_sorensen_20$ColName)
# 删除Value为0的行
beta_div_sorensen_20 <- beta_div_sorensen_20 %>%
  filter(sorensen_20 != 0)

beta_div_sorensen_20 <- beta_div_sorensen_20 %>%
  mutate(RowName = pmin(RowName, ColName), ColName = pmax(RowName, ColName)) %>%
  distinct(RowName, ColName, .keep_all = TRUE)


#####
# 计算 Euclidean 距离作为β多样性的距离
#####
beta_div_euclidean_20 <- vegdist(hsd_vegan_alpha_20, method = "euclidean")
beta_div_euclidean_20_d<- as.data.frame(as.matrix(beta_div_euclidean_20))
# 去掉列名中的 "X" 前缀
colnames(beta_div_euclidean_20_d) <- gsub("^X", "", colnames(beta_div_euclidean_20_d))
# 将行名和列名添加为新列
beta_div_euclidean_20_d$RowName <- rownames(beta_div_euclidean_20_d)
colnames(beta_div_euclidean_20_d) <- make.names(colnames(beta_div_euclidean_20_d))  # 使列名合法

# 将数据框重塑为长格式
library(reshape2)
beta_div_euclidean_20 <- melt(beta_div_euclidean_20_d, id.vars = "RowName")

# 修改列名
colnames(beta_div_euclidean_20) <- c("RowName", "ColName", "euclidean_20")
beta_div_euclidean_20$ColName <- gsub("^X", "", beta_div_euclidean_20$ColName)
# 删除Value为0的行
beta_div_euclidean_20 <- beta_div_euclidean_20 %>%
  filter(euclidean_20 != 0)

beta_div_euclidean_20 <- beta_div_euclidean_20 %>%
  mutate(RowName = pmin(RowName, ColName), ColName = pmax(RowName, ColName)) %>%
  distinct(RowName, ColName, .keep_all = TRUE)

#####
#整理一下β多样性
#####
beta_div <- data.frame(A = beta_div_bray_20$RowName,
                       B = beta_div_bray_20$ColName,
                       bray_20 = beta_div_bray_20$bray_20,
                       euclidean_20 = beta_div_euclidean_20$euclidean_20,
                       jaccard_20 = beta_div_jaccard_20$jaccard_20,
                       sorensen_20 = beta_div_sorensen_20$sorensen_20
                       )

focal_sp <- root_qrl %>%
  rename(scientific.name = Latin) %>%
  select(scientific.name, TagNew)

beta_div <- beta_div %>%
  left_join(focal_sp, by = c("A" = "TagNew")) %>%
  select(A, B, bray_20, euclidean_20, jaccard_20, sorensen_20, A.scientific.name = scientific.name)
beta_div <- beta_div %>%
  left_join(focal_sp, by = c("B" = "TagNew")) %>%
  select(A, B, bray_20, euclidean_20, jaccard_20, sorensen_20, A.scientific.name, B.scientific.name = scientific.name)

#把beta多样性和侵染的插值合并
beta_div <- beta_div %>%
  left_join(root_qrl, by = c("A" = "TagNew")) %>%
  select(A, B, bray_20, euclidean_20, jaccard_20, sorensen_20, A.scientific.name, B.scientific.name,
         A.qr = qr_AM)
beta_div <- beta_div %>%
  left_join(root_qrl, by = c("B" = "TagNew")) %>%
  select(A, B, bray_20, euclidean_20, jaccard_20, sorensen_20, A.scientific.name, B.scientific.name,
         A.qr, B.qr = qr_AM)
beta_div <- beta_div %>%
  mutate(de_qr = abs(A.qr - B.qr)) %>%
  mutate(de_qr = ifelse(is.na(de_qr), 0, de_qr)) %>%
  select(A, B, bray_20, euclidean_20, jaccard_20, sorensen_20, A.scientific.name, B.scientific.name,
         de_qr)
#####
#regression
#####
#这里qr还是要用beta回归
library(betareg)
beta_div$de_qr <- ifelse(beta_div$de_qr <= 0, 0.01, ifelse(beta_div$de_qr >= 1, 0.99, beta_div$de_qr))
beta_model <- betareg(de_qr ~ bray_20 , data = beta_div, link = "probit")
