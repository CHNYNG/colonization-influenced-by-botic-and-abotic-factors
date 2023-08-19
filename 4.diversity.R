# 数据准备
#加载数据
#load("E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl_soil.RData")
####改造一下成为让韦韬师兄不嫌弃的代码，呜呜呜
#首先是数据，用和phylodiversity的一样就可以
#包含 sp_loc（所有focal物种的位置）、hsd_alive_singlebr黑石顶的物种数据


#####
#其他多样性好麻烦哦，还是用vegan包吧。。。

####我的物种在specieslist里
unique_species <- gsub(" ", "_",unique_species)
hsd_vegan_alpha_20 <- data.frame(scientific.name = unique_species)

#一共有unique_species294个，把它们作为第一列
for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 20)
  
  species_counts <- table(hsd_sub$scientific.name) %>%
    as.data.frame() %>%
    rename(scientific.name = Var1, Freq = Freq)
  
  hsd_vegan <- data.frame(scientific.name = unique_species)
  
  hsd_vegan <- hsd_vegan %>%
    left_join(.,species_counts, by = "scientific.name") %>%
    mutate(Freq = ifelse(is.na(Freq), 0, Freq)) 
    
  hsd_vegan_alpha_20[,i+1] <- hsd_vegan$Freq
    
}
#整理一下下
library(tidyr)
rownames(hsd_vegan_alpha_20) <- hsd_vegan_alpha_20$scientific.name
hsd_vegan_alpha_20 <- hsd_vegan_alpha_20[,-1]
colnames(hsd_vegan_alpha_20) <- sp_loc$TagNew
hsd_vegan_alpha_20 <- t(hsd_vegan_alpha_20)

library(vegan)
#shannon
shannon_div_20 <- diversity(hsd_vegan_alpha_20, index = "shannon")
shannon_div_20 <- as.data.frame(shannon_div_20)
#皮尔逊多样性指数（Pielou's Evenness Index）
invsimpson_div_20 <- diversity(hsd_vegan_alpha_20, index = "invsimpson")
invsimpson_div_20 <- as.data.frame(invsimpson_div_20)
#辛普森多样性指数（Simpson Diversity Index）
simpson_div_20 <- diversity(hsd_vegan_alpha_20, index = "simpson")
simpson_div_20 <- as.data.frame(simpson_div_20)
######整理数据为β多样性矩阵(需要的话)
hsd_vegan_beita_20 <- ifelse(hsd_vegan_alpha_20 >= 1, 1, 0)
#计算β多样性指数
# 计算 Jaccard 相似性指数作为β多样性的距离
#beta_div_jaccard_20 <- vegdist(hsd_vegan_alpha_20, method = "jaccard")

# 计算 Bray-Curtis 相似性指数作为β多样性的距离
#beta_div_bray_20 <- vegdist(hsd_vegan_alpha_20, method = "bray")

# 计算 Sørensen-Dice 相似性指数作为β多样性的距离
#beta_div_sorensen_20 <- vegdist(hsd_vegan_alpha_20, method = "kulczynski")

# 计算 Euclidean 距离作为β多样性的距离
#beta_div_euclidean_20 <- vegdist(hsd_vegan_alpha_20, method = "euclidean")

#####
#整理一下
sp_loc <- sp_loc %>% 
  mutate(
    shannon_div_20 = shannon_div_20$shannon_div_20,
    invsimpson_div_20 = invsimpson_div_20$invsimpson_div_20,
    simpson_div_20 = simpson_div_20$simpson_div_20
  )

#删掉一些无用的变量
rm(species_counts,hsd_vegan,hsd_sub,i,invsimpson_div_20,
   shannon_div_20,simpson_div_20,hsd_vegan_alpha_20,
   hsd_vegan_beita_20)

#####
#算一下10m的距离
unique_species <- gsub(" ", "_",unique_species)
hsd_vegan_alpha_10 <- data.frame(scientific.name = unique_species)

#一共有unique_species294个，把它们作为第一列
for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 10)
  
  species_counts <- table(hsd_sub$scientific.name) %>%
    as.data.frame() %>%
    rename(scientific.name = Var1, Freq = Freq)
  
  hsd_vegan <- data.frame(scientific.name = unique_species)
  
  hsd_vegan <- hsd_vegan %>%
    left_join(.,species_counts, by = "scientific.name") %>%
    mutate(Freq = ifelse(is.na(Freq), 0, Freq)) 
  
  hsd_vegan_alpha_10[,i+1] <- hsd_vegan$Freq
  
}
#整理一下下
library(tidyr)
rownames(hsd_vegan_alpha_10) <- hsd_vegan_alpha_10$scientific.name
hsd_vegan_alpha_10 <- hsd_vegan_alpha_10[,-1]
colnames(hsd_vegan_alpha_10) <- sp_loc$TagNew
hsd_vegan_alpha_10 <- t(hsd_vegan_alpha_10)

library(vegan)
#shannon
shannon_div_10 <- diversity(hsd_vegan_alpha_10, index = "shannon")
shannon_div_10 <- as.data.frame(shannon_div_10)
#皮尔逊多样性指数（Pielou's Evenness Index）
invsimpson_div_10 <- diversity(hsd_vegan_alpha_10, index = "invsimpson")
invsimpson_div_10 <- as.data.frame(invsimpson_div_10)
#辛普森多样性指数（Simpson Diversity Index）
simpson_div_10 <- diversity(hsd_vegan_alpha_10, index = "simpson")
simpson_div_10 <- as.data.frame(simpson_div_10)
######整理数据为β多样性矩阵(需要的话)
hsd_vegan_beita_10 <- ifelse(hsd_vegan_alpha_10 >= 1, 1, 0)
#计算β多样性指数
# 计算 Jaccard 相似性指数作为β多样性的距离
#beta_div_jaccard_10 <- vegdist(hsd_vegan_alpha_10, method = "jaccard")

# 计算 Bray-Curtis 相似性指数作为β多样性的距离
#beta_div_bray_10 <- vegdist(hsd_vegan_alpha_10, method = "bray")

# 计算 Sørensen-Dice 相似性指数作为β多样性的距离
#beta_div_sorensen_10 <- vegdist(hsd_vegan_alpha_10, method = "kulczynski")

# 计算 Euclidean 距离作为β多样性的距离
#beta_div_euclidean_10 <- vegdist(hsd_vegan_alpha_10, method = "euclidean")

#####
#整理一下
sp_loc <- sp_loc %>% 
  mutate(
    shannon_div_10 = shannon_div_10$shannon_div_10,
    invsimpson_div_10 = invsimpson_div_10$invsimpson_div_10,
    simpson_div_10 = simpson_div_10$simpson_div_10
  )

#删掉一些无用的变量
rm(species_counts,hsd_vegan,hsd_sub,i,invsimpson_div_10,
   shannon_div_10,simpson_div_10,hsd_vegan_alpha_10,
   hsd_vegan_beita_10)

#####
#算一下10m的距离
unique_species <- gsub(" ", "_",unique_species)
hsd_vegan_alpha_50 <- data.frame(scientific.name = unique_species)

#一共有unique_species294个，把它们作为第一列
for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 50)
  
  species_counts <- table(hsd_sub$scientific.name) %>%
    as.data.frame() %>%
    rename(scientific.name = Var1, Freq = Freq)
  
  hsd_vegan <- data.frame(scientific.name = unique_species)
  
  hsd_vegan <- hsd_vegan %>%
    left_join(.,species_counts, by = "scientific.name") %>%
    mutate(Freq = ifelse(is.na(Freq), 0, Freq)) 
  
  hsd_vegan_alpha_50[,i+1] <- hsd_vegan$Freq
  
}
#整理一下下
library(tidyr)
rownames(hsd_vegan_alpha_50) <- hsd_vegan_alpha_50$scientific.name
hsd_vegan_alpha_50 <- hsd_vegan_alpha_50[,-1]
colnames(hsd_vegan_alpha_50) <- sp_loc$TagNew
hsd_vegan_alpha_50 <- t(hsd_vegan_alpha_50)

library(vegan)
#shannon
shannon_div_50 <- diversity(hsd_vegan_alpha_50, index = "shannon")
shannon_div_50 <- as.data.frame(shannon_div_50)
#皮尔逊多样性指数（Pielou's Evenness Index）
invsimpson_div_50 <- diversity(hsd_vegan_alpha_50, index = "invsimpson")
invsimpson_div_50 <- as.data.frame(invsimpson_div_50)
#辛普森多样性指数（Simpson Diversity Index）
simpson_div_50 <- diversity(hsd_vegan_alpha_50, index = "simpson")
simpson_div_50 <- as.data.frame(simpson_div_50)
######整理数据为β多样性矩阵(需要的话)
hsd_vegan_beita_50 <- ifelse(hsd_vegan_alpha_50 >= 1, 1, 0)
#计算β多样性指数
# 计算 Jaccard 相似性指数作为β多样性的距离
#beta_div_jaccard_50 <- vegdist(hsd_vegan_alpha_50, method = "jaccard")

# 计算 Bray-Curtis 相似性指数作为β多样性的距离
#beta_div_bray_50 <- vegdist(hsd_vegan_alpha_50, method = "bray")

# 计算 Sørensen-Dice 相似性指数作为β多样性的距离
#beta_div_sorensen_50 <- vegdist(hsd_vegan_alpha_50, method = "kulczynski")

# 计算 Euclidean 距离作为β多样性的距离
#beta_div_euclidean_50 <- vegdist(hsd_vegan_alpha_50, method = "euclidean")

#####
#整理一下
sp_loc <- sp_loc %>% 
  mutate(
    shannon_div_50 = shannon_div_50$shannon_div_50,
    invsimpson_div_50 = invsimpson_div_50$invsimpson_div_50,
    simpson_div_50 = simpson_div_50$simpson_div_50
  )

#删掉一些无用的变量
rm(species_counts,hsd_vegan,hsd_sub,i,invsimpson_div_50,
   shannon_div_50,simpson_div_50,hsd_vegan_alpha_50,
   hsd_vegan_beita_50)