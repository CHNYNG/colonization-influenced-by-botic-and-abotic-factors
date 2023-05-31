# neighborselect 数据准备
#加载数据
#load("E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl_soil.RData")
HSD <- HSD_data_0[,c("TagNew","Latin","Qudrat","Species","GX","GY","Status1","Status2","DBH1","DBH2","H","Family","Genus","Species.x","abundance")]
HSD_species <- subset(HSD, select = c(Family, Genus, Species.x, Species, Latin))
HSD_species <- unique(HSD_species)
##1. 筛选出neighbors，这里使用distance＜10
# 创建一个空的数据框用于存储结果
result <- data.frame()

# 遍历d数据集的每一行
for (i in 1:nrow(d)) {
  # 计算当前行的点与HSD数据集中所有点的距离
  distances <- sqrt((d$GX[i] - HSD$GX)^2 + (d$GY[i] - HSD$GY)^2)
  
  # 筛选出距离小于10米的点的索引
  indices <- which(distances < 10 & distances > 0.001)
  
  # 如果存在距离小于10米的点，则将这些点添加到结果数据框中
  if (length(indices) > 0) {
    result <- rbind(result, data.frame(d[i,c("TagNew") ], HSD[indices,c("TagNew") ]))
  }
}

# 改一下列名，把focal的列名改为TagNew_f它的邻体为n
colnames(result) <- c("TagNew_f","TagNew_n")

# 给TagNew_f和TagNew_n带上它们原本的数据
library(dplyr)

neighborselect <- left_join(result, HSD, by = c("TagNew_f" = "TagNew"))

neighborselect <- left_join(neighborselect, HSD, by = c("TagNew_n" = "TagNew"))

# 整理一下列名
colnames(neighborselect) <- c("TagNew_f", "TagNew_n", "Latin_f", "Qudrat_F", "Species_f","GX_f",  "GY_f",  "Status1_f","Status2_F", "DBH1_f","DBH2_f", "H_f","Family_f",  "Genus_f",   "Species.x_f","abundance_f",
                              "Latin_n",   "Qudrat_n",  "Species_n", "GX_n", "GY_n", "Status1_n", "Status2_n", "DBH1_n", "DBH2_n","H_n",   "Family_n",  "Genus_n", "Species.x_n", "abundance_n")

# 3. 把距离打出来
neighborselect$distance <- sqrt((neighborselect$GX_f-neighborselect$GX_n)^2+(neighborselect$GY_f-neighborselect$GY_n)^2)

### 给它们分个家，把同种和异种分出来
neighborselect$Type <- ifelse(neighborselect$Species_f == neighborselect$Species_n,
                              "conspecific", "heterospecific")


### 2. 计算一下每一棵邻居树的截面积(cm²)。基面积是指树木横截面（通常是树干截面）的面积，可以通过测量直径或周长，并应用相应的公式计算得到。
neighborselect$ba1_n <- pi * (neighborselect$DBH1_n/2)^2
neighborselect$ba2_n <- pi * (neighborselect$DBH2_n/2)^2

### 4.5. 对于每棵邻居树，计算其逆距离加权基面积。将邻居树的基面积乘以其对应的逆距离权重，得到加权的基面积。
neighborselect$ba2_dis_n <- neighborselect$ba2_n/(neighborselect$distance+0.001)

### 6. 将所有邻居树的逆距离加权基面积相加，得到种内密度的总和和种间密度的总和

library(dplyr)

heterospecific_dis <- neighborselect %>%
  filter(Type == "heterospecific") %>%
  group_by(TagNew_f) %>%
  summarise(sum_ba2_dis_n = sum(ba2_dis_n, na.rm = TRUE)) %>%
  ungroup()
colnames(heterospecific_dis) <- c("TagNew_f","heterospecific_dis")

conspecific_dis <- neighborselect %>%
  filter(Type == "conspecific") %>%
  group_by(TagNew_f) %>%
  summarise(sum_ba2_dis_n = sum(ba2_dis_n, na.rm = TRUE)) %>%
  ungroup()
colnames(conspecific_dis) <- c("TagNew_f","conspecific_dis")

library(dplyr)

d <- left_join(d, conspecific_dis, by = c("TagNew" = "TagNew_f"))
d <- left_join(d, heterospecific_dis, by = c("TagNew" = "TagNew_f"))
























