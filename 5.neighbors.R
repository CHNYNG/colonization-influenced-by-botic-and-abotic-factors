#####
###计算邻体的DBH和距离的加和
#数据集加一列

#加载前步骤
#source("4.diversity.R")

#####
#下面为不分同种异种的BD
#创建一个dataframe用来存储

BD_20 <- c()
CBD_20 <- c()
HBD_20 <- c()


  
for (i in 1:dim(sp_loc)[1]) {
  
    # 20m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 20)
  hsd_sub <- hsd_sub %>%
    mutate(basal_area = pi * (dbh_multi / 2) ^ 2,
           ba_distance_weighted = ifelse(distance == 0 , 0, basal_area / distance))
  
  #不分同种异种sum之后再相除，影响可能是二次方的关系
  BD_20[i] <- sum(hsd_sub$ba_distance_weighted)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name) %>%
    unique()
  
  if (length(con_species) == 0) {
    # 处理con_species为空的情况，例如将CBD_20[i]设置为0
    CBD_20[i] <- 0
  } else {
    hsd_sub_con <- hsd_sub %>%
      filter(scientific.name %in% con_species)  # 使用%in%来匹配多个scientific.name
    
    CBD_20[i] <- sum(hsd_sub_con$ba_distance_weighted)
  }
  
   #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  HBD_20[i] <- sum(hsd_sub_het$ba_distance_weighted)
}

###下面进行10和50的步骤

BD_10 <- c()
CBD_10 <- c()
HBD_10 <- c()



for (i in 1:dim(sp_loc)[1]) {
  
  # 10m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 10)
  hsd_sub <- hsd_sub %>%
    mutate(basal_area = pi * (dbh_multi / 2) ^ 2,
           ba_distance_weighted = ifelse(distance == 0 , 0, basal_area / distance))
  
  #不分同种异种sum之后再相除，影响可能是二次方的关系
  BD_10[i] <- sum(hsd_sub$ba_distance_weighted)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name) %>%
    unique()
  
  if (length(con_species) == 0) {
    # 处理con_species为空的情况，例如将CBD_10[i]设置为0
    CBD_10[i] <- 0
  } else {
    hsd_sub_con <- hsd_sub %>%
      filter(scientific.name %in% con_species)  # 使用%in%来匹配多个scientific.name
    
    CBD_10[i] <- sum(hsd_sub_con$ba_distance_weighted)
  }
  
  #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  HBD_10[i] <- sum(hsd_sub_het$ba_distance_weighted)
}


BD_50 <- c()
CBD_50 <- c()
HBD_50 <- c()

for (i in 1:dim(sp_loc)[1]) {
  
  # 50m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 50)
  hsd_sub <- hsd_sub %>%
    mutate(basal_area = pi * (dbh_multi / 2) ^ 2,
           ba_distance_weighted = ifelse(distance == 0 , 0, basal_area / distance))
  
  #不分同种异种sum之后再相除，影响可能是二次方的关系
  BD_50[i] <- sum(hsd_sub$ba_distance_weighted)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name) %>%
    unique()
  
  if (length(con_species) == 0) {
    # 处理con_species为空的情况，例如将CBD_50[i]设置为0
    CBD_50[i] <- 0
  } else {
    hsd_sub_con <- hsd_sub %>%
      filter(scientific.name %in% con_species)  # 使用%in%来匹配多个scientific.name
    
    CBD_50[i] <- sum(hsd_sub_con$ba_distance_weighted)
  }
  
  #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  HBD_50[i] <- sum(hsd_sub_het$ba_distance_weighted)
}

hsd_neighbor <- data.frame(sp_loc$TagNew)
colnames(hsd_neighbor) <- "TagNew"

hsd_neighbor <- hsd_neighbor %>%
  mutate(BD_20 =BD_20,
         CBD_20=CBD_20,
         HBD_20=HBD_20,
         BD_50 =BD_50,
         CBD_50=CBD_50,
         HBD_50=HBD_50,
         BD_10 =BD_10,
         CBD_10=CBD_10,
         HBD_10=HBD_10)
rm(BD_20,CBD_20,HBD_20,BD_50,CBD_50,HBD_50,BD_10,CBD_10,HBD_10,
   hsd_sub,hsd_sub_con,hsd_sub_het)
