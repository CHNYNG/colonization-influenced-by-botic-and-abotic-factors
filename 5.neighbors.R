#####
###计算邻体的DBH和距离的加和
#数据集加一列

#加载前步骤
source("4.diversity.R")

#####
#下面为不分同种异种的BD
#创建一个dataframe用来存储

hsd_neighbor <- data.frame(sp_loc$TagNew)
colnames(hsd_neighbor) <- "TagNew"
  
for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 20)
  
  #不分同种异种
  hsd_neighbor$BD_20[i] <- sum(hsd_sub$dbh_multi)/sum(hsd_sub$distance)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name)%>%
    unique()
  
  hsd_sub_con <- hsd_sub %>%
    filter(scientific.name == con_species)
  
    hsd_neighbor$CBD_20[i] <- sum(hsd_sub_con$dbh_multi) / sum(hsd_sub_con$distance)


   #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  hsd_neighbor$HBD_20[i] <- sum(hsd_sub_het$dbh_multi)/sum(hsd_sub_het$distance)
}

###下面进行10和50的步骤

for (i in 1:dim(sp_loc)[1]) {
  
  # 10m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 10)
  
  #不分同种异种
  hsd_neighbor$BD_10[i] <- sum(hsd_sub$dbh_multi)/sum(hsd_sub$distance)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name)%>%
    unique()
  
  hsd_sub_con <- hsd_sub %>%
    filter(scientific.name == con_species)
  
  hsd_neighbor$CBD_10[i] <- sum(hsd_sub_con$dbh_multi) / sum(hsd_sub_con$distance)
  
  
  #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  hsd_neighbor$HBD_10[i] <- sum(hsd_sub_het$dbh_multi)/sum(hsd_sub_het$distance)
}

for (i in 1:dim(sp_loc)[1]) {
  
  # 50m circle
  hsd_sub <- hsd_alive_singlebr %>%
    mutate(
      distance = sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2)
    ) %>%
    filter(distance <= 50)
  
  #不分同种异种
  hsd_neighbor$BD_50[i] <- sum(hsd_sub$dbh_multi)/sum(hsd_sub$distance)
  
  #同种
  con_species <- hsd_sub %>%
    filter(TagNew == sp_loc$TagNew[i]) %>%
    pull(scientific.name)%>%
    unique()
  
  hsd_sub_con <- hsd_sub %>%
    filter(scientific.name == con_species)
  
  hsd_neighbor$CBD_50[i] <- sum(hsd_sub_con$dbh_multi) / sum(hsd_sub_con$distance)
  
  
  #异种
  het_species <- hsd_sub %>%
    filter(TagNew != sp_loc$TagNew[i]) %>%
    pull(scientific.name)
  
  hsd_sub_het <- hsd_sub %>%
    filter(scientific.name %in% het_species)
  
  hsd_neighbor$HBD_50[i] <- sum(hsd_sub_het$dbh_multi)/sum(hsd_sub_het$distance)
}

hsd_neighbor <- hsd_neighbor %>%
  mutate_all(~ ifelse(is.infinite(.), 0, .))
