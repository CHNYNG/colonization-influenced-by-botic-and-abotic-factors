# 系统发育
library(picante)
library(tidyverse)

#数据集是root_qrl，黑石顶的数据集是HSD_data
#这里需要把HSD的数据集match上科属种
library(dplyr)
HSD_data <- left_join(HSD_data,specieslist,
                      by = "Latin")
# 删除Family列中值为NA的行
HSD_data <- HSD_data[!is.na(HSD_data$Family), ]

##导入系统发育树
hsd_phytr <- read.tree("Data/hsd_tree_vphylo.rwk")
alive_status <- c("Alive", "P", "lean", "alive", "snag")
hsd_alive <- HSD_data %>%
  filter(Status2 %in% alive_status, 
         !is.na(TagNew), 
         !is.na(Branch), 
         !is.na(DBH2))

#对hsd_alive数据框按照tagnew列的不同取值进行分组，
#然后在每个分组内计算dbh2列的平方和，
#再取其平方根值，最后得到一个新的数据框，
#其中包含了每个分组的tagnew和对应的dbh_multi值
hsd_multibr <- hsd_alive %>%
  group_by(TagNew) %>%
  summarise(dbh_multi = sqrt(sum(DBH2 ^ 2)))


#从hsd_alive数据框中筛选出branch列中值为0的行，
#然后将该结果与hsd_multibr数据框按照tagnew列进行左连接。
#最后，对连接后的数据框中的scientific.name列进行修改，
#将其中的空格替换为下划线，并得到一个
#新的数据框hsd_alive_singlebr

hsd_alive_singlebr <- hsd_alive %>% 
  filter(Branch == 0) %>% 
  left_join(hsd_multibr, by = "TagNew") %>% 
  mutate(
    scientific.name = str_replace(scientific.name, " ", "_")
  )

##一个看不懂的循环
# loops ----
pd20_unweigh <- c();
mpd20_unweigh <- c(); mpd20_weigh <- c()
mntd20_unweigh <- c(); mntd20_weigh <- c()

for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((gx - sp_loc$gx[i])^2 + (gy - sp_loc$gy[i])^2) <= 20)
  
  # phylo tree with species in 20m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  # abundance-un-weighted phylo 
  # pd
  pd20_unweigh[i] <- sum(tr_sub$edge.length)
  
  # mpd
  tr_sub_dis <- cophenetic(tr_sub)
  tr_sub_dis <- tr_sub_dis[lower.tri(tr_sub_dis, diag = FALSE)]
  mpd20_unweigh[i] <- mean(tr_sub_dis)
  
  # mntd
  tr_sub_dis2 <- cophenetic(tr_sub)
  diag(tr_sub_dis2) <- NA
  mntd20_unweigh[i] <- mean(apply(tr_sub_dis2, 2, min, na.rm = TRUE))
  
  # abundance-weighted phylo
  hsd_sub_sub <- hsd_sub %>% 
    filter(scientific.name %in% tr_sub$tip.label)
  
  tr_sub_community <-
    matrix(table(hsd_sub_sub$scientific.name), 
           nrow = 1, 
           dimnames = list("s1", names(table(hsd_sub_sub$scientific.name))))
  # mpd
  mpd20_weigh[i] <- picante::mpd(tr_sub_community, cophenetic(tr_sub),
                                 abundance.weighted = T)
  # mntd
  mntd20_weigh[i] <- picante::mntd(tr_sub_community, cophenetic(tr_sub),
                                   abundance.weighted = T)
}
