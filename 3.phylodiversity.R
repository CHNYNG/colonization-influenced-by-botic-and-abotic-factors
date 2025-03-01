#####
#加载前步骤
#source("2.krige.R")

# 系统发育

library(picante)
library(tidyverse)

#数据集是root_qrl，黑石顶的数据集是HSD_data
#location of samples

filtered_root_qrl <- root_qrl %>% 
  filter(!is.na(GX),
         Branch ==0 )
filtered_root_qrl <- distinct(filtered_root_qrl)

sp_loc <- data.frame(
  GX = as.numeric(filtered_root_qrl$GX),
  GY = filtered_root_qrl$GY,
  TagNew = filtered_root_qrl$TagNew
)
sp_loc <- distinct(sp_loc, TagNew, .keep_all = TRUE)


##画系统发育树
#library(V.PhyloMaker)
#colnames(specieslist) <- c("species", "Genus", "Family", "Family_number", "Order", "Group")
#hsd_phy <- phylo.maker(sp.list = specieslist)
#hsd_phytr <- hsd_phy$scenario.3
#保存系统发育树
#write.tree(hsd_phytr, "data/hsd_tree_vphylo.rwk")
colnames(specieslist) <- c("Latin", "Genus", "Family", "Family_number", "Order", "Group")
##导入系统发育树
hsd_phytr <- read.tree("data/hsd_tree_vphylo.rwk")

#test <- specieslist==test1
#test1 <- gsub("_", " ",hsd_phytr$tip.label)
#test <- setdiff(specieslist$Latin, test1)
#这里需要把HSD的数据集match上科属种
library(dplyr)
HSD_data <- left_join(HSD_data,specieslist,
                      by = "Latin")
# 删除Family列中值为NA的行
#不能删掉这个，这样会减少物种
#已经核对出减少了
#水杨梅、香桂、菴耳柯、檵木、钩藤、细轴荛花
#往前回溯一下specieslist
#specieslist里都有
#HSD_data <- HSD_data[!is.na(HSD_data$Family), ]

alive_status <- c("Alive", "P", "lean", "alive", "snag")
hsd_alive <- HSD_data %>%
  filter(Status2 %in% alive_status, 
         !is.na(TagNew), 
         !is.na(Branch), 
         !is.na(DBH2)) %>%
  rename(scientific.name = Latin)


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
    scientific.name = gsub(" ", "_",scientific.name),
    GX = as.numeric(as.character(GX))
  )


# 将TagNew列的值尝试转换为数值
hsd_alive_singlebr$TagNew <- as.numeric(hsd_alive_singlebr$TagNew)


# 过滤掉NA值和为0的值
hsd_alive_singlebr <- hsd_alive_singlebr %>%
  filter(!is.na(TagNew)) %>%
  filter(TagNew != 0) %>%
  filter(floor(TagNew) == TagNew)%>%
  distinct(TagNew,.keep_all = TRUE)

hsd_alive_singlebr$TagNew <- as.character(hsd_alive_singlebr$TagNew)

library(stringr)
hsd_alive_singlebr <- hsd_alive_singlebr %>%
  mutate(TagNew = str_pad(TagNew, width = 7, side = "left", pad = "0"))

##一个看不懂的循环
# loops ----

pd100_unweigh <- c();
mpd100_unweigh <- c(); mpd100_weigh <- c()
mntd100_unweigh <- c(); mntd100_weigh <- c()

pd50_unweigh <- c();
mpd50_unweigh <- c(); mpd50_weigh <- c()
mntd50_unweigh <- c(); mntd50_weigh <- c()

pd20_unweigh <- c();
mpd20_unweigh <- c(); mpd20_weigh <- c()
mntd20_unweigh <- c(); mntd20_weigh <- c()

pd10_unweigh <- c();
mpd10_unweigh <- c(); mpd10_weigh <- c()
mntd10_unweigh <- c(); mntd10_weigh <- c()

for (i in 1:dim(sp_loc)[1]) {
  
  #100m circle
  
  hsd_s <- hsd_alive_singlebr %>%
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 100)
                
   # phylo tree with species in 100m circle
   tr_s <- hsd_phytr %>%
     drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_s$scientific.name)))
                
   # abundance-un-weighted phylo 
   # pd
   pd100_unweigh[i] <- sum(tr_s$edge.length)
   
   # mpd
   tr_s_dis <- cophenetic(tr_s)
   tr_s_dis <- tr_s_dis[lower.tri(tr_s_dis, diag = FALSE)]
   mpd100_unweigh[i] <- mean(tr_s_dis)
   
   # mntd
   tr_s_dis2 <- cophenetic(tr_s)
   diag(tr_s_dis2) <- NA
   mntd100_unweigh[i] <- mean(apply(tr_s_dis2, 2, min, na.rm = TRUE))
   
   # abundance-weighted phylo
   hsd_s_weigh <- hsd_s %>% 
     filter(scientific.name %in% tr_s$tip.label)
   
   tr_s_community <-
     matrix(table(hsd_s_weigh$scientific.name), 
            nrow = 1, 
            dimnames = list("s1", names(table(hsd_s_weigh$scientific.name))))
   
   # mpd
   mpd100_weigh[i] <- picante::mpd(tr_s_community, cophenetic(tr_s),
                                  abundance.weighted = TRUE)
   # mntd
   mntd100_weigh[i] <- picante::mntd(tr_s_community, cophenetic(tr_s),
                                    abundance.weighted = TRUE)
   
   
  # 50m circle
  hsd_sub <- hsd_s %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 50)
  
  # phylo tree with species in 50m circle
  tr_sub <- hsd_phytr %>%
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  # abundance-un-weighted phylo 
  # pd
  pd50_unweigh[i] <- sum(tr_sub$edge.length)
  
  # mpd
  tr_sub_dis <- cophenetic(tr_sub)
  tr_sub_dis <- tr_sub_dis[lower.tri(tr_sub_dis, diag = FALSE)]
  mpd50_unweigh[i] <- mean(tr_sub_dis)
  
  # mntd
  tr_sub_dis2 <- cophenetic(tr_sub)
  diag(tr_sub_dis2) <- NA
  mntd50_unweigh[i] <- mean(apply(tr_sub_dis2, 2, min, na.rm = TRUE))
  
  # abundance-weighted phylo
  hsd_sub_weigh <- hsd_sub %>% 
    filter(scientific.name %in% tr_sub$tip.label)
  
  tr_sub_community <-
    matrix(table(hsd_sub_weigh$scientific.name), 
           nrow = 1, 
           dimnames = list("s1", names(table(hsd_sub_weigh$scientific.name))))
  # mpd
  mpd50_weigh[i] <- picante::mpd(tr_sub_community, cophenetic(tr_sub),
                                 abundance.weighted = TRUE)
  # mntd
  mntd50_weigh[i] <- picante::mntd(tr_sub_community, cophenetic(tr_sub),
                                   abundance.weighted = TRUE)


    
  # 20m circle
  hsd_sub_sub <- hsd_sub %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 20)
  
  # phylo tree with species in 20m circle
  tr_sub_sub <- hsd_phytr %>%
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub_sub$scientific.name)))
  
  # abundance-un-weighted phylo 
  # pd
  pd20_unweigh[i] <- sum(tr_sub_sub$edge.length)
  
  # mpd
  tr_sub_sub_dis <- cophenetic(tr_sub_sub)
  tr_sub_sub_dis <- tr_sub_sub_dis[lower.tri(tr_sub_sub_dis, diag = FALSE)]
  mpd20_unweigh[i] <- mean(tr_sub_sub_dis)
  
  # mntd
  tr_sub_sub_dis2 <- cophenetic(tr_sub_sub)
  diag(tr_sub_sub_dis2) <- NA
  mntd20_unweigh[i] <- mean(apply(tr_sub_sub_dis2, 2, min, na.rm = TRUE))
  
  # abundance-weighted phylo
  hsd_sub_sub_weigh <- hsd_sub_sub %>% 
    filter(scientific.name %in% tr_sub_sub$tip.label)
  
  tr_sub_sub_community <-
    matrix(table(hsd_sub_sub_weigh$scientific.name), 
           nrow = 1, 
           dimnames = list("s1", names(table(hsd_sub_sub_weigh$scientific.name))))
  # mpd
  mpd20_weigh[i] <- picante::mpd(tr_sub_sub_community, cophenetic(tr_sub_sub),
                                 abundance.weighted = TRUE)
  # mntd
  mntd20_weigh[i] <- picante::mntd(tr_sub_sub_community, cophenetic(tr_sub_sub),
                                   abundance.weighted = TRUE)
  
  
  
  # 10m circle
  hsd_sub_sub_sub <- hsd_sub_sub %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 20)
  
  tr_sub_sub_sub <- hsd_phytr %>%
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub_sub_sub$scientific.name)))
  
  # abundance-un-weighted phylo 
  # pd
  pd10_unweigh[i] <- sum(tr_sub_sub_sub$edge.length)
  
  # mpd
  tr_sub_sub_sub_dis <- cophenetic(tr_sub_sub_sub)
  tr_sub_sub_sub_dis <- tr_sub_sub_sub_dis[lower.tri(tr_sub_sub_sub_dis, diag = FALSE)]
  mpd10_unweigh[i] <- mean(tr_sub_sub_sub_dis)
  
  # mntd
  tr_sub_sub_sub_dis2 <- cophenetic(tr_sub_sub_sub)
  diag(tr_sub_sub_sub_dis2) <- NA
  mntd10_unweigh[i] <- mean(apply(tr_sub_sub_sub_dis2, 2, min, na.rm = TRUE))
  
  # abundance-weighted phylo
  hsd_sub_sub_sub_weigh <- hsd_sub_sub_sub %>% 
    filter(scientific.name %in% tr_sub_sub_sub$tip.label)
  
  tr_sub_sub_sub_community <-
    matrix(table(hsd_sub_sub_sub_weigh$scientific.name), 
           nrow = 1, 
           dimnames = list("s1", names(table(hsd_sub_sub_sub_weigh$scientific.name))))
  # mpd
  mpd10_weigh[i] <- picante::mpd(tr_sub_sub_sub_community, cophenetic(tr_sub_sub_sub),
                                 abundance.weighted = TRUE)
  # mntd
  mntd10_weigh[i] <- picante::mntd(tr_sub_sub_sub_community, cophenetic(tr_sub_sub_sub),
                                   abundance.weighted = TRUE)
  
}

# data frame ----
sp_loc <- sp_loc %>% 
  mutate(
    pd100_unweigh = pd100_unweigh,
    mpd100_unweigh = mpd100_unweigh,
    mpd100_weigh = mpd100_weigh,
    mntd100_unweigh = mntd100_unweigh,
    mntd100_weigh = mntd100_weigh,
    pd50_unweigh = pd50_unweigh,
    mpd50_unweigh = mpd50_unweigh,
    mpd50_weigh = mpd50_weigh,
    mntd50_unweigh = mntd50_unweigh,
    mntd50_weigh = mntd50_weigh,
    pd20_unweigh = pd20_unweigh,
    mpd20_unweigh = mpd20_unweigh,
    mpd20_weigh = mpd20_weigh,
    mntd20_unweigh = mntd20_unweigh,
    mntd20_weigh = mntd20_weigh,
    pd10_unweigh = pd10_unweigh,
    mpd10_unweigh = mpd10_unweigh,
    mpd10_weigh = mpd10_weigh,
    mntd10_unweigh = mntd10_unweigh,
    mntd10_weigh = mntd10_weigh
  )

######
#删掉一些中间变量
#####

rm(hsd_s,hsd_s_weigh,tr_s, tr_s_dis,tr_s_dis2,tr_s_community,
  hsd_sub,hsd_sub_weigh,
   tr_sub,tr_sub_community,tr_sub_dis2,i,
   mntd50_unweigh,mntd50_weigh,mpd50_unweigh,
   mpd50_weigh,pd50_unweigh,tr_sub_dis,
   hsd_sub_sub,hsd_sub_sub_weigh,
   tr_sub_sub,tr_sub_sub_community,tr_sub_sub_dis2,
   mntd20_unweigh,mntd20_weigh,mpd20_unweigh,
   mpd20_weigh,pd20_unweigh,tr_sub_sub_dis,
   hsd_sub_sub_sub,hsd_sub_sub_sub_weigh,
   tr_sub_sub_sub,tr_sub_sub_sub_community,tr_sub_sub_sub_dis2,
   mntd10_unweigh,mntd10_weigh,mpd10_unweigh,
   mpd10_weigh,pd10_unweigh,tr_sub_sub_sub_dis)

