# loops 2 ----
# this time we will compute:

# totpd: the sum of the phylogenetic distances between a focal tree and its heterospecific neighbors.
# totpd:focal与其异种邻近树之间的系统发育距离的总和
# avepd: the average of the phylogenetic distances between a focal tree and its heterospecific neighbors.
# avepd: focal与其异种邻近树之间系统发育距离的平均值。
# min同理

# apd: the deviation of observed average phylogenetic
# 观察到的平均系统发育的偏差
## distance between a focal tree and its heterospecific
## 也就是focal与其异种树之间的距离
## neighbors from that expected under a null model
## 在零模型下的期望值

# ntpd: the equivalent deviation of observed phylogenetic distance
# 观察到的系统发育距离的等效偏差
# between a focal seedling and its most closely
# 在一个focal和它最密切的之间异种邻居
# related heterospecific neighbor

library(picante)
library(dplyr)

focal_sp<- sp_loc %>%
  left_join(., hsd_alive_singlebr, "TagNew") %>%
  select(GX.x, GY.x, scientific.name) %>%
  rename(
    GX = GX.x,
    GY = GY.x
  )


  # null model of phy tree
  null.phy <- function(phytree, n, focal_label, neigh_label) {
    avepd_null <- c()
    minpd_null <- c()
    
    for(j in 1:n) {
      null_tree <- ladderize(phytree)
      
      null_label <- sample(phytree$tip.label)
      
      null_tree$tip.label <- null_label
      
      phy_dist_null <- cophenetic(null_tree)
      
      focal_sp_in_null_phydist <- 
        which(null_tree$tip.label == focal_label)
      
      neigh_sp_in_null_phydist <- 
        which(null_tree$tip.label %in% neigh_label)
      
      
      avepd_null[j] <- 
        mean(phy_dist_null[focal_sp_in_null_phydist, neigh_sp_in_null_phydist])
      
      minpd_null[j] <- 
        min(phy_dist_null[focal_sp_in_null_phydist, neigh_sp_in_null_phydist])
    }
    
    phy_null_dat <- data.frame(
      avepd_null,
      minpd_null
    )
    
    return(phy_null_dat)
  }


# compute indicies
totpd_20 <- c()
avepd_20 <- c()
minpd_20 <- c()
apd_20 <- c() 
ntpd_20 <- c()

for (i in 1:dim(focal_sp)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - focal_sp$GX[i])^2 + (GY - focal_sp$GY[i])^2) <= 20)
  
  # phylo tree with species in 20m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  phy_dist <- cophenetic(tr_sub)
  
  # focal species and neighbor species in phy tree
  focal_sp_label <- focal_sp$scientific.name[i]
  
  neigh_sp_label <- tr_sub$tip.label[tr_sub$tip.label != focal_sp_label]
  
  focal_sp_in_phydist <- which(tr_sub$tip.label == focal_sp_label)
  
  neigh_sp_in_phydist <- which(tr_sub$tip.label %in% neigh_sp_label)
  
  
  # compute phy indicies 
  totpd_20[i] <- sum(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  avepd_20[i] <- mean(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  minpd_20[i] <- min(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  #
  null_dat_sub <- null.phy(tr_sub, 999, focal_sp_label, neigh_sp_label)
  
  apd_20[i] <- 
    (avepd_20[i] - mean(null_dat_sub$avepd_null)) / sd(null_dat_sub$avepd_null)
  
  ntpd_20[i] <- 
    (minpd_20[i] - mean(null_dat_sub$minpd_null)) / sd(null_dat_sub$minpd_null)
  
}

totpd_10 <- c()
avepd_10 <- c()
minpd_10 <- c()
apd_10 <- c() 
ntpd_10 <- c()

for (i in 1:dim(focal_sp)[1]) {
  
  # 10m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - focal_sp$GX[i])^2 + (GY - focal_sp$GY[i])^2) <= 10)
  
  # phylo tree with species in 10m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  phy_dist <- cophenetic(tr_sub)
  
  # focal species and neighbor species in phy tree
  focal_sp_label <- focal_sp$scientific.name[i]
  
  neigh_sp_label <- tr_sub$tip.label[tr_sub$tip.label != focal_sp_label]
  
  focal_sp_in_phydist <- which(tr_sub$tip.label == focal_sp_label)
  
  neigh_sp_in_phydist <- which(tr_sub$tip.label %in% neigh_sp_label)
  
  
  # compute phy indicies 
  totpd_10[i] <- sum(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  avepd_10[i] <- mean(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  minpd_10[i] <- min(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  #
  null_dat_sub <- null.phy(tr_sub, 999, focal_sp_label, neigh_sp_label)
  
  apd_10[i] <- 
    (avepd_10[i] - mean(null_dat_sub$avepd_null)) / sd(null_dat_sub$avepd_null)
  
  ntpd_10[i] <- 
    (minpd_10[i] - mean(null_dat_sub$minpd_null)) / sd(null_dat_sub$minpd_null)
  
}

totpd_50 <- c()
avepd_50 <- c()
minpd_50 <- c()
apd_50 <- c() 
ntpd_50 <- c()

for (i in 1:dim(focal_sp)[1]) {
  
  # 50m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - focal_sp$GX[i])^2 + (GY - focal_sp$GY[i])^2) <= 50)
  
  # phylo tree with species in 50m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  phy_dist <- cophenetic(tr_sub)
  
  # focal species and neighbor species in phy tree
  focal_sp_label <- focal_sp$scientific.name[i]
  
  neigh_sp_label <- tr_sub$tip.label[tr_sub$tip.label != focal_sp_label]
  
  focal_sp_in_phydist <- which(tr_sub$tip.label == focal_sp_label)
  
  neigh_sp_in_phydist <- which(tr_sub$tip.label %in% neigh_sp_label)
  
  
  # compute phy indicies 
  totpd_50[i] <- sum(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  avepd_50[i] <- mean(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  minpd_50[i] <- min(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  #
  null_dat_sub <- null.phy(tr_sub, 999, focal_sp_label, neigh_sp_label)
  
  apd_50[i] <- 
    (avepd_50[i] - mean(null_dat_sub$avepd_null)) / sd(null_dat_sub$avepd_null)
  
  ntpd_50[i] <- 
    (minpd_50[i] - mean(null_dat_sub$minpd_null)) / sd(null_dat_sub$minpd_null)
  
}

totpd_100 <- c()
avepd_100 <- c()
minpd_100 <- c()
apd_100 <- c() 
ntpd_100 <- c()

for (i in 1:dim(focal_sp)[1]) {
  
  # 100m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - focal_sp$GX[i])^2 + (GY - focal_sp$GY[i])^2) <= 100)
  
  # phylo tree with species in 100m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  phy_dist <- cophenetic(tr_sub)
  
  # focal species and neighbor species in phy tree
  focal_sp_label <- focal_sp$scientific.name[i]
  
  neigh_sp_label <- tr_sub$tip.label[tr_sub$tip.label != focal_sp_label]
  
  focal_sp_in_phydist <- which(tr_sub$tip.label == focal_sp_label)
  
  neigh_sp_in_phydist <- which(tr_sub$tip.label %in% neigh_sp_label)
  
  
  # compute phy indicies 
  totpd_100[i] <- sum(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  avepd_100[i] <- mean(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  minpd_100[i] <- min(phy_dist[focal_sp_in_phydist, neigh_sp_in_phydist])
  
  #
  null_dat_sub <- null.phy(tr_sub, 999, focal_sp_label, neigh_sp_label)
  
  apd_100[i] <- 
    (avepd_100[i] - mean(null_dat_sub$avepd_null)) / sd(null_dat_sub$avepd_null)
  
  ntpd_100[i] <- 
    (minpd_100[i] - mean(null_dat_sub$minpd_null)) / sd(null_dat_sub$minpd_null)
  
}



pd_ind_all <- data.frame(
#  totpd_20, avepd_20, minpd_20, apd_20, ntpd_20,
  totpd_10, avepd_10, minpd_10, apd_10, ntpd_10,
#  totpd_50, avepd_50, minpd_50, apd_50, ntpd_50,
#  totpd_100, avepd_100, minpd_100, apd_100, ntpd_100
)

rm(focal_sp, focal_sp_in_phydist, null_dat_sub, phy_dist, tr_sub, apd_10, apd_20, apd_50,apd_100,
   focal_sp_label,i, minpd_10, minpd_20, minpd_50, minpd_100, neigh_sp_in_phydist, neigh_sp_label, ntpd_10, ntpd_20, ntpd_50, ntpd_100,
   totpd_10, totpd_20, totpd_50, totpd_100, avepd_10, avepd_20, avepd_50, avepd_100)
