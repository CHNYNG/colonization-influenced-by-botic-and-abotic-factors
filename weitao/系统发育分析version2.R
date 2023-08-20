# 
library(picante)
library(tidyverse)

# data preparation ----
# location of 3552 samples
sp_loc <- data.frame(
  gx = as.numeric(prd.loc$GX),
  gy = prd.loc$GY
)

# phylo tree data
hsd_phytr <- read.tree("Data/trephylo 29rhmc tre by zhouwen.nwk")

# hsd 2nd census data
# combining different branches
# hsd_alive_singlebr is the final data frame we should use
alive_status <- c("Alive", "P", "lean", "alive", "snag")
hsd_alive <- hsd %>%
  filter(status2 %in% alive_status, 
         !is.na(tagnew), 
         !is.na(branch), 
         !is.na(dbh2))


hsd_multibr <- hsd_alive %>%
  group_by(tagnew) %>%
  summarise(dbh_multi = sqrt(sum(dbh2 ^ 2)))

hsd_alive_singlebr <- hsd_alive %>% 
  filter(branch == 0) %>% 
  left_join(hsd_multibr, by = "tagnew") %>% 
  mutate(
    scientific.name = str_replace(scientific.name, " ", "_")
  )

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

# data frame ----
sp_loc <- sp_loc %>% 
  mutate(
    pd20_unweigh = pd20_unweigh,
    mpd20_unweigh = mpd20_unweigh,
    mpd20_weigh = mpd20_weigh,
    mntd20_unweigh = mntd20_unweigh,
    mntd20_weigh = mntd20_weigh
  )

