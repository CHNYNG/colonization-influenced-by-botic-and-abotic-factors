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
pd20 <- c(); mpd20 <- c(); mntd20 <- c()
for (i in 1:dim(sp_loc)[1]) {
  
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((gx - sp_loc$gx[i])^2 + (gy - sp_loc$gy[i])^2) <= 20)
  
  # phylo tree with species in 20m circle
  tr_sub <- hsd_phytr %>% 
    drop.tip(setdiff(hsd_phytr$tip.label, unique(hsd_sub$scientific.name)))
  
  # pd
  pd20[i] <- sum(tr_sub$edge.length)
  
  # mpd
  tr_sub_dis <- cophenetic(tr_sub)
  tr_sub_dis <- tr_sub_dis[lower.tri(tr_sub_dis, diag = FALSE)]
  mpd20[i] <- mean(tr_sub_dis)
  
  # mntd
  tr_sub_dis2 <- cophenetic(tr_sub)
  diag(tr_sub_dis2) <- NA
  mntd20[i] <- mean(apply(tr_sub_dis2, 2, min, na.rm = TRUE))
}

# data frame ----
sp_loc <- sp_loc %>% 
  mutate(
    pd20 = pd20,
    mpd20= mpd20,
    mntd20 = mntd20
  )
