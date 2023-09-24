# data generation ----
focal_sp20 <- data.frame(
  gx = as.numeric(reg_AM_20$GX),
  gy = as.numeric(reg_AM_20$GY),
  latin_name = reg_AM_20$Latin
) %>%
  mutate(latin_name = gsub(" ", "_", latin_name))

hsd <- hsd %>%
  mutate(scientific.name = gsub(" ", "_", scientific.name))
hsd$scientific.name[which(hsd$scientific.name == "symplocos_congesta")] <-
  "Symplocos_congesta"

diff_name20 <- setdiff(focal_sp20$latin_name, hsd$scientific.name)
focal_sp20$latin_name[which(focal_sp20$latin_name == diff_name20[1])] <-
  "Cinnamomum_subavenium"
focal_sp20$latin_name[which(focal_sp20$latin_name == diff_name20[2])] <-
  "Wikstroemia_nutans"
focal_sp20$latin_name[which(focal_sp20$latin_name == diff_name20[3])] <-
  "Adina_pilulifera"
focal_sp20$latin_name[which(focal_sp20$latin_name == diff_name20[4])] <-
  "Lithocarpus_haipinii"


focal_sp20_neigh <- list()

for (i in 1:dim(focal_sp20)[1]) {
  focalSpecies <- focal_sp20$latin_name[i]
  focal_x <- focal_sp20$gx[i]
  focal_y <- focal_sp20$gy[i]
  
  # covert neighbor data to a one row matrix
  neigh_sub <- hsd %>%
    filter(sqrt((gx - focal_x) ^ 2 + (gy - focal_y) ^ 2) <= 20 &
             sqrt((gx - focal_x) ^ 2 + (gy - focal_y) ^ 2) > 0)
  
  neigh_sub_sp <- table(neigh_sub$scientific.name)
  
  #
  focal_sp20_neigh[[i]] <- data.frame(neigh_sub_sp) %>%
    column_to_rownames("Var1") %>%
    rename(Abundance = Freq) %>%
    t() %>%
    as.data.frame()
}

focal_sp20_neigh_spcom <- bind_rows(focal_sp20_neigh) %>%
  replace(is.na(.), 0)

rownames(focal_sp20_neigh_spcom) <-
  1:length(focal_sp20_neigh)

# beta calculation ----
focal_bray20 <- as.matrix(vegdist(focal_sp20_neigh_spcom, "bray"))
rownames(focal_bray20) <- str_c("focal", 1:2504)
colnames(focal_bray20) <- str_c("focal", 1:2504)

focal_bray20_dat <-
  melt(focal_bray20,
       varnames = c("focali", "focalj"),
       value.name = "beta") %>%
  subset(focali != focalj)

# distance of AM_20
am20_rate_dis <- as.matrix(dist(reg_AM_20$qr_AM,
                                diag = T,
                                method = "euclidean"))
rownames(am20_rate_dis) <- str_c("focal", 1:2504)
colnames(am20_rate_dis) <- str_c("focal", 1:2504)

am20_rate_dis_dat <-
  melt(am20_rate_dis,
       varnames = c("focali", "focalj"),
       value.name = "am_rate_diff") %>%
  subset(focali != focalj) %>%
  left_join(focal_bray20_dat)

cor.test(~ am20_rate_dis_dat$beta + am20_rate_dis_dat$am_rate_diff)

# mantel test ----
# easier way
mantel1 <- mantel(
  vegdist(focal_sp20_neigh_spcom, method = "bray"),
  dist(reg_AM_20$qr_AM, method = "euclidean"),
  permutations = 999
)



#####
#谱系距离
#####
#成对物种间的边长
library(ape)
edge_lengths <- cophenetic(hsd_phytr)
edge_lengths <- as.data.frame(edge_lengths)
edge_lengths$RowName <- rownames(edge_lengths)
edge_lengths <- melt(edge_lengths, id.vars = "RowName")
colnames(edge_lengths) <- c("RowName", "ColName", "Value")
# 删除Value为0的行
edge_lengths <- edge_lengths %>%
  filter(Value != 0)
edge_lengths <- edge_lengths %>%
  mutate(
    RowName = as.character(RowName),
    ColName = as.character(ColName),
    RowName = ifelse(RowName < ColName, RowName, ColName),
    ColName = ifelse(RowName < ColName, ColName, RowName)
  ) %>% 
  distinct(RowName, ColName, .keep_all = TRUE)
