# new beta regression by weitao 
# required packages ----
library(tidyverse)
library(glmmTMB)
library(piecewiseSEM)
library(performance)
library(MuMIn)

load("D:/Users/65205/Documents/WeChat Files/wxid_ulzapgmya6hv22/FileStorage/File/2024-06/cy 6.7.RData")
# data preparation
# origninal data set is reg_scaled
# we can exclude ecm species and re-scaled data
# try how to get the pcoa
#first obtain the distance matrix of the species phylo
library(ape)
library(phangorn)
species_to_keep <- am_phylo_pcoa_axes$Latin
filtered_tree <- keep.tip(hsd_phytr, species_to_keep)
dist_matrix <- cophenetic(filtered_tree)
#dist_matrix <- scale(dist_matrix)

#second, play the pcoa
library(vegan)
#pcoa_result <- pcoa(dist_matrix)
pcoa_result <- cmdscale(dist_matrix) 


am_phylo_pcoa_axes_try <- as.data.frame(as.data.frame(pcoa_result))
am_phylo_pcoa_axes_try <- decostand(am_phylo_pcoa_axes_try, "normalize", MARGIN = 2)*50
am_phylo_pcoa_axes_try <- as.data.frame(am_phylo_pcoa_axes_try)
summary(am_phylo_pcoa_axes_try)
summary(am_phylo_pcoa_axes)
#pcoa_result[["Value"]]
#把pcoa的结果中心化一下？

am_beta_dat <- reg_scaled %>%
  filter(!is.na(am)) %>%
  select(
    # species identity
    Latin,
    Genus,
    Family,
    # resource diversity
    RDi, RRi, REi,
    # size
    DBH2,
    # root traits
    AD,
    SRL,
    SRA,
    # topography
    elevation,
    slope,
    aspect,
    convexity,
    # phylogenetic
    minpd_10,
    avepd_10,
    totpd_10,
    # neighbor effect
    CBD_10,
    BD_10,
    invsimpson_div_10,
  ) %>%
  left_join(am_phylo_pcoa_axes) %>%
  mutate_if(is.numeric, scale) %>%
  mutate(
    rd_is = scale(RDi * invsimpson_div_10),
    rr_is = scale(RRi * invsimpson_div_10),
    re_is = scale(REi * invsimpson_div_10),
    am_beta_dat %>% filter(!is.na(am)) %>% select(am)
  )

# I only calculate the pcoa axes of 161 species, and put them in: am_phylo_pcoa_axes

# beta regression and model selection
# based on your manuscript, we can set up the full model like this:
library(betareg)
options(na.action = "na.omit")
am_beta_mod_ful <- betareg(
  am ~ minpd_10 + avepd_10 + totpd_10 + AD + SRL + SRA + DBH2 + CBD_10 + BD_10 + 
    elevation + invsimpson_div_10 + RDi + pcoa1 + pcoa2 + rd_is,
  data = am_beta_dat
)
summary(am_beta_mod_ful)

# R can automatically choose optimal model for us based on AIC
options(na.action = "na.fail")
am_beta_mod_sub <- dredge(am_beta_mod_ful)


subset(am_beta_mod_sub, delta < 2) # optimal model and those with delta AICc < 2
am_beta_mod_minaic <- summary(get.models(am_beta_mod_sub, 1)[[1]]) # optimal model
am_beta_mod_ave <- summary(model.avg(am_beta_mod_sub, subset = delta < 2)) # averaging method
# am_beta_mod_ave is another choice if you like, see more in Burnham, K. P. and Anderson, D. R. 2002

# if you want to make everything simple, use am_beta_mod_minaic
# and based on am_beta_mod_minaic:
am_beta_mod_optimal <- betareg(
  am ~ minpd_10 + avepd_10 + DBH2 + CBD_10 + invsimpson_div_10 + RDi + pcoa1 + pcoa2 + rd_is,
  data = am_beta_dat
)

# you will find it's wired that in am_beta_mod_minaic, the RDi is excluded
# but in am_beta_mod_optimal, I included it!
# The reason is that if we keep the interaction term in the model,
# it's better to keep both of the original variables that create the interaction term.


# so, basically am_beta_mod_optimal is the final model you should report in your article
# it contains some variables that are not statistically significant,
# it doesn't matter
# constructing your psem with all the variables in am_beta_mod_optimal


#psem
library(piecewiseSEM)
# 筛选出特定列并删除带有NA的行

# 加载 dplyr 包
library(dplyr)

# 筛选出 am 列不为空值的行，并选择特定列
reg_sc_psem <- reg_scaled %>%
  filter(!is.na(am)) %>%
  select(am, minpd_10, avepd_10, totpd_10, SRA, DBH2, CBD_10, invsimpson_div_10, RDi, soc, tn, tp, ap, ph, moisture,
         Order, Family, Genus, Latin)

library(lme4)

#####直接再重新做一个psem
##玩我…
#把pd做一个复合函数
PD1 <- lm(am ~ minpd_10 + avepd_10 + totpd_10 , reg_sc_psem)
summary(PD1)
coefs(Neighbor, standardize = 'scale')
beta_minpd_10 <- summary(PD1)$coefficients[2,1]
beta_avepd_10 <- summary(PD1)$coefficients[3,1]
beta_totpd_10 <- summary(PD1)$coefficients[4,1]
PD <- beta_minpd_10*reg_sc_psem$minpd_10 + beta_avepd_10*reg_sc_psem$avepd_10 +
  beta_totpd_10*reg_sc_psem$totpd_10
reg_sc_psem$PD <- PD
summary(lm(am ~ PD, reg_sc_psem))
coefs(lm(am ~ PD, reg_sc_psem))

#其他的都分别来
modelList <- psem(
  glm(am ~ PD + CBD_10 + invsimpson_div_10 + DBH2 + SRA +  RDi , data = reg_sc_psem,family = gaussian(link = "identity")),
  lm(DBH2  ~ PD + CBD_10 + RDi + SRA  , data = reg_sc_psem),
  lm(CBD_10 ~ RDi + PD + invsimpson_div_10 + SRA, data = reg_sc_psem),
  lm(RDi ~ PD + invsimpson_div_10 + SRA, data = reg_sc_psem),
  lm(PD ~ invsimpson_div_10 + SRA, data = reg_sc_psem),
  SRA %~~% invsimpson_div_10
)
summary(modelList)
plot(modelList)



#啊啊啊啊，检验之后不行哇！！！！完全不行哇！！！
#glm_10 <- glmmTMB(am ~ minpd_10 + avepd_10 + totpd_10 + SRA + DBH2 + CBD_10 + invsimpson_div_10 * RDi +(1|Family),
#                 reg_scaled, family = ordbeta)
#summary(glm_10)
