##sem
#psem

# 筛选出特定列并删除带有NA的行
# 用am_beta_dat
# 加载包
library(dplyr)
library(piecewiseSEM)
library(betareg)
# 筛选出 am 列不为空值的行，并选择特定列

library(lme4)
library(dplyr)


set.seed(123)  # 确保结果可复现
am_beta_dat$random_effect <- sample(1:100, nrow(am_beta_dat), replace = TRUE)
am_beta_dat <- am_beta_dat %>%
  mutate(log_CBD_10 = log(CBD_10),
         log_DBH2 = log(DBH2))

# consider spatial random effect
# we construct 50m * 50m cells
# you can try other cells like 20 * 20
am_beta_dat$GX[am_beta_dat$GX == 1000] <- 999.99
am_beta_dat$GY[am_beta_dat$GY == 500] <- 499.99
am_beta_dat$GX[am_beta_dat$GX == 0] <- 0.01
am_beta_dat$GY[am_beta_dat$GY == 0] <- 0.01

am_beta_dat$cellx50 <- am_beta_dat$GX %/% 50 
am_beta_dat$celly50 <- am_beta_dat$GY %/% 50

am_beta_dat$c50 <- factor(paste((am_beta_dat$celly50 + 1001), 
                                sep = "_", 
                                (am_beta_dat$cellx50 + 1001)))

modelList <- psem(
  glmmTMB(am ~ minpd_10 + avepd_10 + DBH2 + CBD_10 + invsimpson_div_10 +
            RDi + pcoa1 + pcoa2 + rd_is + (1|c50),
    data = am_beta_dat,
    family = beta_family(link = "logit") ),
  lm(CBD_10 ~ minpd_10 + avepd_10 + invsimpson_div_10 + RDi +
       pcoa1 + pcoa2 + rd_is, data = am_beta_dat),
  lm(DBH2 ~ minpd_10 + avepd_10 + CBD_10 + invsimpson_div_10 +
       RDi + pcoa1 + pcoa2 + rd_is, data = am_beta_dat),
  CBD_10%~~%DBH2
)
summary(modelList)
plot(modelList)
am_col_psem <- psem(
  glmmTMB(
    am ~ minpd_10 + avepd_10 + DBH2 + CBD_10 + invsimpson_div_10 +
      RDi + pcoa1 + pcoa2 + (1 | c50),
    data = am_beta_dat,
    family = beta_family
  ),
  
  lm(CBD_10 ~ RDi + pcoa1 + pcoa2 + minpd_10 + avepd_10, data = am_beta_dat),
  
  lm(DBH2 ~ RDi + CBD_10 + pcoa1 + pcoa2 + minpd_10 + avepd_10, data = am_beta_dat),
  
  lm(invsimpson_div_10 ~ RDi + minpd_10 + avepd_10 , data = am_beta_dat),
  
  lm(minpd_10 ~ RDi + pcoa1 + pcoa2, data = am_beta_dat),
  
  lm(avepd_10 ~ RDi + pcoa1 + pcoa2, data = am_beta_dat),
  
  avepd_10 %~~% minpd_10,
  
  invsimpson_div_10 %~~% minpd_10,
  
  invsimpson_div_10 %~~% avepd_10,
  
  invsimpson_div_10 %~~% CBD_10
)
summary(am_col_psem)
plot(am_col_psem)
