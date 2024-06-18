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

am_beta_dat <- am_beta_dat %>%
  mutate(across(c(rd_is, rr_is, re_is), as.numeric))

set.seed(123)  # 确保结果可复现
am_beta_dat$random_effect <- sample(1:100, nrow(am_beta_dat), replace = TRUE)
am_beta_dat <- am_beta_dat %>%
  mutate(log_CBD_10 = log(CBD_10),
         log_DBH2 = log(DBH2))

modelList <- psem(
  glmmTMB(am ~ minpd_10 + avepd_10 + DBH2 + CBD_10 + invsimpson_div_10 +
            RDi + pcoa1 + pcoa2 + rd_is ,
    data = am_beta_dat,
    family = beta_family(link = "logit") ),
  lm(CBD_10 ~ minpd_10 + avepd_10 + invsimpson_div_10 + RDi +
       pcoa1 + pcoa2 + rd_is, data = am_beta_dat),
  lm(DBH2 ~ minpd_10 + avepd_10 + CBD_10 + invsimpson_div_10 +
       RDi + pcoa1 + pcoa2 + rd_is, data = am_beta_dat),
  CBD_10%~~%DBH2
)
summary(modelList)
