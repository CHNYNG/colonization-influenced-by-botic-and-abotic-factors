#####
#回归#
#####




#####
#加载之前数据
#####

#load("data/data_for_reg.RData")

#####
#未标准化的数据
#####
library(ggplot2)
library(dplyr)
library(tidyverse)


library(betareg)
beta <- betareg(am ~ 1+ sptype2, reg_sc, link = "logit")

#numeric_vars$qr_AM <- ifelse(numeric_vars$qr_AM <= 0, 0.01, ifelse(numeric_vars$qr_AM >= 1, 0.99, numeric_vars$qr_AM))
#numeric_vars$qr_EM <- ifelse(numeric_vars$qr_EM >= 1, 0.99, numeric_vars$qr_EM)

#library(betareg)
#AM_20
#reg$qr_AM <- ifelse(reg$qr_AM <= 0, 0.01, ifelse(reg$qr_AM >= 1, 0.99, reg$qr_AM))
#beita_model_forest <- betareg(qr_AM ~
 #                                mntd20_unweigh * mntd20_weigh * pd20_unweigh + ap * tp * soc * tn + SRL * SRA * AD + CBD_20, data = reg, link = "logit")
#summary(beita_model_forest)
library(glmmTMB)
glm_20 <- glmmTMB(am ~ pd20_unweigh + minpd_20 + apd_20 + avepd_20 + totpd_20 + ntpd_20 + shannon_div_20 + SRL + AD + mntd20_unweigh + sptype2, reg_sc, family=beta_family)
glm_20 <- glmmTMB(am ~  mntd20_unweigh + abundance, reg_sc, family=beta_family)
#glm_20 <- glmmTMB(am ~  apd_20+ ap + tp + soc + tn + SRL + SRA + AD + CBD_20 + (1|resource), reg_sc, family=beta_family)
glm_20b <- update(glm_20, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
all.equal(fixef(glm_20), fixef(glm_20b))
###例子
m1 <- glmmTMB(count~ mined, family=poisson, data=Salamanders)
m1B <- update(m1, control=glmmTMBControl(optimizer=optim,
                                         optArgs=list(method="BFGS")))
## estimates are *nearly* identical:
all.equal(fixef(m1), fixef(m1B))

#EM_10
beita_model_forest <- betareg(qr_EM ~
                         soc + mpd10_unweigh + tn + pd10_unweigh + DBH1 + shannon_div_10 + invsimpson_div_10 + GX, data = numeric_vars, link = "logit")

summary(beita_model_forest)
vif_values <- vif(beita_model_forest)
print(vif_values)


beita_model_cor <- betareg(qr_AM ~ 
                              AD + SRL + soc+ tp + ap + ph + pd50_unweigh +mpd50_unweigh + mpd50_weigh + mntd50_unweigh + mntd50_weigh + shannon_div_50 + BD_50 + CBD_50 , data = numeric_vars, link = "probit")
beita_model_cor <- betareg(qr_EM ~
                             DBH2 + AD + SRL + soc+ tp + ap + ph + pd10_unweigh +mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 + growth_rate, data = numeric_vars, link = "logit")

summary(beita_model_cor)

vif_values <- vif(beita_model_cor)
print(vif_values)



library(lme4)

reg_AM_10$qr_AM <- ifelse(reg_AM_10$qr_AM <= 0, 0.01, ifelse(reg_AM_10$qr_AM >= 1, 0.99, reg_AM_10$qr_AM))
reg_AM_10$Order <- as.factor(reg_AM_10$Order)
model <- betareg(qr_AM ~ AD + SRL + soc + tp + ap + ph + pd10_unweigh + mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 ,
                 data = reg_AM_10)
glmer(qr_AM ~ AD + SRL + soc+ tp + ap + ph + pd10_unweigh + mpd10_unweigh + mpd10_weigh + mntd10_unweigh + mntd10_weigh + shannon_div_10 + BD_10 + CBD_10 + Order , data = reg_AM_10, family = betabinomial(link = "logit"))

# 假设 "Order" 是一个因子型变量
reg_AM_10$Order <- factor(reg_AM_10$Order)

# 创建虚拟变量
dummy_vars <- model.matrix(~ Order - 1, data = reg_AM_10)

# 合并虚拟变量到数据框
reg_AM_10 <- cbind(reg_AM_10, dummy_vars)

# 构建Beta回归模型，包括虚拟变量
model <- betareg(qr_AM ~ AD * SRL * soc * tp * ap * ph * pd10_unweigh * mpd10_unweigh * mpd10_weigh * mntd10_unweigh * mntd10_weigh * shannon_div_10 * BD_10 * CBD_10 * Order, 
                 link = "cloglog", data = reg_AM_10)

# 查看模型摘要
summary(model)

####看一下自相关的距离
data <- reg 
coordinates(data) <- c("GX", "GY")
spatial_weights <- dnearneigh(data, d1 = 0, d2 = 50)
weights <- nb2listw(spatial_weights, style = "W", zero.policy = TRUE)
moran_result <- moran.test(data$soil_pc1, weights)
print(moran_result)

###整理每个

glm_20 <- glmmTMB(am ~  mntd20_unweigh + sptype2/resource, reg_sc, family=beta_family)
glm_20b <- update(glm_20, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))