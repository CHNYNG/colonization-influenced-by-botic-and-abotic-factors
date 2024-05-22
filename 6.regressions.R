#####
#回归#
#####

#####
#am
#####

#####
#10#####
library(glmmTMB)
#apd_10 + RDi +  ntpd_10  +  avepd_10 + totpd_10
#RDi + totpd_10 + moisture + avepd_10 + mntd10_unweigh + elevation + shannon_div_10 + SRL + mpd10_weigh + mntd10_weigh
#avepd_10 + totpd_10 + RRi + moisture + elevation + shannon_div_10 + mntd10_unweigh + REi + pd10_unweigh + mntd10_weigh
#从逻辑顺序上来添加变量
###使用ntpd + apd
#totpd_10 + avepd_10 + minpd_10
#pd10_unweigh + mpd10_unweigh + mntd10_unweigh
#mpd10_weigh + mntd10_weigh
#minpd_10 + SRA + DBH2 + CBD_10 + shannon_div_10 * RDi
glm_10 <- glmmTMB(am ~ minpd_10 + avepd_10 + totpd_10 + SRA + DBH2 + CBD_10 + invsimpson_div_10* RDi,
                  reg_scaled, family = beta_family)
summary(glm_10)
library(lme4)

glm_10_glm <- glm(am ~ minpd_10 + avepd_10 + totpd_10 + SRA + DBH2 + CBD_10 + + invsimpson_div_10 + RDi + invsimpson_div_10 : RDi,
          data = reg_scaled, family = gaussian)
summary(glm_10_glm)
library(glmm.hp)
hp_results <- glmm.hp(glm_10_glm) 
partition_data <- hp_results$hierarchical.partitioning
variable_names <- hp_results$variables
library(ggplot2)

ggplot(partition_data, aes(x = reorder(variable_names, -Individual), y = Individual)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "darkblue") +  # 更改填充色和边框色
#  geom_text(aes(label = round(Individual, 2)), vjust = -0.3, color = "black", size = 3) +  # 添加柱形图上的数据标签
  labs(title = "Individual Predictor Contribution to R-squared", 
       x = "Predictor", y = "Proportion of Explained Variance") +
  scale_x_discrete(labels = c("DBH", "minpd", "CBD", "RDi", "avepd", "SRA", "invsimpon:RDi", "invsimpson", "totpd")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1)) +  # 调整横坐标标签的大小、颜色和角度
  theme(axis.text.y = element_text(size = 10, color = "black")) +  # 调整纵坐标标签的大小和颜色
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +  # 调整标题的大小和位置
  theme(plot.subtitle = element_text(size = 12, hjust = 0.5)) +  # 调整副标题的大小和位置
  theme(legend.position = "none")  # 移除图例

plot(glmm.hp(glm_10_glm))

library(piecewiseSEM)
# 筛选出特定列并删除带有NA的行
reg_sc_psem <- na.omit(reg_scaled[, c("am", "minpd_10", "avepd_10", "totpd_10", "SRA", "DBH2", "CBD_10", "invsimpson_div_10", "RDi", "soc", "tn", "tp","ap","ph","moisture")])

####构建第一个复合变量Neighbor
Neighbor <- lm(am ~ minpd_10 + avepd_10 + totpd_10 + CBD_10 + invsimpson_div_10, reg_sc_psem)
summary(Neighbor)
coefs(Neighbor, standardize = 'scale')
beta_minpd_10 <- summary(Neighbor)$coefficients[2,1]
beta_avepd_10 <- summary(Neighbor)$coefficients[3,1]
beta_totpd_10 <- summary(Neighbor)$coefficients[4,1]
beta_CBD_10 <- summary(Neighbor)$coefficients[5,1]
beta_invsimpson_div_10 <- summary(Neighbor)$coefficients[6,1]
Neighbor1 <- beta_minpd_10*reg_sc_psem$minpd_10 + beta_avepd_10*reg_sc_psem$avepd_10 +
  beta_totpd_10*reg_sc_psem$totpd_10 +
  beta_CBD_10*reg_sc_psem$CBD_10 + beta_invsimpson_div_10*reg_sc_psem$invsimpson_div_10
reg_sc_psem$Neighbor1 <- Neighbor1
summary(lm(am ~ Neighbor1, reg_sc_psem))
coefs(lm(am ~ Neighbor1, reg_sc_psem))
#####构建第二个复合变量Host
Host <- lm(am ~ DBH2 + SRA, reg_sc_psem)
summary(Host)
coefs(Host, standardize = 'scale')
beta_DBH2 <- summary(Host)$coefficients[2,1]
beta_SRA <- summary(Host)$coefficients[3,1]
Host1 <- beta_DBH2*reg_sc_psem$DBH2 + beta_SRA*reg_sc_psem$SRA
reg_sc_psem$Host1 <- Host1
summary(lm(am ~ Host1, reg_sc_psem))
coefs(lm(am ~ Host1, reg_sc_psem))
#####构建第三个复合变量Environment
#先尝试一下只用RDi，之后再尝试使用土壤因子
#用RDi好像不行，那就把各种元素都加进来试试
#我有这些个元素"soc", "tn", "tp","ap","ph","moisture" 
Soil <- lm(am ~ soc + tn + tp +ap + ph + moisture, reg_sc_psem)
summary(Soil)
coefs(Soil, standardize = 'scale')
beta_soc <- summary(Soil)$coefficients[2,1]
beta_tn <- summary(Soil)$coefficients[3,1]
beta_tp <- summary(Soil)$coefficients[4,1]
beta_ap <- summary(Soil)$coefficients[5,1]
beta_ph <- summary(Soil)$coefficients[6,1]
beta_moisture <- summary(Soil)$coefficients[7,1]
Soil1 <- beta_soc*reg_sc_psem$soc + beta_tn*reg_sc_psem$tn +
  beta_tp*reg_sc_psem$tp + beta_ap*reg_sc_psem$ap + beta_ph*reg_sc_psem$ph +
  beta_moisture*reg_sc_psem$moisture
reg_sc_psem$Soil1 <- Soil1
summary(lm(am ~ Soil1, reg_sc_psem))
coefs(lm(am ~ Soil1, reg_sc_psem))
#####拟合复合变量间的多元回归模型
modelList <- psem(
  glm(am ~ Neighbor1 + Host1 +  Soil1 , data = reg_sc_psem),
  lm(Host1 ~  Soil1 ,data = reg_sc_psem),
  lm(Neighbor1 ~ Soil1 , data = reg_sc_psem)
)

summary(modelList, .progressBar = FALSE)
pdf("pic/sem.pdf")
plot(modelList)
dev.off()

#啊啊啊啊，检验之后不行哇！！！！完全不行哇！！！
#glm_10 <- glmmTMB(am ~ minpd_10 + avepd_10 + totpd_10 + SRA + DBH2 + CBD_10 + invsimpson_div_10 * RDi +(1|Family),
 #                 reg_scaled, family = ordbeta)
#summary(glm_10)
library(DHARMa)
glm_10_simres <- simulateResiduals(glm_10)
plot(glm_10_simres)
library(car)
Anova(glm_10) ## default type II
Anova(glm_10,type="III")
library(effects)
ae <- allEffects(glm_10)
plot(ae)
glm_10_sr <- glmmTMB(am ~ shannon_div_10 * RDi,
                     reg_sc, family = ordbeta)

plot(allEffects(glm_10_sr))
#画图

library(ggplot2)
coef_data <- summary(glm_10_2)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "invsimpson_div_10:RDi", "RDi", "invsimpson_div_10", "CBD_10", "DBH2", "SRA", "avepd_10", "totpd_10", "minpd_10", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), color = "#bcbddc", height = 0.2, linewidth = 2) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "#756bb1", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20)) + #调整所有文字的大小
  scale_y_discrete(labels = c( "IS:RDi", "RDi", "IS", "CBD", "DBH", "SRA", "avepd", "totpd", "minpd", "(Intercept)"))
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围


glm <- glmmTMB(am ~ avepd_10 + totpd_10 + RRi + moisture + elevation + shannon_div_10 + mntd10_unweigh + REi + pd10_unweigh + mntd10_weigh + (1| sptype2), reg_sc, family = beta_family)
glm_10_M_1 <- glmmTMB(am ~  apd_10 + RDi +  ntpd_10  +  avepd_10 + totpd_10  + (1|sptype2*Genus), reg_sc, family=beta_family)
glm_10_M_1b <- update(glm_10_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_1)
glm_10_M_2 <- glmmTMB(am ~  mntd10_unweigh + AD + aspect + DBH2 + totpd_10 + pd10_unweigh + avepd_10 + CBD_10 + minpd_10 + (1|resource), reg_sc, family=beta_family)
glm_10_M_2b <- update(glm_10_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_2b)
glm_10_M_3 <- glmmTMB(am ~  mntd10_unweigh + AD + aspect + DBH2 + totpd_10 + pd10_unweigh + avepd_10 + CBD_10 + minpd_10 + (1|sptype2), reg_sc, family=beta_family)
glm_10_M_3b <- update(glm_10_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_3b)
glm_10_M_4 <- glmmTMB(am ~  mntd10_unweigh + AD + aspect + DBH2 + totpd_10 + pd10_unweigh + avepd_10 + CBD_10 + minpd_10, reg_sc, family=beta_family)
glm_10_M_4b <- update(glm_10_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_M_4b)

#mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10
glm_10_P_1 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10_P_1b <- update(glm_10_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_1b)
glm_10_P_2 <- glmmTMB(am ~ mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|resource), reg_sc, family=beta_family)
glm_10_P_2b <- update(glm_10_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_2b)
glm_10_P_3 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|sptype2), reg_sc, family=beta_family)
glm_10_P_3b <- update(glm_10_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_3b)
glm_10_P_4 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|sptype2), reg_sc, family=beta_family)
glm_10_P_4b <- update(glm_10_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_4b)
glm_10_P_5 <- glmmTMB(am ~  mntd10_unweigh + aspect + DBH2 + AD + pd10_unweigh + totpd_10 + avepd_10 + minpd_10 + (1|Genus), reg_sc, family=beta_family)
glm_10_P_5b <- update(glm_10_P_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10_P_5)
#####
#20#####

#apd_20 + minpd_20 +  ap + soc + tp + pd20_unweigh + moisture + mntd20_unweigh
glm_20 <- glmmTMB(am ~ minpd_20 + avepd_20 + totpd_20 + SRA + DBH2 + CBD_20 + shannon_div_20 * RDi, reg_sc, family = beta_family)
summary(glm_20)

#画图
library(ggplot2)
coef_data <- summary(glm_20)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_20:RDi", "RDi", "shannon_div_20", "CBD_20", "DBH2", "SRA", "avepd_20", "totpd_20", "minpd_20", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "blue", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围

glm_20_M_1 <- glmmTMB(am ~ apd_20 + minpd_20 +  ap + soc + tp + pd20_unweigh + moisture + mntd20_unweigh + (1|sptype2), reg_sc, family=beta_family)
glm_20_M_1b <- update(glm_20_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_1)
glm_20_M_2 <- glmmTMB(am ~  apd_20 + minpd_20 + CBD_20 + RDi + totpd_20  + (1|resource), reg_sc, family=beta_family)
glm_20_M_2b <- update(glm_20_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_2b)
glm_20_M_3 <- glmmTMB(am ~  apd_20 + minpd_20 + CBD_20 + RDi + totpd_20  + (1|sptype2), reg_sc, family=beta_family)
glm_20_M_3b <- update(glm_20_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_3)
glm_20_M_4 <- glmmTMB(am ~  apd_20 + minpd_20 + CBD_20 + RDi + totpd_20 + (1|Genus) , reg_sc, family=beta_family)
glm_20_M_4b <- update(glm_20_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_M_4)
#mntd20_unweigh + AD + CBD_20 + aspect + avepd_20 + pd20_unweigh + DBH2 + minpd_20
glm_20_P_1 <- glmmTMB(am ~  mntd20_unweigh + AD + aspect + CBD_20 + avepd_20 + DBH2 + minpd_20 + pd20_unweigh + apd_20 + ap + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20_P_1b <- update(glm_20_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_1b)
glm_20_P_2 <- glmmTMB(am ~  mntd20_unweigh + AD + aspect + CBD_20 + avepd_20 + DBH2 + minpd_20 + pd20_unweigh + apd_20 + ap + (1|resource), reg_sc, family=beta_family)
glm_20_P_2b <- update(glm_20_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_2b)
glm_20_P_3 <- glmmTMB(am ~  mntd20_unweigh + AD + aspect + CBD_20 + avepd_20 + DBH2 + minpd_20 + pd20_unweigh + apd_20 + ap + (1|sptype2), reg_sc, family=beta_family)
glm_20_P_3b <- update(glm_20_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_3b)
glm_20_P_4 <- glmmTMB(am ~  mntd20_unweigh + AD + aspect + CBD_20 + avepd_20 + DBH2 + minpd_20 + pd20_unweigh + apd_20 + ap , reg_sc, family=beta_family)
glm_20_P_4b <- update(glm_20_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20_P_4b)

#####
#50#####
#apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect
glm_50 <- glmmTMB(am ~ minpd_50 + avepd_50 + totpd_50 + SRA + DBH2 + CBD_50 + shannon_div_50 * RDi, reg_sc, family = beta_family)
summary(glm_50)

#画图
library(ggplot2)
coef_data <- summary(glm_50)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_50:RDi", "RDi", "shannon_div_50", "CBD_50", "DBH2", "SRA", "avepd_50", "totpd_50", "minpd_50", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "blue", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围


glm_50_M_1 <- glmmTMB(am ~  apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect  + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50_M_1b <- update(glm_50_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_1b)
glm_50_M_2 <- glmmTMB(am ~  apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect + (1|resource), reg_sc, family=beta_family)
glm_50_M_2b <- update(glm_50_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_2b)
glm_50_M_3 <- glmmTMB(am ~  apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect + (1|sptype2), reg_sc, family=beta_family)
glm_50_M_3b <- update(glm_50_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_3b)
glm_50_M_4 <- glmmTMB(am ~  apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect , reg_sc, family=beta_family)
glm_50_M_4b <- update(glm_50_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_4b)
glm_50_M_5 <- glmmTMB(am ~  apd_50 + CBD_50 + ap + AD + minpd_50 + tp + aspect + (1|Type * Genus), reg_sc, family=beta_family)
glm_50_M_5b <- update(glm_50_M_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_M_5)
#apd_50 + AD + CBD_50 + aspect + minpd_50 + ap + mntd50_unweigh + DBH2  
glm_50_P_1 <- glmmTMB(am ~  apd_50 + AD + CBD_50 + aspect + minpd_50 + ap + mntd50_unweigh + DBH2     + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50_P_1b <- update(glm_50_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_1b)
glm_50_P_2 <- glmmTMB(am ~  apd_50 + AD + CBD_50 + aspect + minpd_50 + ap + mntd50_unweigh + DBH2   + (1|resource), reg_sc, family=beta_family)
glm_50_P_2b <- update(glm_50_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_2b)
glm_50_P_3 <- glmmTMB(am ~  apd_50 + AD + CBD_50 + aspect + minpd_50 + ap + mntd50_unweigh + DBH2  + (1|sptype2), reg_sc, family=beta_family)
glm_50_P_3b <- update(glm_50_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_3b)
glm_50_P_4 <- glmmTMB(am ~  apd_50 + AD + CBD_50 + aspect + minpd_50 + ap + mntd50_unweigh + DBH2 , reg_sc, family=beta_family)
glm_50_P_4b <- update(glm_50_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50_P_4b)

#####
#100#####
#CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100

glm_100 <- glmmTMB(am ~ minpd_100 + avepd_100 + totpd_100 + SRA + DBH2 + CBD_100 + shannon_div_100 * RDi, reg_sc, family = beta_family)
summary(glm_100)

#画图
library(ggplot2)
coef_data <- summary(glm_100)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_100:RDi", "RDi", "shannon_div_100", "CBD_100", "DBH2", "SRA", "avepd_100", "totpd_100", "minpd_100", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "blue", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围

glm_100_M_1 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100   + (1|sptype2/resource), reg_sc, family=beta_family)
glm_100_M_1b <- update(glm_100_M_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_M_1b)
glm_100_M_2 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100  + (1|resource), reg_sc, family=beta_family)
glm_100_M_2b <- update(glm_100_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_M_2b)
glm_100_M_3 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100  + (1|sptype2), reg_sc, family=beta_family)
glm_100_M_3b <- update(glm_100_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_M_3b)
glm_100_M_4 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100 , reg_sc, family=beta_family)
glm_100_M_4b <- update(glm_100_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_M_4b)
glm_100_M_5 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + ap + tn + tp + totpd_100 + ntpd_100 + (1| Genus), reg_sc, family=beta_family)
summary(glm_100_M_5)
#CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope  
glm_100_P_1 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope + (1|sptype2/resource), reg_sc, family=beta_family)
glm_100_P_1b <- update(glm_100_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_1b)
glm_100_P_2 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope + (1|resource), reg_sc, family=beta_family)
glm_100_P_2b <- update(glm_100_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_2b)
glm_100_P_3 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope  + (1|sptype2), reg_sc, family=beta_family)
glm_100_P_3b <- update(glm_100_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_3b)
glm_100_P_4 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope + (1|sptype2), reg_sc, family=beta_family)
glm_100_P_4b <- update(glm_100_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_4b)
glm_100_P_5 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope + (1|Genus), reg_sc, family=beta_family)
glm_100_P_5b <- update(glm_100_P_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_5)
glm_100_P_6 <- glmmTMB(am ~  CBD_100 + avepd_100 + AD + aspect + ntpd_100 + DBH2 + ap + slope + (1|Type * Genus), reg_sc, family=beta_family)
glm_100_P_6b <- update(glm_100_P_6, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100_P_6b)

glm_10_1 <- glmmTMB(am ~  + (1|Type * Genus), reg_sc, family=beta_family)




#####
#em
#####

#####
#10#####

glm_10 <- glmmTMB(em ~ minpd_10 + SRA + DBH2 + CBD_10 + shannon_div_10 * RDi, reg_sc, family = beta_family)
summary(glm_10)

#画图
library(ggplot2)
coef_data <- summary(glm_10)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_10:RDi", "RDi", "shannon_div_10", "CBD_10", "DBH2", "SRA", "minpd_10", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "black", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围


#moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL
glm_10em_M_1 <- glmmTMB(em ~ moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10em_M_1b <- update(glm_10em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_10em_M_1b)
glm_10em_M_2 <- glmmTMB(em ~ moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL + (1|resource), reg_sc, family=beta_family)
glm_10em_M_2b <- update(glm_10em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_2b)
glm_10em_M_3 <- glmmTMB(em ~ moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL + (1|sptype2), reg_sc, family=beta_family)
glm_10em_M_3b <- update(glm_10em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_3b)
glm_10em_M_4 <- glmmTMB(em ~ moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL , reg_sc, family=beta_family)
glm_10em_M_4b <- update(glm_10em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_4b)
glm_10em_M_5 <- glmmTMB(em ~ moisture + tn + totpd_10 + SRA + DBH2 + mpd10_weigh + AD + soil_pc2 + SRL + (1|Type* Genus), reg_sc, family=beta_family)
glm_10em_M_5b <- update(glm_10em_M_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_M_5b)
#moisture + mntd10_weigh + mpd10_unweigh + tn + ap + shannon_div_10 + elevation
glm_10em_P_1 <- glmmTMB(em ~ moisture + mntd10_weigh + mpd10_unweigh + tn + ap + shannon_div_10 + elevation + (1|sptype2/resource), reg_sc, family=beta_family)
glm_10em_P_1b <- update(glm_10em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_1b)
glm_10em_P_2 <- glmmTMB(em ~ moisture + mntd10_weigh + mpd10_unweigh + tn + ap + shannon_div_10 + elevation + (1|resource), reg_sc, family=beta_family)
glm_10em_P_2b <- update(glm_10em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_2b)
glm_10em_P_3 <- glmmTMB(em ~ moisture + mntd10_weigh + mpd10_unweigh + tn + ap + shannon_div_10 + elevation + (1|sptype2), reg_sc, family=beta_family)
glm_10em_P_3b <- update(glm_10em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_3b)
glm_10em_P_4 <- glmmTMB(em ~ moisture + mntd10_weigh + mpd10_unweigh + tn + ap + shannon_div_10 + elevation, reg_sc, family=beta_family)
glm_10em_P_4b <- update(glm_10em_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_10em_P_4b)

#####
#20#####
#moisture + DBH2 + SRA + mpd20_weigh + ntpd_20 + slope + CBD_20 + soil_pc1 + HBD_20   

glm_20 <- glmmTMB(em ~ minpd_20 + SRA + DBH2 + CBD_20 + shannon_div_20 * RDi, reg_sc, family = beta_family)
summary(glm_20)

#画图
library(ggplot2)
coef_data <- summary(glm_20)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_20:RDi", "RDi", "shannon_div_20", "CBD_20", "DBH2", "SRA", "minpd_20", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "black", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围


glm_20em_M_1 <- glmmTMB(em ~ moisture + DBH2 + SRA + mpd20_weigh + ntpd_20 + slope + CBD_20 + soil_pc1 + HBD_20 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20em_M_1b <- update(glm_20em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_20em_M_1b)
glm_20em_M_2 <- glmmTMB(em ~  moisture + DBH2 + SRA + mpd20_weigh + ntpd_20 + slope + CBD_20 + soil_pc1 + HBD_20 + (1|resource), reg_sc, family=beta_family)
glm_20em_M_2b <- update(glm_20em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_2b)
glm_20em_M_3 <- glmmTMB(em ~  moisture + DBH2 + SRA + mpd20_weigh + ntpd_20 + slope + CBD_20 + soil_pc1 + HBD_20 + (1|sptype2), reg_sc, family=beta_family)
glm_20em_M_3b <- update(glm_20em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_3b)
glm_20em_M_4 <- glmmTMB(em ~ moisture + DBH2 + SRA + mpd20_weigh + ntpd_20 + slope + CBD_20 + soil_pc1 + HBD_20 , reg_sc, family=beta_family)
glm_20em_M_4b <- update(glm_20em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_M_4b)
#mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi 
glm_20em_P_1 <- glmmTMB(em ~  mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi + (1|sptype2/resource), reg_sc, family=beta_family)
glm_20em_P_1b <- update(glm_20em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_1b)
glm_20em_P_2 <- glmmTMB(em ~ mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi + (1|resource), reg_sc, family=beta_family)
glm_20em_P_2b <- update(glm_20em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_2b)
glm_20em_P_3 <- glmmTMB(em ~  mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi + (1|sptype2), reg_sc, family=beta_family)
glm_20em_P_3b <- update(glm_20em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_3b)
glm_20em_P_4 <- glmmTMB(em ~ mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi , reg_sc, family=beta_family)
glm_20em_P_4b <- update(glm_20em_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_4b)
glm_20em_P_5 <- glmmTMB(em ~ mntd20_weigh + moisture + tn + mpd20_unweigh + slope + soc + ap + convexity + elevation + RDi + (1|history), reg_sc, family=beta_family)
glm_20em_P_5b <- update(glm_20em_P_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_20em_P_5)

#####
#50#####

#RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50

glm_50 <- glmmTMB(em ~ minpd_50 + SRA + DBH2 + CBD_50 + shannon_div_50 * RDi, reg_sc, family = beta_family)
summary(glm_50)

#画图
library(ggplot2)
coef_data <- summary(glm_50)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_50:RDi", "RDi", "shannon_div_50", "CBD_50", "DBH2", "SRA", "minpd_50", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "black", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围

glm_50em_M_1 <- glmmTMB(em ~ RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50em_M_1b <- update(glm_50em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_50em_M_1b)
glm_50em_M_2 <- glmmTMB(em ~ RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50 + (1|resource), reg_sc, family=beta_family)
glm_50em_M_2b <- update(glm_50em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_2b)
glm_50em_M_3 <- glmmTMB(em ~ RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50 + (1|sptype2), reg_sc, family=beta_family)
glm_50em_M_3b <- update(glm_50em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_3b)
glm_50em_M_4 <- glmmTMB(em ~ RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50, reg_sc, family=beta_family)
glm_50em_M_4b <- update(glm_50em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_4b)
glm_50em_M_5 <- glmmTMB(em ~ RDi + DBH2 + tn + ntpd_50 + SRA + CBD_50 + mntd50_unweigh + mntd50_weigh + shannon_div_50 + (1|Type* Genus), reg_sc, family=beta_family)
glm_50em_M_5b <- update(glm_50em_M_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_M_5b)
#moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity
glm_50em_P_1 <- glmmTMB(em ~ moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity + (1|sptype2/resource), reg_sc, family=beta_family)
glm_50em_P_1b <- update(glm_50em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_1b)
glm_50em_P_2 <- glmmTMB(em ~ moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity + (1|resource), reg_sc, family=beta_family)
glm_50em_P_2b <- update(glm_50em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_2b)
glm_50em_P_3 <- glmmTMB(em ~ moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity + (1|sptype2), reg_sc, family=beta_family)
glm_50em_P_3b <- update(glm_50em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_3b)
glm_50em_P_4 <- glmmTMB(em ~ moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity, reg_sc, family=beta_family)
glm_50em_P_4b <- update(glm_50em_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_4b)
glm_50em_P_5 <- glmmTMB(em ~ moisture + mntd50_weigh + ntpd_50 + tn + mntd50_unweigh + slope + CBD_50 + elevation + ap + convexity + (1|Type/Genus), reg_sc, family=beta_family)
glm_50em_P_5b <- update(glm_50em_P_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_50em_P_5b)

#####
#100#####
#moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100

glm_100 <- glmmTMB(em ~ minpd_100 + SRA + DBH2 + CBD_100 + shannon_div_100 * RDi, reg_sc, family = beta_family)
summary(glm_100)

#画图
library(ggplot2)
coef_data <- summary(glm_100)$coefficients$cond
coef_data <- as.data.frame(coef_data)


# 创建一个数据框
coef_data <- data.frame(
  Variable = rownames(coef_data),
  Estimate = coef_data$Estimate,
  Std_Error = coef_data$`Std. Error`,
  Pr_value = coef_data$`Pr(>|z|)`
)

# 根据 Pr_value 创建一个新列，用于标记显著性水平
coef_data$Significance <- ifelse(coef_data$Pr_value < 0.01, "***",
                                 ifelse(coef_data$Pr_value < 0.05, "**",
                                        ifelse(coef_data$Pr_value < 0.1, "*",
                                               ifelse(coef_data$Pr_value < 0.5, ".", ""))))

# 计算线段的起点和终点
coef_data$lower <- coef_data$Estimate - coef_data$Std_Error / 2
coef_data$upper <- coef_data$Estimate + coef_data$Std_Error / 2

# 计算标记位置
coef_data$label_x <- coef_data$upper + 0.05  # 调整标记位置的 x 坐标

# 对 Variable 进行排序
desired_order <- c( "shannon_div_100:RDi", "RDi", "shannon_div_100", "CBD_100", "DBH2", "SRA", "minpd_100", "(Intercept)")
coef_data$Variable <- factor(coef_data$Variable, levels = desired_order)

# 绘制线段图和标记
ggplot(coef_data, aes(x = Estimate, y = Variable)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.1, linewidth = 0.8) + #调整线段粗细
  geom_text(aes(x = label_x, label = Significance), color = "black", hjust =0, vjust = -0.5, size = 6) + # 调整文字标记的大小
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # 添加 Estimate=0 的虚线
  labs(title = "Estimates with Significance Levels",
       x = "Estimate",
       y = "Variable") +
  theme_minimal() +
  theme(text = element_text(size = 20))# + #调整所有文字的大小
#  coord_cartesian(xlim = c(-2, 2))  # 调整 x 轴坐标范围

glm_100em_M_1 <- glmmTMB(em ~ moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_100em_M_1b <- update(glm_100em_M_1, control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(glm_100em_M_1b)
glm_100em_M_2 <- glmmTMB(em ~ moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100 + (1|resource), reg_sc, family=beta_family)
glm_100em_M_2b <- update(glm_100em_M_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_M_2b)
glm_100em_M_3 <- glmmTMB(em ~ moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100 + (1|sptype2), reg_sc, family=beta_family)
glm_100em_M_3b <- update(glm_100em_M_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_M_3b)
glm_100em_M_4 <- glmmTMB(em ~ moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100, reg_sc, family=beta_family)
glm_100em_M_4b <- update(glm_100em_M_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_M_4b)
glm_100em_M_5 <- glmmTMB(em ~ moisture + DBH2 + CBD_100 + tn + SRA + soil_pc2 + slope + totpd_100 + (1|Type/Genus), reg_sc, family=beta_family)
glm_100em_M_5b <- update(glm_100em_M_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_M_5b)
#mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2
glm_100em_P_1 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2 + (1|sptype2/resource), reg_sc, family=beta_family)
glm_100em_P_1b <- update(glm_100em_P_1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_1b)
glm_100em_P_2 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2 + (1|resource), reg_sc, family=beta_family)
glm_100em_P_2b <- update(glm_100em_P_2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_2b)
glm_100em_P_3 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2 + (1|sptype2), reg_sc, family=beta_family)
glm_100em_P_3b <- update(glm_100em_P_3, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_3b)
glm_100em_P_4 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2, reg_sc, family=beta_family)
glm_100em_P_4b <- update(glm_100em_P_4, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_4b)
glm_100em_P_5 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2 + (1|Type/Genus), reg_sc, family=beta_family)
glm_100em_P_5b <- update(glm_100em_P_5, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_5b)
glm_100em_P_6 <- glmmTMB(em ~ mntd100_weigh + moisture + tn + CBD_100 + slope +  convexity + ap + DBH2 + (1|Genus), reg_sc, family=beta_family)
glm_100em_P_6b <- update(glm_100em_P_6, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(glm_100em_P_6)


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