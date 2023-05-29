
###先加载数据

load("E:/黑石顶测菌根/菌根侵染率/数据整理/tmp/For_git_Rstudio/root_qrl_soil.RData")


###计算一下乔木生长速率
d$gr_rate <- (log(d$DBH2)-log(d$DBH1))

####筛选出AM的个体
d_AM <- subset(d, !is.na(qr_AM))
#对样本总体分析
#####
#1.相关性分析

####进行boxcox转换
library(MASS)
# 对d_AM$qr_AM列使用Box-Cox转换
boxcox(qr_AM ~ 1, data = d_AM,lambda = seq(4,5, length.out = 10))

#画出转换后的形状
hist(d_AM$qr_AM^(4.6))
#进行lm
r2=lm(gr_rate~qr_AM^4.6, data=d_AM)
#画出残差的分布
hist(r2$residuals)
summary(r2)

#画出一列变换后的数据，方便画图
d_AM$AM_4.6 = (d_AM$qr_AM)^4.6

###使用ggplot2画图
plot(d_AM$gr_rate~d_AM$AM_4.6)

####试一下结果的变化
cor.test(d_AM$qr_AM, d_AM$gr_rate)
cor.test(d_AM$AM_4.6, d_AM$gr_rate)


# 拟合二项式回归模型

fit <- glm(sc_gr_rate ~ sc_qr_AM, data = d_AM, family = binomial(link = "logit"))

# 绘制散点图和拟合曲线
library(ggplot2)
ggplot(d, aes(x = qr_AM, y = gr_rate)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +
  ggtitle(paste0("Logistic Regression (AIC = ", round(AIC(fit), 2), ")"))



######
#2.主成分分析
######
# 画散点图
ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2])) +
  geom_point(aes(x=PC1, y=PC2))

#####
#3.线性回归分析
#####
# 画散点图和回归线
# calculate slope and intercept
m <- coef(lm(sc_gr_rate ~ sc_qr_AM, data = d_AM))[2]
b <- coef(lm(sc_gr_rate ~ sc_qr_AM, data = d_AM))[1]

p <- ggplot(d_AM, aes(x=sc_qr_AM, y=sc_gr_rate)) + 
  geom_point(color="#e5f5f9", shape=16, size=3) + 
  geom_smooth(method="lm", se=TRUE, color="#99d8c9", size=1.2) +
  xlab("AM侵染率") + 
  ylab("生长速率") +
  theme_minimal()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 18))

# calculate correlation coefficient and p-value
r <- cor(d_AM$sc_qr_AM, d_AM$sc_gr_rate)
pval <- summary(lm(sc_gr_rate ~ sc_qr_AM, data = d_AM))$coef[2, 4]

# add annotation for correlation coefficient and p-value
p <- p +
  annotate("text", x = max(d_AM$sc_qr_AM), y = min(d_AM$sc_gr_rate), 
           label = paste0("R = ", round(r, 2), "\n p = ", format(pval, scientific = TRUE, digits = 2)),
           hjust = 1) +
  geom_text(x = max(d_AM$sc_qr_AM), y = max(d_AM$sc_gr_rate), 
            label = paste0("Slope = ", round(m, 2), "\n Intercept = ", round(b, 2)), 
            hjust = 1, vjust = 1)

p

#4.方差分析
#####
# 画箱线图
ggplot(d, aes(x=Species.x, y=qr_AM)) + 
  geom_boxplot() +
  xlab("Species.x") +
  ylab("qr_AM") +
  theme_minimal()

# 进行多重比较
TukeyHSD(model)

#对不同物种
#########
# 使用ggplot2包绘制箱线图，以乔木树种作为分组变量
library(ggplot2)

ggplot(d, aes(x=Species.x, y=gr_rate, fill=物种)) + 
  geom_boxplot() +
  xlab("Species.x") +
  ylab("gr_rate") +
  theme_minimal()

# 使用facet_grid()函数绘制多图合并，以乔木树种作为分组变量
ggplot(d, aes(x=qr_AM, y=gr_rate)) + 
  geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  facet_grid(~物种) +
  xlab("qr_AM") + 
  ylab("gr_rate") +
  theme_minimal()

# 使用tapply()函数计算不同物种的平均菌根侵染率
mean_by_species <- tapply(d$qr_AM, d$物种, mean)

# 使用barplot()函数绘制不同物种的平均菌根侵染率柱状图
barplot(mean_by_species, 
        main="不同物种的平均菌根侵染率", 
        xlab="物种", 
        ylab="平均菌根侵染率")

