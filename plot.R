
#pca
library(ggplot2)

# 示例数据框
data <- data.frame(
  Sample = c("Sample1", "Sample2", "Sample3", "Sample4"),
  PC1 = c(0.1, -0.2, 0.3, -0.4),
  PC2 = c(0.4, -0.3, 0.2, -0.1),
  SizeVariable = c(10, 15, 20, 25)  # 用于控制点大小的额外变量
)

pca_data <- reg[, c("soc","tn","tp","ap","ph","moisture")]
pca_result <- prcomp(pca_data, scale = TRUE)
pca_result <- pca_result$x
pca_results <- data.frame(PC1 = pca_result[, 1], PC2 = pca_result[, 2], reg[,c("am","sptype2")])
pca_results$qr_AM <- ifelse(pca_results$qr_AM < 0.25, 0.125,
                     ifelse(pca_results$qr_AM < 0.5, 0.375,
                     ifelse(pca_results$qr_AM < 0.75, 0.625, 0.875)))
pca_results$qr_AM <- pca_results$qr_AM/10
# 创建PCA/PCoA图，并根据SizeVariable调整点的大小
pca_plot <- ggplot(pca_results, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sptype2)) +  # size = am,根据SizeVariable调整点的大小
  labs(x = "PC1", y = "PC2") +            # 添加坐标标签
  ggtitle("PCA Plot")+                # 添加标题
  scale_color_manual(values = c("dominant" = "#cb181d", "common" = "#9e9ac8", "rare" = "#fd8d3c"))#M_Type

# 打印图形
print(pca_plot)


library(ggplot2)
###一次
# 使用 ggplot2 创建散点图
ggplot(reg, aes(x = am, y = tn, color = resource)) +
  geom_point() +
  facet_grid(. ~ sptype2) +
  labs(x = "am", y = "tn") +
  ggtitle("关系图")





###二次
# 创建散点图，执行二次拟合
# 创建散点图，执行二次拟合

# 选择数据并删除包含 NA 的行
reg_am <- reg %>%  
  filter(!is.na(am)) 
  
reg_am$sptype2 <- factor(reg_am$sptype2, levels = c("dominant", "common", "rare"))

# 创建散点图，执行二次拟合
ggplot(reg_am, aes(x = mntd20_unweigh, y = am, color = resource)) +
  geom_point() +
  facet_grid(. ~ reorder(sptype2, c("dominant", "common", "rare"))) +
  labs(x = "mntd20_unweigh", y = "am") +
  ggtitle("am-mntd20_unweigh") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  theme(text = element_text(size = 24))

reg_sc$sptype2 <- factor(reg_sc$sptype2, levels = c("dominant", "common", "rare"))

summary(glm_20b)
ggplot(reg_am, aes(x = tn, y = am, color = resource)) +
  geom_point() +
  facet_grid(. ~ sptype2) +
  labs(x = "tn", y = "am") +
  ggtitle("am-tn") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(text = element_text(size = 20))


####画大图
# 创建 ggplot 对象并指定数据和映射
library(ggplot2)

# 创建 ggplot 对象并指定数据和映射
ggplot(reg_am, aes(x = RDi, y = am, color = sptype2)) +
  
  # 添加散点图层
  geom_point() +
  
  # 自定义颜色映射
  scale_color_manual(values = c("common" = "#9e9ac8", "dominant" = "#cb181d", "rare" = "#fd8d3c")) +
  
  # 添加 x 和 y 轴标签
  labs(x = "RDi", y = "am") +
  
  # 添加图表标题
  ggtitle("RDi vs. am") +
  
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  
  # 添加所有点的拟合曲线
  geom_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2) ,color = "black", linetype = "dashed")

# 拟合 "common" sptype2 的线
model_common <- glm(am ~  RDi, data = subset(reg_am, sptype2 == "common"), family = quasibinomial(link = "logit"))
# 拟合 "dominant" sptype2 的线
model_dominant <- glm(am ~  RDi, data = subset(reg_am, sptype2 == "dominant"), family = quasibinomial(link = "logit"))

# 拟合 "rare" sptype2 的线
model_rare <- glm(am ~  RDi, data = subset(reg_am, sptype2 == "rare"), family = quasibinomial(link = "logit"))


# 拟合所有数据的总线
model_total <- glm(am ~  RDi, data = reg_am, family = quasibinomial(link = "logit"))


