
#pca
library(ggplot2)

# 示例数据框
data <- data.frame(
  Sample = c("Sample1", "Sample2", "Sample3", "Sample4"),
  PC1 = c(0.1, -0.2, 0.3, -0.4),
  PC2 = c(0.4, -0.3, 0.2, -0.1),
  SizeVariable = c(10, 15, 20, 25)  # 用于控制点大小的额外变量
)
pca_data <- reg[, c("soc","tn","tp","ap","ph")]
pca_result <- prcomp(pca_data, scale = TRUE)
pca_result <- pca_result$x
pca_results <- data.frame(PC1 = pca_result[, 1], PC2 = pca_result[, 2], reg[,c("qr_AM","M_Type")])
pca_results$qr_AM <- ifelse(pca_results$qr_AM < 0.25, 0.125,
                     ifelse(pca_results$qr_AM < 0.5, 0.375,
                     ifelse(pca_results$qr_AM < 0.75, 0.625, 0.875)))
pca_results$qr_AM <- pca_results$qr_AM/10
# 创建PCA/PCoA图，并根据SizeVariable调整点的大小
pca_plot <- ggplot(pca_results, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = M_Type)) +  # 根据SizeVariable调整点的大小
  labs(x = "PC1", y = "PC2") +            # 添加坐标标签
  ggtitle("PCA Plot")+                # 添加标题
  scale_color_manual(values = c("AM" = "#ffffcc", "AM-EM" = "#a1dab4", "EM" = "#41b6c4", "NA" = "#225ea8"))

# 打印图形
print(pca_plot)
