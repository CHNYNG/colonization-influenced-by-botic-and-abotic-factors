# for the results
##### caculate the colnization rate of each species ####
library(dplyr)
Species <- reg %>%
  group_by(Latin) %>%
  summarize(mean_am = mean(qr_AM, na.rm = TRUE),
            sd_am = sd(qr_AM, na.rm = TRUE)^2) %>%
  arrange(desc(mean_am))
###### 做密度图 ####
library(dplyr)
library(ggplot2)
species_rate <- reg %>%
  subset(, select = c(Latin, qr_AM))
ggplot(species_rate, aes(x = qr_AM, fill = Latin)) +
  geom_density(alpha = 0.3, show.legend = FALSE) +  # 隐藏图例以简化视觉
  theme_minimal() +
  labs(title = "Density Plot of Colonization Rates for Multiple Species",
       x = "Colonization Rate", y = "Density") +
  scale_fill_viridis_d()  # 使用渐变调色板

##### 做热图 ####
install.packages("pheatmap")
library(pheatmap)

# 创建一个数据矩阵，假设列是侵染率、标准差，行是物种
# 你可以从原始数据中提取并生成类似矩阵

data_matrix <- as.matrix(cbind(Species$mean_am,Species$sd_am))
rownames(data_matrix) <- Species$Latin
colnames(data_matrix) <- c("ColonizationRate", "SD") 


# 使用pheatmap生成热图
p <- pheatmap(data_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(50),  # 定义颜色渐变
         cluster_rows = TRUE,   # 按行聚类
         cluster_cols = FALSE,  # 不按列聚类
         main = "Heatmap of Colonization Rates for 161 Species")

# 为了避免breaks值重复问题，可以微调两个seq之间的边界
# 使用na.rm = TRUE 来处理缺失值
my_breaks <- c(seq(min(data_matrix, na.rm = TRUE), 
                   quantile(data_matrix, 0.75, na.rm = TRUE), length.out = 30),  # 前75%部分
               seq(quantile(data_matrix, 0.75, na.rm = TRUE) + 1e-6,  # 添加一个微小值避免重复
                   max(data_matrix, na.rm = TRUE), length.out = 20))  # 后25%部分

# 生成热图
p <- pheatmap(data_matrix, 
              color = colorRampPalette(c("blue", "white", "red"))(50),  
              breaks = my_breaks,  # 使用自定义的颜色区间
              cluster_rows = TRUE,  
              cluster_cols = FALSE,  
              main = "Heatmap with Enhanced Red Area for High Values")
install.packages("ComplexHeatmap")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
Heatmap(data_matrix, 
        name = "ColonizationRate",  # 自定义图例名称
        heatmap_legend_param = list(
          title = "ColonizationRate",   # 图例的标题
          title_gp = gpar(fontsize = 12), # 图例标题字体大小
          labels_gp = gpar(fontsize = 10), # 图例标签字体大小
          direction = "horizontal"  # 设置图例的方向为横向
        ),
        show_column_names = TRUE,
        show_row_names = TRUE)
ggsave(
  "pic/heat.png",
  p,
  width = 2000,
  height = 8000,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

#总结物种的组间存在显著差异
kruskal_test <- kruskal.test(am ~ Latin, data = reg)
kruskal_test
#两两的组间比价
model <- lm(am ~ Latin, data = reg)
tukey_result <- TukeyHSD(aov(model))
# 筛选出 p < 0.05 的部分
significant_pairs <- tukey_result$Latin[tukey_result$Latin[, "p adj"] < 0.05, ]
# 筛选出 p = 1 的部分
non_significant_pairs <- tukey_result$Latin[tukey_result$Latin[, "p adj"] == 1, ]



Genus <- reg %>%
  group_by(Genus)%>%
  summarize(mean_am = mean(am, na.rm = TRUE),
            sd_am = sd(am, na.rm = TRUE)^2) %>%
  arrange(desc(mean_am))
#总结物种的组间存在显著差异
kruskal_test <- kruskal.test(am ~ Genus, data = reg)
kruskal_test
#两两的组间比价
model <- lm(am ~ Genus, data = reg)
tukey_result <- TukeyHSD(aov(model))
# 筛选出 p < 0.05 的部分
significant_pairs <- tukey_result$Genus[tukey_result$Genus[, "p adj"] < 0.05, ]
# 筛选出 p = 1 的部分
non_significant_pairs <- tukey_result$Genus[tukey_result$Genus[, "p adj"] == 1, ]

Family <- reg %>%
  group_by(Family)%>%
  summarize(mean_am = mean(am, na.rm = TRUE),
            sd_am = sd(am, na.rm = TRUE)^2) %>%
  arrange(desc(mean_am))
