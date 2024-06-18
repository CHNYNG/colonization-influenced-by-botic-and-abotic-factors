library(ggplot2)


elev <- read.csv("data/HSD+env5m.csv", header = TRUE, fileEncoding = "GBK")
ggplot(elev, aes(x = gx, y = gy, z = elev)) +
  geom_contour() +
  labs(title = "Contour Plot with Height Information",
       x = "GX",
       y = "GY",
       z = "Height")


ggplot(data = reg, aes(x = GX, y = GY)) +
  geom_point(color = "lightblue", size = 5, shape = 16) +  # 浅色大点
  geom_point(color = "darkblue", size = 1) +  # 深色小点
  theme_minimal() +  # 设置主题为简洁风格
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5)
  ) +
  labs(
    title = "Sample Sites",
    x = "GX",
    y = "GY"
  )+
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 500))  # 限制坐标轴范围

library(ggplot2)

plot <- ggplot() +
  geom_contour(data = elev, aes(x = gx, y = gy, z = elev), color = "#636363", linewidth = 1, show.legend = FALSE) +  # 设置等高线颜色为黑色
  geom_point(data = reg, aes(x = GX, y = GY, size = am*0.01 ), shape = 16, color = "#636380") +  # 根据 am 大小确定点的大小，根据 Latin 分类确定颜色
  labs(x = "GX",
       y = "GY",
       z = "Height") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 16, hjust = 0.5),
#    plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5)
    legend.position = "none"
  ) +
  scale_size_continuous(range = c(1, 5)) +  # 设置点的大小范围
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 500))  # 限制坐标轴范围
print(plot)

library(ggplot2)
combined_plot <- ggplot() +
  geom_contour_filled(data = elev, aes(x = gx, y = gy, z = elev), color = "#636363", linewidth = 1, show.legend = FALSE)+  # 设置等高线颜色为#636363
#  scale_fill_gradientn(colours = terrain.colors(14))+
#  scale_fill_distiller(palette = "Spectral") +
  discrete_scale('fill', 'myscale', colorRampPalette(c("#440154", "#481d6f","#453581","#3d4d8a","#34618d", "#2b748e", "#24878e","#1f998a","#25ac82","#40bc72","#63cd54", "#97d83f","#d9f0a3","#f7fcb9"))) +
    geom_point(data = reg, aes(x = GX, y = GY), color = "#bcbddc", size = 1, shape = 16) +  # 设置散点颜色为浅蓝色，大小为3，形状为实心圆
#  geom_point(data = reg, aes(x = GX, y = GY, size = am * 0.001), color = "#a6bddb", ) +  # 设置散点颜色为深蓝色，大小为1
  labs(x = "GX",
       y = "GY",
       z = "Height") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 24),
    axis.title = element_text(color = "black", size = 24),
    plot.title = element_text(color = "black", size = 24, hjust = 2),
    plot.subtitle = element_text(color = "black", size = 24, hjust = 2)
  ) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 500))  # 限制坐标轴范围

print(combined_plot)

library(dplyr)

# 计算 reg$am 的均值并按其排序
Latin_mean <- reg %>%
  group_by(Latin) %>%
  summarize(mean_am = mean(am, na.rm = TRUE)) %>%
  arrange(desc(mean_am))


# 按 reg$Latin 的均值排序后绘制柱状图
top_6 <- Latin_mean %>%
  arrange(desc(mean_am)) %>%
  head(6) %>%
  pull(Latin)

bottom_6 <- Latin_mean %>%
  arrange(mean_am) %>%
  head(6) %>%
  pull(Latin)

# Combine top 5 and bottom 5 species
selected_species <- c(top_6, bottom_6)

# Filter the original data to include only these species
filtered_reg <- reg %>%
  mutate(group = case_when(
    Latin %in% top_6 ~ "Top 5",
    Latin %in% bottom_6 ~ "Bottom 5",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))
filtered_reg <- filtered_reg %>%
  filter(Latin %in% selected_species)
filtered_reg <- filtered_reg %>%
  mutate(Latin = factor(Latin, levels = c(top_6, bottom_6)))


# Create the violin plot
ggplot(filtered_reg, aes(x = reorder(Latin, -am), y = am, fill = Latin)) +
  geom_violin() +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +  # Add jitter for better visualization of points
  labs(x = "Species", y = "Colonization Intensity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic', size = 8),
    panel.grid.major = element_line(color = "gray"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5),
    legend.position = "none"  # Remove legend if not needed
  ) +
  coord_cartesian(ylim = c(0.3, 1))  # Set y-axis limits

