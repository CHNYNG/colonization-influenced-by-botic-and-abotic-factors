library(betareg)
library(ggplot2)
library(dplyr)
library(tidyr)
#############selet tn tp and tn/tp and colonization rate############ 

df <- reg_scaled %>%
  dplyr::select(tn, tp, qr_AM, soc, ap, ph, moisture, tndtp, Latin) %>%       # 选择需要的列
  tidyr::drop_na() 

# 查看处理后的数据框
head(df)

##### tn #########
am_tn <- betareg(qr_AM ~ tn, data = df)
summary(am_tn) #0.882
p_value <- coef(summary(am_tn))$mean["tn", "Pr(>|z|)"]
tn_coef <- coef(summary(am_tn))$mean["tn", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", tn_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#tn fitted values
am_tn_fitted <- predict(am_tn, type = "response", na.action = na.exclude)

am_tn_dat <- df %>% 
  dplyr::select(qr_AM, tn, Latin) %>% 
  dplyr::mutate(
    fitted = am_tn_fitted)

am_tn_p <- ggplot(am_tn_dat, aes(x = tn, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "TN", y = "AMCR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tn), y = max(df$qr_AM), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tn_p)

##### tp #########
am_tp <- betareg(qr_AM ~ tp, data = df)
summary(am_tp) #0.882
p_value <- coef(summary(am_tp))$mean["tp", "Pr(>|z|)"]
tp_coef <- coef(summary(am_tp))$mean["tp", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", tp_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#tp fitted values
am_tp_fitted <- predict(am_tp, type = "response", na.action = na.exclude)

am_tp_dat <- df %>% 
  dplyr::select(qr_AM, tp, Latin) %>% 
  dplyr::mutate(
    fitted = am_tp_fitted)

am_tp_p <- ggplot(am_tp_dat, aes(x = tp, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "TP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tp), y = max(df$qr_AM), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tp_p)

#### tn/tp ####
df <- df[!is.na(df$tndtp) & df$tndtp != max(df$tndtp, na.rm = TRUE), ]
am_tndtp <- betareg(qr_AM ~ tndtp, data = df)
summary(am_tndtp) #0.882
p_value <- coef(summary(am_tndtp))$mean["tndtp", "Pr(>|z|)"]
tndtp_coef <- coef(summary(am_tndtp))$mean["tndtp", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", tndtp_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#tndtp fitted values
am_tndtp_fitted <- predict(am_tndtp, type = "response", na.action = na.exclude)

am_tndtp_dat <- df %>% 
  dplyr::select(qr_AM, tndtp, Latin) %>% 
  dplyr::mutate(
    fitted = am_tndtp_fitted)

am_tndtp_p <- ggplot(am_tndtp_dat, aes(x = tndtp, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "TN/TP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tndtp), y = max(df$qr_AM), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tndtp_p)



##### soc #########
am_soc <- betareg(qr_AM ~ soc, data = df)
summary(am_soc) #0.882
p_value <- coef(summary(am_soc))$mean["soc", "Pr(>|z|)"]
soc_coef <- coef(summary(am_soc))$mean["soc", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", soc_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#soc fitted values
am_soc_fitted <- predict(am_soc, type = "response", na.action = na.exclude)

am_soc_dat <- df %>% 
  dplyr::select(qr_AM, soc, Latin) %>% 
  dplyr::mutate(
    fitted = am_soc_fitted)

am_soc_p <- ggplot(am_soc_dat, aes(x = soc, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "SOC", y = "AMCR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$soc), y = max(df$qr_AM), 
           label =  sprintf("italic(beta) == '%0.3f' ~ italic(p)=='%0.3f'", soc_coef, p_value), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_soc_p)

### ap ####

am_ap <- betareg(qr_AM ~ ap, data = df)
summary(am_ap) #0.882
p_value <- coef(summary(am_ap))$mean["ap", "Pr(>|z|)"]
ap_coef <- coef(summary(am_ap))$mean["ap", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", ap_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#ap fitted values
am_ap_fitted <- predict(am_ap, type = "response", na.action = na.exclude)

am_ap_dat <- df %>% 
  dplyr::select(qr_AM, ap, Latin) %>% 
  dplyr::mutate(
    fitted = am_ap_fitted)

am_ap_p <- ggplot(am_ap_dat, aes(x = ap, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "AP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$ap), y = max(df$qr_AM), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_ap_p)

### ph ####

am_ph <- betareg(qr_AM ~ ph, data = df)
summary(am_ph) #0.882
p_value <- coef(summary(am_ph))$mean["ph", "Pr(>|z|)"]
ph_coef <- coef(summary(am_ph))$mean["ph", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", ph_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#ph fitted values
am_ph_fitted <- predict(am_ph, type = "response", na.action = na.exclude)

am_ph_dat <- df %>% 
  dplyr::select(qr_AM, ph, Latin) %>% 
  dplyr::mutate(
    fitted = am_ph_fitted)

am_ph_p <- ggplot(am_ph_dat, aes(x = ph, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "pH", y = "AMCR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$ph), y = max(df$qr_AM), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_ph_p)

### moisture ####

am_moisture <- betareg(qr_AM ~ moisture, data = df)
summary(am_moisture) #0.882
p_value <- coef(summary(am_moisture))$mean["moisture", "Pr(>|z|)"]
moisture_coef <- coef(summary(am_moisture))$mean["moisture", "Estimate"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# 设置系数标签
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", moisture_coef))

# 根据 p 值调整线型
line_type <- ifelse(p_value > 0.01, "dashed", "solid")
#moisture fitted values
am_moisture_fitted <- predict(am_moisture, type = "response", na.action = na.exclude)

am_moisture_dat <- df %>% 
  dplyr::select(qr_AM, moisture, Latin) %>% 
  dplyr::mutate(
    fitted = am_moisture_fitted)

am_moisture_p <- ggplot(am_moisture_dat, aes(x = moisture, y = qr_AM)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "moisture", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$moisture), y = max(df$qr_AM), 
           label =  sprintf("italic(beta) == '%0.3f' ~ italic(p)=='%0.3f'", moisture_coef, p_value), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_moisture_p)

#### output the pictures ####

library(cowplot)
lmplots2 <- plot_grid(
  plot_grid(am_SRA_p, am_SRL_p, am_AD_p, ncol = 3),   # 第一行3个图
  plot_grid(am_BD_p,am_CBD_p, am_invsimpson_div_p, ncol = 3),   # 第二行3个图
  plot_grid(am_minpd_p, am_avepd_p, am_totpd_p,ncol = 3),        # 第三行3个图
  plot_grid(am_pcoa1_p,am_pcoa2_p ,ncol = 3), 
  plot_grid(am_dbh_p,am_RD_p ,NULL,  ncol = 3), 
  plot_grid(am_tn_p,am_tp_p, am_tndtp_p,  ncol = 3),
  plot_grid(am_soc_p, am_ap_p, NULL, ncol = 3),
  plot_grid(am_ph_p, am_moisture_p, NULL, ncol = 3),
  ncol = 1  # 总体按列排列
)+ 
  theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))

#lmplots2

ggsave(
    "pic/figure S2_2.png",
    lmplots2,
    width = 4000,
    height = 7200,
    units = "px",
    dpi = 300,
    bg = "#FFFFFF"
   )

ggsave(
  "pic/tn.png",
  am_tn_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/tp.png",
  am_tp_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/tndtp.png",
  am_tndtp_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/soc.png",
  am_soc_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/ap.png",
  am_ap_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/ph.png",
  am_ph_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)

ggsave(
  "pic/moisture.png",
  am_moisture_p,
  width = 1500,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)
##### for the data ring ####

####### first try #######
library(ape)
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
library(dplyr)
#BiocManager::install("ggtreeExtra")
library(ggtreeExtra)
# in put the tree
phylo_tree <- hsd_phytr
phylo_tree$tip.label <- gsub("_", " ", phylo_tree$tip.label)

colonization_data <- reg %>%
  select(Latin, am) %>%
  rename(species = Latin,
         colonization_rate = am)
colonization_data$species <- gsub("_", " ", colonization_data$species)
colonization_mean <- colonization_data %>%
  group_by(species) %>%
  summarize(colonization_rate = mean(colonization_rate, na.rm = TRUE)) %>%
  as.data.frame() %>%
  drop_na()


# 按照系统发育树的顺序重新排列 colonization_mean$species
colonization_mean <- colonization_mean[match(selected_phylo_tree$tip.label, colonization_mean$species), ]

# 现在 colonization_mean 数据框的物种顺序应该与 selected_phylo_tree$tip.label 相同了
selected_phylo_tree <- keep.tip(phylo_tree, tip = intersect(phylo_tree$tip.label, colonization_mean$species))

p <- ggtree(selected_phylo_tree, layout = "circular") + 
  geom_fruit(
    data = colonization_data, 
    geom = geom_bar, 
    mapping = aes(y = species, x = colonization_rate, fill = colonization_rate), 
    orientation = "y", 
    stat = "identity", 
    offset = 0,   # 增大柱状图偏移，使其占据第二圈
    width = 0.15    # 调整柱状图的宽度
  ) +
  geom_tiplab(aes(label = label), align = TRUE, linesize = 0.5, offset = 50, fontface = "italic") +  # 增大物种名的偏移到最外圈
  scale_fill_gradient2(low = "#99ff4d", mid = "#22c32e", high = "#00a15c", midpoint = 0) +
  theme(legend.position = "none")
p

ggsave(
  "pic/figure S1.png",
  p,
  width = 6000,
  height = 6000,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)


##### second try #####

library(phytools)
library(ape)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(dplyr)
library(tidyr)

#### 看下还能不能和phylo tree match上####
selected_phylo_tree$node.label <- NULL

# 找出热点区域，假设高于某个阈值为热点
# 可以根据实际数据调整阈值
library(caper)

comp_data <- comparative.data(selected_phylo_tree, colonization_mean, species)

pgls_model <- pgls(colonization_rate ~ 1, data = comp_data)

# 获取残差
residuals <- residuals(pgls_model)

# 判断残差显著性，用于标记hotspot物种
hotspot_threshold <- 1 * sd(residuals)  # 设定阈值为2个标准差

hotspots <- ifelse(abs(residuals) > hotspot_threshold, "hotspot", "non-hotspot")
colonization_mean$hotspots <- hotspots

# 绘制基础系统发育树
p <- ggtree(selected_phylo_tree, layout = "circular") +
  theme(legend.position = "none") +
  # 绘制侵染率的图，并根据hotspots调整颜色
  geom_fruit(
    data = colonization_mean %>% 
      mutate(fillcol = ifelse(hotspots == "hotspot", 
                             "#feb24c", "#004529" )),
    geom = geom_bar,
    mapping = aes(y = species, x = colonization_rate, 
                  fill = fillcol),
    orientation = "y",
    stat = "identity",
    offset = 0,   # 增大柱状图偏移，使其占据第二圈
    width = 0.3
  ) +
  geom_tiplab(aes(label = label), align = FALSE, linesize = 0.5, offset = 60, fontface = "italic") +  # 增大物种名的偏移到最外圈
  scale_fill_identity()
  # scale_fill_gradient2(low = "#99ff4d", mid = "#22c32e", high = "#00a15c", midpoint = 0, na.value = "gray") +
 # scale_fill_manual(values = c("hotspot" = "red", "non-hotspot" = "#00a15c"))

# 添加侵染率的可视化
print(p)

ggsave(
  "pic/figure S1.png",
  p,
  width = 6000,
  height = 6000,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)
