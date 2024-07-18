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

ggsave(
  "pic/figure1.png",
  combined_plot,
  width = 3000,
  height = 1600,
  units = "px",
  dpi = 300
  
)
library(dplyr)
library(ggplot2)
library(ggpubr)
# 计算 reg$am 的均值并按其排序
Latin_mean <- reg %>%
  group_by(Latin) %>%
  summarize(mean_am = mean(am, na.rm = TRUE)) %>%
  arrange(desc(mean_am))
Latin_mean <- reg %>%
  group_by(Latin) %>%
  summarize(
    mean_am = mean(am, na.rm = TRUE),
    sd_am = sd(am, na.rm = TRUE)
  )
write.csv(Latin_mean, "Appendix.Table S1.csv")
# 手动选择要显示的物种
# 选定的物种名称
selected_species <- c("Exbucklandia_tonkinensis", "Garcinia_multiflora", "Ardisia_elegans", 
                      "Cryptocarya_concinna", "Vitex_quinata", "Rhododendron_simsii", 
                      "Castanopsis_fissa", "Ixonanthes_chinensis", "Craibiodendron_stellatum", 
                      "Lithocarpus_lohangwu")

# 设置颜色
species_colors <- c("Exbucklandia tonkinensis" = "#004529", 
                    "Garcinia multiflora" = "#004529", 
                    "Ardisia elegans" = "#004529", 
                    "Cryptocarya concinna" = "#004529", 
                    "Vitex quinata" = "#004529", 
                    "Rhododendron simsii" = "#addd8e", 
                    "Castanopsis fissa" = "#addd8e", 
                    "Ixonanthes chinensis" = "#b2df8a", 
                    "Craibiodendron stellatum" = "#addd8e", 
                    "Lithocarpus lohangwu" = "#addd8e")

# 构建数据框，并将拉丁名中的下划线替换为空格
filtered_reg <- reg %>%
  filter(Latin %in% selected_species) %>%
  mutate(Latin = gsub("_", " ", Latin)) %>%
  mutate(Latin = factor(Latin, levels = gsub("_", " ", selected_species)))

am_stats <- filtered_reg %>%
  group_by(Latin) %>%
  summarize(
    mean_am = mean(am, na.rm = TRUE),
    sd_am = sd(am, na.rm = TRUE)
)

am.plot <- ggplot(am_stats, aes(x = Latin, y = mean_am, fill = Latin)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_am - sd_am, ymax = mean_am + sd_am), 
                width = 0.2, position = position_dodge(0.9)) +
#  geom_point(data = filtered_reg, aes(x = Latin, y = am), 
#             position = position_jitter(width = 0.2), alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(
    c("Exbucklandia tonkinensis", "Rhododendron simsii"),
    c("Garcinia multiflora","Castanopsis fissa" ),
    c("Ardisia elegans","Ixonanthes chinensis" ),
    c("Cryptocarya concinna", "Craibiodendron stellatum"),
    c("Vitex quinata", "Lithocarpus lohangwu")
  )) +
  labs(x = NULL, y = "Colonization Rates") +
  scale_fill_manual(values = species_colors) +  # 手动设置颜色
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic', size = 14),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_line(color = "gray"),
    axis.title = element_text(color = "black", size = 14),
    plot.title = element_text(color = "black", size = 16, hjust = 0.5),
    plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5),
    legend.position = "none"  # 移除图例（如果不需要）
  ) +
  coord_cartesian(ylim = c(0.3, 1))  # 设置 y 轴范围

print(am.plot)


ggsave(
  "pic/figure2.png",
  am.plot,
  width = 3000,
  height = 2200,
  units = "px",
  dpi = 300,
  bg = "#ffffff"
)

# 每个变量单独作图 ----
# DBH2
am_dbh <- betareg(am ~ DBH2, data = am_beta_dat)
summary(am_dbh) #显著，留下
p_value <- coef(summary(am_dbh))$mean["DBH2", "Pr(>|z|)"]
# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

#dbh fitted values
am_dbh_fitted <- predict(am_dbh, type = "response")

am_dbh_dat <- am_beta_dat %>% 
  select(am, DBH2, Latin) %>% 
  mutate(
    fitted = am_dbh_fitted)

am_dbh_p <- ggplot(am_dbh_dat, aes(x = DBH2, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "DBH", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(am_dbh_dat$DBH2), y = max(am_dbh_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_dbh_p)


# Add fitted values
am_dbh_fitted <- predict(am_dbh, type = "response")
am_dbh_lower <- predict(am_dbh, type = "quantile", at = 0.025)
am_dbh_upper <- predict(am_dbh, type = "quantile", at = 0.975)

am_dbh_dat <- am_beta_dat %>% 
  select(am, DBH2, Latin) %>% 
  mutate(
    fitted = am_dbh_fitted,
    lower = am_dbh_lower,
    upper = am_dbh_upper
  )

am_dbh_p <- ggplot(am_dbh_dat, aes(x = DBH2, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "DBH", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_dbh_p

###做其他变量的图
###### SRA
am_beta_dat_SRA <- am_beta_dat[!is.na(am_beta_dat$SRA), ]
am_SRA <- betareg(am ~ SRA, data = am_beta_dat_SRA)
summary(am_SRA) #p = 0.0881
# Add fitted values
am_SRA_fitted <- predict(am_SRA, type = "response")

am_SRA_dat <- am_beta_dat_SRA %>% 
  select(am, SRA, Latin) %>% 
  mutate(
    fitted = am_SRA_fitted)

am_SRA_p <- ggplot(am_SRA_dat, aes(x = SRA, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "SRA", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_SRA_p

###### SRL
am_beta_dat_SRL <- am_beta_dat[!is.na(am_beta_dat$SRL), ]
am_SRL <- betareg(am ~ SRL, data = am_beta_dat_SRL)
summary(am_SRL) #不显著，丢掉
# Add fitted values
am_SRL_fitted <- predict(am_SRL, type = "response")

am_SRL_dat <- am_beta_dat_SRL %>% 
  select(am, SRL, Latin) %>% 
  mutate(
    fitted = am_SRL_fitted)

am_SRL_p <- ggplot(am_SRL_dat, aes(x = SRL, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "SRL", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_SRL_p


#AD
am_beta_dat_AD <- am_beta_dat[!is.na(am_beta_dat$AD), ]
am_AD <- betareg(am ~ AD, data = am_beta_dat_AD)
summary(am_AD) #p=0.0294

p_value <- coef(summary(am_AD))$mean["AD", "Pr(>|z|)"]

# 设置 p 值标签，保留三位小数
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# Add fitted values
am_AD_fitted <- predict(am_AD, type = "response")

am_AD_dat <- am_beta_dat_AD %>% 
  select(am, AD, Latin) %>% 
  mutate(
    fitted = am_AD_fitted)

am_AD_p <- ggplot(am_AD_dat, aes(x = AD, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "AD", y =NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )+
  annotate("text", x = max(am_AD_dat$AD), y = max(am_AD_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)
am_AD_p


#BD
am_beta_dat_BD <- am_beta_dat[!is.na(am_beta_dat$BD_10), ]
am_BD <- betareg(am ~ BD_10, data = am_beta_dat_BD)
summary(am_BD) #p=0.767
# BDd fitted values
am_BD_fitted <- predict(am_BD, type = "response")

am_BD_dat <- am_beta_dat_BD %>% 
  select(am, BD_10, Latin) %>% 
  mutate(
    fitted = am_BD_fitted)

am_BD_p <- ggplot(am_BD_dat, aes(x = BD_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "BD", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_BD_p


#CBD
am_beta_dat_CBD <- am_beta_dat[!is.na(am_beta_dat$CBD_10), ]
am_CBD <- betareg(am ~ CBD_10, data = am_beta_dat_CBD)
summary(am_CBD) #p=0.0105
p_value <- coef(summary(am_CBD))$mean["CBD_10", "Pr(>|z|)"]
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
# CBDd fitted values
am_CBD_fitted <- predict(am_CBD, type = "response")

am_CBD_dat <- am_beta_dat_CBD %>% 
  select(am, CBD_10, Latin) %>% 
  mutate(
    fitted = am_CBD_fitted)

am_CBD_p <- ggplot(am_CBD_dat, aes(x = CBD_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "con.BD", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )+
  annotate("text", x = max(am_CBD_dat$CBD_10), y = max(am_CBD_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)
am_CBD_p


#invsimpson_div
am_beta_dat_invsimpson_div <- am_beta_dat[!is.na(am_beta_dat$invsimpson_div_10), ]
am_invsimpson_div <- betareg(am ~ invsimpson_div_10, data = am_beta_dat_invsimpson_div)
summary(am_invsimpson_div) #p=0.17
# invsimpson_divd fitted values
am_invsimpson_div_fitted <- predict(am_invsimpson_div, type = "response")

am_invsimpson_div_dat <- am_beta_dat_invsimpson_div %>% 
  select(am, invsimpson_div_10, Latin) %>% 
  mutate(
    fitted = am_invsimpson_div_fitted)

am_invsimpson_div_p <- ggplot(am_invsimpson_div_dat, aes(x = invsimpson_div_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "2D", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_invsimpson_div_p

#minpd
am_beta_dat_minpd <- am_beta_dat[!is.na(am_beta_dat$minpd_10), ]
am_minpd <- betareg(am ~ minpd_10, data = am_beta_dat_minpd)
summary(am_minpd) #p=0.17
p_value <- coef(summary(am_minpd))$mean["minpd_10", "Pr(>|z|)"]
# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

#dbh fitted values
am_minpd_fitted <- predict(am_minpd, type = "response")

am_minpd_dat <- am_beta_dat %>% 
  select(am, minpd_10, Latin) %>% 
  mutate(
    fitted = am_minpd_fitted)

am_minpd_p <- ggplot(am_minpd_dat, aes(x = minpd_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "min.Pd", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(am_minpd_dat$minpd_10), y = max(am_minpd_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_minpd_p)

#avepd
am_beta_dat_avepd <- am_beta_dat[!is.na(am_beta_dat$avepd_10), ]
am_avepd <- betareg(am ~ avepd_10, data = am_beta_dat_avepd)
summary(am_avepd) #p=0.4
# avepdd fitted values
am_avepd_fitted <- predict(am_avepd, type = "response")

am_avepd_dat <- am_beta_dat_avepd %>% 
  select(am, avepd_10, Latin) %>% 
  mutate(
    fitted = am_avepd_fitted)

am_avepd_p <- ggplot(am_avepd_dat, aes(x = avepd_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "ave.Pd", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_avepd_p

#totpd
am_beta_dat_totpd <- am_beta_dat[!is.na(am_beta_dat$totpd_10), ]
am_totpd <- betareg(am ~ totpd_10, data = am_beta_dat_totpd)
summary(am_totpd) #p=0.651
# totpdd fitted values
am_totpd_fitted <- predict(am_totpd, type = "response")

am_totpd_dat <- am_beta_dat_totpd %>% 
  select(am, totpd_10, Latin) %>% 
  mutate(
    fitted = am_totpd_fitted)

am_totpd_p <- ggplot(am_totpd_dat, aes(x = totpd_10, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "sum.Pd", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_totpd_p

#pcoa1
am_beta_dat_pcoa1 <- am_beta_dat[!is.na(am_beta_dat$pcoa1), ]
am_pcoa1 <- betareg(am ~ pcoa1, data = am_beta_dat_pcoa1)
summary(am_pcoa1) #p=0.00786
p_value <- coef(summary(am_pcoa1))$mean["pcoa1", "Pr(>|z|)"]
# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

# pcoa2d fitted values
am_pcoa1_fitted <- predict(am_pcoa1, type = "response")

am_pcoa1_dat <- am_beta_dat_pcoa1 %>% 
  select(am, pcoa1, Latin) %>% 
  mutate(
    fitted = am_pcoa1_fitted)

am_pcoa1_p <- ggplot(am_pcoa1_dat, aes(x = pcoa1, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "phylo_PCoA1", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(am_pcoa1_dat$pcoa1), y = max(am_pcoa1_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_pcoa1_p)

#pcoa2
am_beta_dat_pcoa2 <- am_beta_dat[!is.na(am_beta_dat$pcoa2), ]
am_pcoa2 <- betareg(am ~ pcoa2, data = am_beta_dat_pcoa2)
summary(am_pcoa2) #p=0.00786

p_value <- coef(summary(am_pcoa2))$mean["pcoa2", "Pr(>|z|)"]
# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", format(p_value, digits = 3))
}

# pcoa2d fitted values
am_pcoa2_fitted <- predict(am_pcoa2, type = "response")

am_pcoa2_dat <- am_beta_dat_pcoa2 %>% 
  select(am, pcoa2, Latin) %>% 
  mutate(
    fitted = am_pcoa2_fitted)

am_pcoa2_p <- ggplot(am_pcoa2_dat, aes(x = pcoa2, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "phylo_PCoA2", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(am_pcoa2_dat$pcoa2), y = max(am_pcoa2_dat$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_pcoa2_p)

#RD
am_beta_dat_RD <- am_beta_dat[!is.na(am_beta_dat$RDi), ]
am_RD <- betareg(am ~ RDi, data = am_beta_dat_RD)
summary(am_RD) #p=0.207
# RDd fitted values
am_RD_fitted <- predict(am_RD, type = "response")

am_RD_dat <- am_beta_dat_RD %>% 
  select(am, RDi, Latin) %>% 
  mutate(
    fitted = am_RD_fitted)

am_RD_p <- ggplot(am_RD_dat, aes(x = RDi, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "RD", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )
am_RD_p

library(cowplot)
lmplot <- plot_grid(am_dbh_p, am_pcoa1_p, am_pcoa2_p,
                            am_minpd_p,am_CBD_p, am_AD_p,
                            nrow = 2,
                            ncol = 3)
lmplot


ggsave(
  "pic/figure5.png",
  lmplot,
  width = 2000,
  height = 1100,
  units = "px",
  dpi = 300,
  bg = "#FFFFFF"
)




# 对比linear的结果
ggplot(am_dbh_dat, aes(x = DBH2, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  geom_smooth(method = "lm", se = F, linetype = "dashed") +
  labs(x = "DBH", y = "AMF colonization rates") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

ggsave(
  "test.png",
  am_dbh_p,
  dpi = 600
)