library(betareg)
library(ggplot2)
library(dplyr)
library(tidyr)
#############selet tn tp and tn/tp and colonization rate############ 

df <- reg_scaled %>%
  select(tn, tp, am, soc, ap, ph, moisture, tndtp, Latin) %>%       # 选择需要的列
  drop_na() 

# 查看处理后的数据框
head(df)

##### tn #########
am_tn <- betareg(am ~ tn, data = df)
summary(am_tn) #显著，留下
p_value <- coef(summary(am_tn))$mean["tn", "Pr(>|z|)"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

#tn fitted values
am_tn_fitted <- predict(am_tn, type = "response", na.action = na.exclude)

am_tn_dat <- df %>% 
  select(am, tn, Latin) %>% 
  mutate(
    fitted = am_tn_fitted)

am_tn_p <- ggplot(am_tn_dat, aes(x = tn, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "TN", y = "AMFR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tn), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tn_p)

#### tp ####
am_tp <- betareg(am ~ tp, data = df)
summary(am_tp) #显著，留下
p_value <- coef(summary(am_tp))$mean["tp", "Pr(>|z|)"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

#tp fitted values
am_tp_fitted <- predict(am_tp, type = "response", na.action = na.exclude)

am_tp_dat <- df %>% 
  select(am, tp, Latin) %>% 
  mutate(
    fitted = am_tp_fitted)

am_tp_p <- ggplot(am_tp_dat, aes(x = tp, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "TP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tp), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tp_p)

### tn/tp ####

###### fillter ####

df <- df %>%
  filter(tndtp != max(tndtp, na.rm = TRUE))
am_tndtp <- betareg(am ~ tndtp, data = df)
summary(am_tndtp) #显著，留下
p_value <- coef(summary(am_tndtp))$mean["tndtp", "Pr(>|z|)"]

###### 设置 p 值标签 ######
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

###### tndtp fitted values #####
am_tndtp_fitted <- predict(am_tndtp, type = "response", na.action = na.exclude)

am_tndtp_dat <- df %>% 
  select(am, tndtp, Latin) %>% 
  mutate(
    fitted = am_tndtp_fitted)

am_tndtp_p <- ggplot(am_tndtp_dat, aes(x = tndtp, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "TN/TP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tndtp), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tndtp_p)

### soc ####

###### fillter ####

am_soc <- betareg(am ~ soc, data = df)
summary(am_soc) #显著，留下
p_value <- coef(summary(am_soc))$mean["soc", "Pr(>|z|)"]

###### 设置 p 值标签 ######
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

###### soc fitted values #####
am_soc_fitted <- predict(am_soc, type = "response", na.action = na.exclude)

am_soc_dat <- df %>% 
  select(am, soc, Latin) %>% 
  mutate(
    fitted = am_soc_fitted)

am_soc_p <- ggplot(am_soc_dat, aes(x = soc, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "SOC", y = "AMFR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$soc), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_soc_p)

### ap ####

###### fillter ####

am_ap <- betareg(am ~ ap, data = df)
summary(am_ap) #显著，留下
p_value <- coef(summary(am_ap))$mean["ap", "Pr(>|z|)"]

###### 设置 p 值标签 ######
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

###### ap fitted values #####
am_ap_fitted <- predict(am_ap, type = "response", na.action = na.exclude)

am_ap_dat <- df %>% 
  select(am, ap, Latin) %>% 
  mutate(
    fitted = am_ap_fitted)

am_ap_p <- ggplot(am_ap_dat, aes(x = ap, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "AP", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$ap), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_ap_p)

### ph ####

###### fillter ####

am_ph <- betareg(am ~ ph, data = df)
summary(am_ph) #显著，留下
p_value <- coef(summary(am_ph))$mean["ph", "Pr(>|z|)"]

###### 设置 p 值标签 ######
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

###### ph fitted values #####
am_ph_fitted <- predict(am_ph, type = "response", na.action = na.exclude)

am_ph_dat <- df %>% 
  select(am, ph, Latin) %>% 
  mutate(
    fitted = am_ph_fitted)

am_ph_p <- ggplot(am_ph_dat, aes(x = ph, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "pH", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$ph), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_ph_p)

### moisture ####

###### fillter ####

am_moisture <- betareg(am ~ moisture, data = df)
summary(am_moisture) #显著，留下
p_value <- coef(summary(am_moisture))$mean["moisture", "Pr(>|z|)"]

###### 设置 p 值标签 ######
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

###### moisture fitted values #####
am_moisture_fitted <- predict(am_moisture, type = "response", na.action = na.exclude)

am_moisture_dat <- df %>% 
  select(am, moisture, Latin) %>% 
  mutate(
    fitted = am_moisture_fitted)

am_moisture_p <- ggplot(am_moisture_dat, aes(x = moisture, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    almoisturea = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2) +
  labs(x = "moisture", y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$moisture), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_moisture_p)
#### output the pictures ####

library(cowplot)
lmplots2 <- plot_grid(am_soc_p, am_ap_p, am_ph_p, am_moisture_p,
                    am_tn_p,am_tp_p, am_tndtp_p,
                    nrow = 2,
                    ncol = 4)

lmplots2

ggsave(
    "pic/figure S2.png",
    lmplots2,
    width = 2500,
    height = 1100,
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

