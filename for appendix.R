library(betareg)
library(ggplot2)
library(dplyr)
library(tidyr)
#############selet tn tp and tn/tp and colonization rate############ 

df <- reg %>%
  select(tn, tp, am, Latin) %>%       # 选择需要的列
  mutate(tndtp = tn / tp) %>%     # 创建新变量 tndtp
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

##### tp #######
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
  labs(x = "TP", y = "AMFR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tp), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tp_p)

##### tn/tp #####
am_tndtp <- betareg(am ~ tndtp, data = df)
summary(am_tndtp) #显著，留下
p_value <- coef(summary(am_tndtp))$mean["tndtp", "Pr(>|z|)"]

# 设置 p 值标签
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}

#tndtp fitted values
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
  labs(x = "TN/TP", y = "AMFR") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  annotate("text", x = max(df$tndtp), y = max(df$am), 
           label = p_label, parse = TRUE, hjust = 1, vjust = 1, size = 5)

print(am_tndtp_p)

##### output the pictures #####

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
