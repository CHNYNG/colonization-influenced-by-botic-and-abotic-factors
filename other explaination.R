###试一下会不会有哪个科、属的情况不同
# 加载 dplyr 包
library(dplyr)

# 选择特定的 Family 作为新的数据框
Fabaceae <- am_beta_dat %>% filter(Family == "Fabaceae")
library(tidyr)
library(ggplot2)
library(betareg)


####dbh
##### RAD #####
Fabaceae_AD <- Fabaceae[!is.na(Fabaceae$AD), ]
Fabaceae_AD <- Fabaceae_AD %>%
  filter(AD!=max(AD))
am_AD <- betareg(am ~ AD, data = Fabaceae_AD)
summary(am_AD) #p=0.0294

p_value <- coef(summary(am_AD))$mean["AD", "Pr(>|z|)"]
AD_coef <- coef(summary(am_AD))$mean["AD", "Estimate"]
# 设置 p 值标签，保留三位小数
if (p_value < 0.001) {
  p_label <- "italic(p) < 0.001"
} else {
  p_label <- paste0("italic(p) == ", sprintf("%.3f", p_value))
}
coef_label <- paste0("italic(beta) == ", sprintf("%.3f", AD_coef))
# Add fitted values
am_AD_fitted <- predict(am_AD, type = "response")

am_AD_dat <- Fabaceae_AD %>% 
  dplyr::select(am, AD, Latin) %>% 
  dplyr::mutate(
    fitted = am_AD_fitted)
line_type <- ifelse(p_value > 0.01, "dashed", "solid")

am_AD_p <- ggplot(am_AD_dat, aes(x = AD, y = am)) +
  geom_point(
    size = 2,
    color = "darkgrey",
    alpha = 0.7
  ) +
  geom_line(aes(y = fitted),
            color = "#004529",
            linewidth = 1.2,
            linetype = line_type) +
  labs(x = "RAD", y =NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )+
  annotate("text", x = max(Fabaceae_AD$AD), y = max(Fabaceae_AD$am), 
           label =  paste(coef_label, p_label, sep = '~~~'), 
           parse = TRUE, hjust = 1, vjust = 1, size = 5)
am_AD_p

#### 直接出每个科的β和p ####
library(dplyr)
library(betareg)

# 创建一个空的数据框来存储结果
results <- data.frame(
  Family = character(),
  SpeciesCount = integer(),
  Beta_Coefficient = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# 循环遍历每个 Family
for(family in unique(am_beta_dat_AD$Family)) {
  # 筛选出当前 Family 的数据
  family_data <- am_beta_dat_AD %>% filter(Family == family)
  
  # 计算物种数
  species_count <- nrow(family_data)
  
  # 检查是否有足够的数据点进行回归
  if(nrow(family_data) > 5) {
    # 进行 beta 回归
    model <- betareg(am ~ AD, data = family_data)
    
    # 提取回归系数和 p 值
    beta_coef <- coef(summary(model))$mean["AD", "Estimate"]
    p_value <- coef(summary(model))$mean["AD", "Pr(>|z|)"]
  } else {
    # 如果数据不足，则填 NA
    beta_coef <- NA
    p_value <- NA
  }
  
  # 将结果添加到数据框中
  results <- rbind(results, data.frame(
    Family = family,
    SpeciesCount = species_count,
    Beta_Coefficient = beta_coef,
    P_Value = p_value
  ))
}

# 查看结果
print(results)

