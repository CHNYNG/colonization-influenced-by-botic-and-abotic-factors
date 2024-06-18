# for the results
#paragraph 1
library(dplyr)
Species <- reg %>%
  group_by(Latin) %>%
  summarize(mean_am = mean(am, na.rm = TRUE),
            sd_am = sd(am, na.rm = TRUE)^2) %>%
  arrange(desc(mean_am))
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
