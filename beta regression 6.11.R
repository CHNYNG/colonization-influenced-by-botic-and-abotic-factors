# new beta regression by weitao 
# required packages ----
library(tidyverse)
library(glmmTMB)
library(piecewiseSEM)
library(performance)
library(MuMIn)

load("D:/Users/65205/Documents/WeChat Files/wxid_ulzapgmya6hv22/FileStorage/File/2024-06/cy 6.7.RData")
# data preparation
# origninal data set is reg_scaled
# we can exclude ecm species and re-scaled data
# try how to get the pcoa
#first obtain the distance matrix of the species phylo
library(ape)
library(phangorn)
library(vegan)
species_to_keep <- am_phylo_pcoa_axes$Latin
filtered_tree <- keep.tip(hsd_phytr, species_to_keep)
dist_matrix <- cophenetic(filtered_tree)
dist_matrix <- sqrt(dist_matrix)
decostand(dist_matrix, "range")

#second, play the pcoa
library(vegan)
#pcoa_result <- pcoa(dist_matrix)
pcoa_result <- cmdscale(dist_matrix) 


am_phylo_pcoa_axes_try <- as.data.frame(as.data.frame(pcoa_result))
am_phylo_pcoa_axes_try <- decostand(am_phylo_pcoa_axes_try, "standardize", MARGIN = 2)
am_phylo_pcoa_axes_try <- as.data.frame(am_phylo_pcoa_axes_try)
summary(am_phylo_pcoa_axes_try)
summary(am_phylo_pcoa_axes)
#pcoa_result[["Value"]]
#把pcoa的结果中心化一下？

am_beta_dat <- reg_scaled %>%
  filter(!is.na(am)) %>%
  select(
    # species identity
    Latin,
    Genus,
    Family,
    # resource diversity
    RDi, RRi, REi,
    # size
    DBH2,
    # root traits
    AD,
    SRL,
    SRA,
    # topography
    elevation,
    slope,
    aspect,
    convexity,
    # phylogenetic
    minpd_10,
    avepd_10,
    totpd_10,
    # neighbor effect
    CBD_10,
    BD_10,
    invsimpson_div_10,
  ) %>%
  left_join(am_phylo_pcoa_axes, by = "Latin") %>%
#  mutate_if(is.numeric, scale) %>%
  mutate(
    rd_is = scale(RDi * invsimpson_div_10),
    rr_is = scale(RRi * invsimpson_div_10),
    re_is = scale(REi * invsimpson_div_10),
    am_beta_dat %>% filter(!is.na(am)) %>% select(am)
  )

# I only calculate the pcoa axes of 161 species, and put them in: am_phylo_pcoa_axes

# beta regression and model selection
# based on your manuscript, we can set up the full model like this:
library(betareg)
options(na.action = "na.omit")
am_beta_mod_ful <- betareg(
  am ~ minpd_10 + avepd_10 + totpd_10 + AD + SRL + SRA + DBH2 + CBD_10 + BD_10 + 
    elevation + invsimpson_div_10 + RDi + pcoa1 + pcoa2 + rd_is,
  data = am_beta_dat
)
summary(am_beta_mod_ful)

# R can automatically choose optimal model for us based on AIC
options(na.action = "na.fail")
am_beta_mod_sub <- dredge(am_beta_mod_ful)


subset(am_beta_mod_sub, delta < 2) # optimal model and those with delta AICc < 2
am_beta_mod_minaic <- summary(get.models(am_beta_mod_sub, 1)[[1]]) # optimal model
am_beta_mod_ave <- summary(model.avg(am_beta_mod_sub, subset = delta < 2)) # averaging method
# am_beta_mod_ave is another choice if you like, see more in Burnham, K. P. and Anderson, D. R. 2002

# if you want to make everything simple, use am_beta_mod_minaic
# and based on am_beta_mod_minaic:
am_beta_mod_optimal <- betareg(
  am ~ minpd_10 + avepd_10 + DBH2 + CBD_10 + invsimpson_div_10 + RDi + pcoa1 + pcoa2 + rd_is,
  data = am_beta_dat
)
summary(am_beta_mod_optimal)

if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)

# you will find it's wired that in am_beta_mod_minaic, the RDi is excluded
# but in am_beta_mod_optimal, I included it!
# The reason is that if we keep the interaction term in the model,
# it's better to keep both of the original variables that create the interaction term.


# so, basically am_beta_mod_optimal is the final model you should report in your article
# it contains some variables that are not statistically significant,
# it doesn't matter
# constructing your psem with all the variables in am_beta_mod_optimal

#画图
# 安装并加载所需的包
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")
library(ggplot2)
library(gridExtra)

# 提取估计值和标准误差
summary_obj <- summary(am_beta_mod_optimal)

estimates <- as.vector(summary_obj$coefficients$mean[,"Estimate"])
std_errors <- as.vector(summary_obj$coefficients$mean[, "Std. Error"])

# 去掉(Intercept)的行
estimates <- estimates[-1]
std_errors <- std_errors[-1]

variables <- c("min.Pd", "ave.Pd", "DBH", "con.BD", "IS", "RD", "pcoa1", "pcoa2", "RD*IS")

# 计算估计值的绝对值
abs_estimates <- abs(estimates)

# 计算总的估计值绝对值
total_abs_estimate <- sum(abs_estimates)

# 计算每个变量的贡献值
contributions <- abs_estimates / total_abs_estimate

# 创建数据框
data <- data.frame(
  Variable = factor(variables, levels = c("RD*IS", "RD", "DBH", "pcoa2", "pcoa1", "ave.Pd", "min.Pd", "IS", "con.BD")),
  Estimate = estimates,
  StdError = std_errors,
  Contribution = contributions,
  CumulativeContribution = cumsum(contributions)
)
colors <- c("con.BD" = "#004529", "IS" = "#004529", "min.Pd" = "#006837", "ave.Pd" = "#006837", 
            "pcoa1" = "#41ab5d", "pcoa2" = "#41ab5d", "DBH" = "#addd8e", "RD" = "#feb24c", "RD*IS" = "#ffeda0")

y_cols <- c("#ffeda0" , "#feb24c", "#addd8e","#41ab5d", "#41ab5d", "#006837", "#006837","#00452c", "#004529")
# 创建左侧的估计值和标准误差图
# 创建左侧的估计值和标准误差图
p1 <- ggplot(data, aes(x = Estimate, y = Variable, color = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Estimate - StdError, xmax = Estimate + StdError), height = 0.1, size = 2) +
  geom_vline(xintercept =   0, linetype = "dashed", color = "#bdbdbd", size = 1) + # 添加竖线
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(legend.position = "none",# 去掉图例
        axis.text.y = element_text(colour = y_cols, size = 20),
        axis.text.x = element_text(size = 20, colour = "#252525"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
       # panel.grid.major = element_blank(),# 去掉主网格线
        panel.grid.minor = element_blank()   # 去掉次网格线
  )
p1
# 创建右侧的累计贡献值柱状图
# 创建右侧的累计贡献值柱状图，按反序排列变量
p2 <- ggplot(data, aes(x = 1, y = Contribution, fill = factor(Variable, levels = rev(levels(data$Variable))))) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = colors) +
#  labs(x = "Contribution", y = NULL) +
  scale_y_continuous(labels = scales::percent, position = "right") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),#element_text(size = 20, colour = "#252525"),
    axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20, colour = "#252525"),
        legend.position = "none",
    panel.grid.major = element_blank(),# 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次网格线
    )
p2

# 将两个图放在一起
#pdf(file = "output.pdf", width = 6, height = 4)
library(cowplot)
plot_grid(p1, p2, rel_widths = c(5, 1))
#grid.arrange(p1, p2, ncol = 2)
#导出到ppt里改吧。。。
library(eoffice)
library(ggplotify)
topptx(p, filename = "~/pic/eoffice.pptx")
