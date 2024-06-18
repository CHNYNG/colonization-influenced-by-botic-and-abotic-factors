#####画物种-侵染率图
library(dplyr)
library(ggplot2)
library(gridExtra)
library(eoffice)
# 使用ggplot2包绘制箱线图
ggplot(reg, aes(x = Latin, y = am)) +
  geom_boxplot() +
  labs(title = "Infection Rate by Species", x = "Species", y = "Infection Rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# 选择特定物种的数据
selected_speciesa <- c( "Engelhardtia_fenzlii", 
                      "Sinosideroxylon_wightianum", "Lithocarpus_uvariifolius", "Xanthophyllum_hainanense", 
                      "Symplocos_lancifolia")

filtered_dataa <- reg[reg$Latin %in% selected_speciesa,]
filtered_dataa$Latin <- gsub("_", " ", filtered_dataa$Latin)
# 定义物种名称的顺序
species_ordera <- c( "Xanthophyllum hainanense", "Engelhardtia fenzlii",
                     "Lithocarpus uvariifolius", "Sinosideroxylon wightianum",
                     "Symplocos lancifolia")#黄叶树 少叶黄杞 革叶铁榄 水杨梅 南酸枣 光叶山矾


# 将物种名称转换为因子变量，并按照指定的顺序排序
filtered_dataa$Latin <- factor(filtered_dataa$Latin, levels = species_ordera)

# 使用ggplot2包绘制筛选后的箱线图，旋转横坐标上的文字
plot1 <- ggplot(filtered_dataa, aes(x = am , y = Latin)) +
  geom_boxplot() +
  labs(title = "Arbuscular Mycorrhiza", y = "Species",  x = "Mycorrhizal colonization rate") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle =0, vjust = 00, hjust= 0, face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(hjust = 1),
        plot.title = element_blank())+
       # plot.title = element_text(hjust = -0.5)) +
  coord_cartesian(xlim = c(0.6, 0.9)) # 设置纵坐标范围
print(plot1)
topptx(plot1,filename = "mycorrhizal colonization rate.pptx")
####em
ggplot(reg, aes(x = Latin, y = em)) +
  geom_boxplot() +
  labs(title = "Infection Rate by Species", x = "Species", y = "Infection Rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

selected_speciese <- c("Cyclobalanopsis_chungii", "Cyclobalanopsis_fleuryi", "Cyclobalanopsis_hui", 
                      "Cyclobalanopsis_bambusaefolia", "Castanopsis_carlesii", "Lithocarpus_litseifolius", 
                      "Castanopsis_fordii")#福建青冈 饭甑青冈 雷公青冈 竹叶青冈 米锥 甜茶椆 毛椎
# 定义物种名称的顺序
species_ordere <- c( "Cyclobalanopsis_chungii", "Cyclobalanopsis_hui",  
                    "Cyclobalanopsis_bambusaefolia","Cyclobalanopsis_fleuryi", "Castanopsis_carlesii", "Lithocarpus_litseifolius", 
                    "Castanopsis_fordii")


filtered_datae <- reg[reg$Latin %in% selected_speciese,]
# 将物种名称转换为因子变量，并按照指定的顺序排序
filtered_datae$Latin <- factor(filtered_datae$Latin, levels = species_ordere)

# 使用ggplot2包绘制筛选后的箱线图，旋转横坐标上的文字
plot2 <- ggplot(filtered_datae, aes(y = Latin, x = em)) +
  geom_boxplot() +
  labs(title = "Ecomycorrhiza", y = "Species",  x = "Mycorrhizal colonization rate") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle =0, vjust = 0.5, hjust=0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_text(hjust = -0.8),
        plot.title = element_text(hjust = -0.5)) +
  coord_cartesian(xlim = c(0.6, 0.9)) # 设置纵坐标范围
print(plot2)

#拼接
grid.arrange(plot1, plot2, nrow = 2)
