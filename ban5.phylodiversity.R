#首先说明我的物种
#####
head(specieslist)

#####
####把这些物种按照数据库建树
library("V.PhyloMaker")
library(ape)
#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggtree")
library(ggtree)
library(ggplot2)
library(stringr)
specieslist <- na.omit(specieslist)
sp_tree_v.p <- phylo.maker(sp.lis=specieslist)
sp_tr <- sp_tree_v.p$scenario.3
pdf(file = "pic/tree_of_vphylo.pdf", width = 8.27, height = 50)
plot(sp_tr)
dev.off()

#####
#把韦韬师兄给我的谱系树拉出来
sp_tree_hsd <- read.tree("data/hsd_phylo_tree.nwk")
pdf(file = "pic/tree_of_hsd.pdf", width = 8.27, height = 40)
plot(sp_tree_hsd)
dev.off()

#####
#来！算PD！
#算PD前要把每一个小群落提出来，先去4.算出小群落


