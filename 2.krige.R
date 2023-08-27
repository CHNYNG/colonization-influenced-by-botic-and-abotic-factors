##### 
#加载前步骤
source("1.colonization_data_caculate.R")

##### 
#土壤差值
soil <- read.csv("data/HSD_Soil_Nutrition_Liu.csv",header = T,fileEncoding = "GBK")
####把数据和坐标合并一下
soil_to_liu <- read.csv("data/Ru Wang/HSD_Soil_To Liu.csv",header = T,fileEncoding = "GBK")
soil <- left_join(soil,soil_to_liu,by = "SampleID")
soil <- na.omit(soil)
##整理已知土壤数据
library(tidyverse)
soil <- soil %>% 
  rename(
    sampleid = SampleID,
    soc = SOC..mg.g.,
    tn = TN..mg.g.,
    tp = TP..mg.g.,
    cn = C.N,
    cp = C.P,
    np = N.P,
    ap = AP..mg.kg.,
    ph = pH,
    gx = GX,
    gy = GY,
    tagnew = TagNew
  )
###变量包括 "soc","tn","tp","cn","cp","np","ap","ph"
gx_column <- as.numeric(which(names(soil) == "gx"))
gy_column <- as.numeric(which(names(soil) == "gy"))

####整理需要预测的点的数据集prd_loc
library(geoR)
# 使用complete.cases检查root_qrl中每一行是否包含NA值
# 然后用subset函数筛选掉包含NA值的行
prd_loc <- subset(root_qrl,
                  complete.cases(GX, GY, TagNew),
                  select = c(GX, GY, TagNew))

prd_loc$GX <- as.numeric(prd_loc$GX)


#####
#其余的一起做
#####
#以soc为例，创建soc的已知数据的geodata
soc_g <- as.geodata(soil, 
                    coords.col = gx_column:gy_column,  
                    data.col = as.numeric(which(names(soil) == "soc")))

##plot(soc_g)
##plot(soc_g, trend = "1st")
##plot(soc_g, trend = "2nd") # it seems that we should consider 2nd trend

#师兄说，这一行不用知道为啥，做就行
soc_g_maxdis <- summary(soc_g)$distances.summary[2] / 2 # we usually do this

#这里的trend要选择一个上面认为表现最好的
soc_variog <- variog(
  soc_g,
  trend = "2nd",
  max.dist = soc_g_maxdis
)


#####
#tn
#####
tn_g <- as.geodata(soil, 
                    coords.col = gx_column:gy_column,  
                    data.col = as.numeric(which(names(soil) == "tn")))

#plot(tn_g)
#plot(tn_g, trend = "1st")
#plot(tn_g, trend = "2nd")

#师兄说，这一行不用知道为啥，做就行
tn_g_maxdis <- summary(tn_g)$distances.summary[2] / 2 # we usually do this

#这里的trend要选择一个上面认为表现最好的
tn_variog <- variog(
  tn_g,
  trend = "2nd",
  max.dist = tn_g_maxdis
)
#####
#tp
tp_g <- as.geodata(soil, 
                   coords.col = gx_column:gy_column,  
                   data.col = as.numeric(which(names(soil) == "tp")))

#plot(tp_g)
#plot(tp_g, trend = "1st")
#plot(tp_g, trend = "2nd")


tp_g_maxdis <- summary(tp_g)$distances.summary[2] / 2 # we usually do this

#这里的trend要选择一个上面认为表现最好的
tp_variog <- variog(
  tp_g,
  trend = "2nd",
  max.dist = tp_g_maxdis
)
#####
#ap
ap_g <- as.geodata(soil, 
                   coords.col = gx_column:gy_column,  
                   data.col = as.numeric(which(names(soil) == "ap")))

#plot(ap_g)
#plot(ap_g, trend = "1st")
#plot(ap_g, trend = "2nd")


ap_g_maxdis <- summary(ap_g)$distances.summary[2] / 2 # we usually do this

#这里的trend要选择一个上面认为表现最好的
ap_variog <- variog(
  ap_g,
  trend = "2nd",
  max.dist = ap_g_maxdis
)

#####
#ph
ph_g <- as.geodata(soil, 
                   coords.col = gx_column:gy_column,  
                   data.col = as.numeric(which(names(soil) == "ph")))

#plot(ph_g)
#plot(ph_g, trend = "1st")
#plot(ph_g, trend = "2nd")


ph_g_maxdis <- summary(ph_g)$distances.summary[2] / 2 # we usually do this

#这里的trend要选择一个上面认为表现最好的
ph_variog <- variog(
  ph_g,
  trend = "cte",
  max.dist = ph_g_maxdis
)
#####
#函数
#####
###师兄把模型筛选的过程集合成了一个函数，先跑这个函数
auto_kriging <- function(variogram_mod, geo_dat) {
  # fitting covariance models on a variogram
  exp_mod <-
    variofit(variogram_mod,
             cov.model = "exp",
             fix.nugget = T)
  
  exp_mod_nofix <-
    variofit(variogram_mod,
             cov.model = "exp",
             fix.nugget = F)
  
  sph_mod <-
    variofit(variogram_mod,
             cov.model = "sph",
             fix.nugget = T)
  
  sph_mod_nofix <-
    variofit(variogram_mod,
             cov.model = "sph",
             fix.nugget = F)
  
  gau_mod <- 
    variofit(variogram_mod,
             cov.model = "gaussian",
             fix.nugget = T)
  
  gau_mod_nofix <- 
    variofit(variogram_mod,
             cov.model = "gaussian",
             fix.nugget = F)
  
  variofit_list <- list(
    exp_mod,
    exp_mod_nofix,
    sph_mod,
    sph_mod_nofix,
    gau_mod,
    gau_mod_nofix
  )
  
  # model selection: sph_soc_nofix ----
  exp_mod_ssq <- summary(exp_mod)$sum.of.squares
  exp_mod_nofix_ssq <- summary(exp_mod_nofix)$sum.of.squares
  sph_mod_ssq <- summary(sph_mod)$sum.of.squares
  sph_mod_nofix_ssq <- summary(sph_mod_nofix)$sum.of.squares
  gau_mod_ssq <- summary(gau_mod)$sum.of.squares
  gau_mod_nofix_ssq <- summary(gau_mod_nofix)$sum.of.squares
  
  mod_num <- which.min(list(
    exp_mod_ssq,
    exp_mod_nofix_ssq,
    sph_mod_ssq,
    sph_mod_nofix_ssq,
    gau_mod_ssq,
    gau_mod_nofix_ssq
  ))
  
  # prediction
  soil_prd <- krige.conv(
    geo_dat,
    loc = prd_loc,
    krige = krige.control(obj.model = variofit_list[[mod_num]])
  )
  
  # data frame
  prd_dat <- data.frame(
    gx = prd_loc$GX,
    gy = prd_loc$GY,
    prediction = soil_prd$predict
  )
  
  #
  return(prd_dat)
}

######
#然后把每个变量过一遍
#####
soc_prd_dat <- auto_kriging(soc_variog, soc_g)
tn_prd_dat <- auto_kriging(tn_variog, tn_g)
tp_prd_dat <- auto_kriging(tp_variog, tp_g)
ap_prd_dat <- auto_kriging(ap_variog, ap_g)
ph_prd_dat <- auto_kriging(ph_variog, ph_g)

##之后把它们合在一起
soil_pred <- data.frame(
  GX = prd_loc$GX,
  GY = prd_loc$GY,
  TagNew = prd_loc$TagNew,
  soc = soc_prd_dat$prediction,
  tn = tn_prd_dat$prediction,
  tp = tp_prd_dat$prediction,
  ap = ap_prd_dat$prediction,
  ph = ph_prd_dat$prediction
)
#####
#删除一些过程变量
#####
rm(ap_g,ap_prd_dat,ap_variog,ph_g,ph_prd_dat,ph_variog,
   soc_g,soc_prd_dat,soc_variog,tn_g,tn_prd_dat,tn_variog,
   tp_g,tp_prd_dat,tp_variog,ap_g_maxdis,gx_column,gy_column,
   ph_g_maxdis,soc_g_maxdis,tn_g_maxdis,tp_g_maxdis)
