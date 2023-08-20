# required packages ----
library(tidyverse)
library(geoR)

# data preparation ----
soil_cy <- Soil %>% 
  rename(
    soc = SOC..mg.g.,
    tn = TN..mg.g.,
    tp = TP..mg.g.,
    cn = C.N,
    cp = C.P,
    np = N.P,
    ap = AP..mg.kg.,
    ph = pH,
    gx = GX,
    gy = GY
  )

gx_column <- which(names(soil) == "gx")
gy_column <- which(names(soil) == "gy")

prd_loc <- data.frame(
  x = as.numeric(prd.loc$GX),
  y = prd.loc$GY
)

# krging ----
# taken soc as an example
soc_g <- as.geodata(soil_cy, 
                    coords.col = gx_column:gy_column, 
                    data.col = which(names(soil_cy) == "soc"))
plot(soc_g)
plot(soc_g, trend = "1st")
plot(soc_g, trend = "2nd") # it seems that we should consider 2nd trend

soc_g_maxdis <- summary(soc_g)$distances.summary[2] / 2 # we usually do this

soc_variog <- variog(
  soc_g,
  trend = "2nd",
  max.dist = soc_g_maxdis
)

# plot(soc_variog)

# fitting covariance models on a variogram
exp_soc <-
  variofit(soc_variog,
           cov.model = "exp",
           fix.nugget = T)

exp_soc_nofix <-
  variofit(soc_variog,
           cov.model = "exp",
           fix.nugget = F)

sph_soc <-
  variofit(soc_variog,
           cov.model = "sph",
           fix.nugget = T)

sph_soc_nofix <-
  variofit(soc_variog,
           cov.model = "sph",
           fix.nugget = F)

gau_soc <- 
  variofit(soc_variog,
           cov.model = "gaussian",
           fix.nugget = T)

gau_soc_nofix <- 
  variofit(soc_variog,
           cov.model = "gaussian",
           fix.nugget = F)

# you do not need to run the codes below
plot(soc_variog, pch = 16)
lines(exp_soc,
      col = "green",
      lwd = 4,
      lty = 1)

plot(soc_variog, pch = 16)
lines(exp_soc_nofix,
      col = "red",
      lwd = 4,
      lty = 2)

plot(soc_variog, pch = 16)
lines(sph_soc,
      col = "blue",
      lwd = 4,
      lty = 3) 

plot(soc_variog, pch = 16)
lines(sph_soc_nofix,
      col = "orange",
      lwd = 4,
      lty = 4)

plot(soc_variog, pch = 16)
lines(gau_soc,
      col = "lightblue",
      lwd = 4,
      lty = 3) 

plot(soc_variog, pch = 16)
lines(gau_soc_nofix,
      col = "purple",
      lwd = 4,
      lty = 4)

# model selection: sph_soc_nofix ----
exp_soc_ssq <- summary(exp_soc)$sum.of.squares
exp_soc_nofix_ssq <- summary(exp_soc_nofix)$sum.of.squares
sph_soc_ssq <- summary(sph_soc)$sum.of.squares
sph_soc_nofix_ssq <- summary(sph_soc_nofix)$sum.of.squares
gau_soc_ssq <- summary(gau_soc)$sum.of.squares
gau_soc_nofix_ssq <- summary(gau_soc_nofix)$sum.of.squares

which.min(list(
  exp_soc_ssq,
  exp_soc_nofix_ssq,
  sph_soc_ssq,
  sph_soc_nofix_ssq, # seems that we should choose this model
  gau_soc_ssq,
  gau_soc_nofix_ssq
))

# final kriging, we just use simple kriging
soc_prd <- krige.conv(
  soc_g,
  loc = prd_loc,
  krige = krige.control(obj.model = sph_soc_nofix)
)

# create data frame
soc_prd_dat <- data.frame(
  gx = prd_loc$x,
  gy = prd_loc$y,
  soc = soc_prd$predict
)

# performing cross validation of results, you don't need always run this
soc_xv <- xvalid(soc_g, model = sph_soc_nofix)

soc_xv_df <- data.frame(soc_xv$data, soc_xv$predicted)

ggplot(data = soc_xv_df, aes(x = soc_xv.data, y = soc_xv.predicted)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  xlab("Observed soc") +
  ylab("Predicted soc") +
  labs(title = "Linear regression model for predicted elevation")
# seems the model performs badly, but we can still use the results

