# global variable
prd_loc <- data.frame(
  x = as.numeric(prd.loc$GX),
  y = prd.loc$GY
)

# defining function to kriging
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
    gx = prd_loc$x,
    gy = prd_loc$y,
    prediction = soil_prd$predict
  )
  
  #
  return(prd_dat)
}

# generating prediction data
soc_prd_dat <- auto_kriging(soc_variog, soc_g)
tn_prd_dat <- auto_kriging(tn_variog, tn_g)
tp_prd_dat <- auto_kriging(tp_variog, tp_g)
ap_prd_dat <- auto_kriging(ap_variog, ap_g)
ph_prd_dat <- auto_kriging(ph_variog, ph_g)

soil_cy_pred <- data.frame(
  gx = prd_loc$x,
  gy = prd_loc$y,
  soc = soc_prd_dat$prediction,
  tn = tn_prd_dat$prediction,
  tp = tp_prd_dat$prediction,
  ap = ap_prd_dat$prediction,
  ph = ph_prd_dat$prediction
)



################################################################################
# not run
ph_g <- as.geodata(soil_cy, 
                   coords.col = gx_column:gy_column, 
                   data.col = which(names(soil_cy) == "ph"))
plot(ph_g) # better
plot(ph_g, trend = "1st")
plot(ph_g, trend = "2nd") 

ph_g_maxdis <- summary(ph_g)$distances.summary[2] / 2 # we usually do this

ph_variog <- variog(
  ph_g,
  # trend = "1st",
  max.dist = ph_g_maxdis
)

