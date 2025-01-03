##### summary the base data of mycorrhiza colonization #####
summary(reg$qr_AM)
mean(reg$qr_AM, na.rm = TRUE)
sd(reg$qr_AM, na.rm = TRUE)
#Rhododendron simsii
mean(reg$qr_AM[reg$Latin == "Rhododendron simsii"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Rhododendron simsii"], na.rm = TRUE)
#Castanopsis fissa
mean(reg$qr_AM[reg$Latin == "Castanopsis fissa"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Castanopsis fissa"], na.rm = TRUE)
#Ixonanthes chinensis
mean(reg$qr_AM[reg$Latin == "Ixonanthes chinensis"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Ixonanthes chinensis"], na.rm = TRUE)
#Craibiodendron stellatum
mean(reg$qr_AM[reg$Latin == "Craibiodendron stellatum"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Craibiodendron stellatum"], na.rm = TRUE)
#Exbucklandia tonkinensis
mean(reg$qr_AM[reg$Latin == "Exbucklandia tonkinensis"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Exbucklandia tonkinensis"], na.rm = TRUE)
#Garcinia multiflora
mean(reg$qr_AM[reg$Latin == "Garcinia multiflora"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Garcinia multiflora"], na.rm = TRUE)
#Ardisia elegans 
mean(reg$qr_AM[reg$Latin == "Ardisia elegans"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Ardisia elegans"], na.rm = TRUE)
#Cryptocarya concinna
mean(reg$qr_AM[reg$Latin == "Cryptocarya concinna"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Cryptocarya concinna"], na.rm = TRUE)
#Vitex quinata
mean(reg$qr_AM[reg$Latin == "Vitex quinata"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Vitex quinata"], na.rm = TRUE)
#Craibiodendron stellatum
mean(reg$qr_AM[reg$Latin == "Lithocarpus lohangwu"], na.rm = TRUE)
sd(reg$qr_AM[reg$Latin == "Lithocarpus lohangwu"], na.rm = TRUE)

##### beta-regression models #####
library(glmmTMB)
library(piecewiseSEM)
summary(am_beta_mod_optimal)
contributions

##### discussion #####
####### paragraph 1 ########
summary(reg)
length(unique(reg$Latin))
length(am_beta_dat$Latin)

####### section 1 ######
