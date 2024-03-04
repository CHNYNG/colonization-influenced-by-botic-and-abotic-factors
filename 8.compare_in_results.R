####关于模型的结果探究

###不对，这里它们参数不同，不能比较
#这样的方法应该用于比较glm_10内部的不同
#比较glm10和20之间的差异
lr_test_10_20 <- anova(glm_10, glm_20, test = "Chisq")
print(lr_test_10_20)
#比较glm10和50之间的差异
lr_test_10_50 <- anova(glm_10, glm_50, test = "Chisq")
print(lr_test_10_50)
#比较50和20
lr_test_20_50 <- anova(glm_20, glm_50, test = "Chisq")
print(lr_test_20_50)
lr_test_50_20 <- anova(glm_50, glm_20, test = "Chisq")
print(lr_test_50_20)

####再试试
#将结果的矩阵提出来
coef_10 <- summary(glm_10)$coefficients$cond
coef_10 <- as.data.frame(coef_10)

coef_20 <- summary(glm_20)$coefficients$cond
coef_20 <- as.data.frame(coef_20)

coef_50 <- summary(glm_50)$coefficients$cond
coef_50 <- as.data.frame(coef_50)

#rownames(coef_10)
#[1] "(Intercept)"        "minpd_10"           "avepd_10"           "totpd_10"          
#[5] "SRA"                "DBH2"               "CBD_10"             "shannon_div_10"    
#[9] "RDi"                "shannon_div_10:RDi"
#再试试
#比较它们
#intercept
t_test_10_20_intercept <- t.test(c(coef_10["(Intercept)","Estimate"], coef_20["(Intercept)","Estimate"]), mu = 0)
t_test_10_50_intercept <- t.test(c(coef_10["(Intercept)","Estimate"], coef_50["(Intercept)","Estimate"]), mu = 0)
t_test_20_50_intercept <- t.test(c(coef_20["(Intercept)","Estimate"], coef_50["(Intercept)","Estimate"]), mu = 0)
print(t_test_10_20_intercept)
print(t_test_10_50_intercept)
print(t_test_20_50_intercept)



#minpd
t_test_10_20_minpd <- t.test(c(coef_10["minpd_10","Estimate"], coef_20["minpd_20","Estimate"]), mu = 0)
t_test_10_50_minpd <- t.test(c(coef_10["minpd_10","Estimate"], coef_50["minpd_50","Estimate"]), mu = 0)
t_test_20_50_minpd <- t.test(c(coef_20["minpd_20","Estimate"], coef_50["minpd_50","Estimate"]), mu = 0)
print(t_test_10_20_minpd)
print(t_test_10_50_minpd)
print(t_test_20_50_minpd)

#avepd
t_test_10_20_avepd <- t.test(c(coef_10["avepd_10","Estimate"], coef_20["avepd_20","Estimate"]), mu = 0)
t_test_10_50_avepd <- t.test(c(coef_10["avepd_10","Estimate"], coef_50["avepd_50","Estimate"]), mu = 0)
t_test_20_50_avepd <- t.test(c(coef_20["avepd_20","Estimate"], coef_50["avepd_50","Estimate"]), mu = 0)
print(t_test_10_20_avepd)
print(t_test_10_50_avepd)
print(t_test_20_50_avepd)

#totpd
t_test_10_20_totpd <- t.test(c(coef_10["totpd_10","Estimate"], coef_20["totpd_20","Estimate"]), mu = 0)
t_test_10_50_totpd <- t.test(c(coef_10["totpd_10","Estimate"], coef_50["totpd_50","Estimate"]), mu = 0)
t_test_20_50_totpd <- t.test(c(coef_20["totpd_20","Estimate"], coef_50["totpd_50","Estimate"]), mu = 0)
print(t_test_10_20_totpd)
print(t_test_10_50_totpd)
print(t_test_20_50_totpd)

#SRA
t_test_10_20_SRA <- t.test(c(coef_10["SRA","Estimate"], coef_20["SRA","Estimate"]), mu = 0)
t_test_10_50_SRA <- t.test(c(coef_10["SRA","Estimate"], coef_50["SRA","Estimate"]), mu = 0)
t_test_20_50_SRA <- t.test(c(coef_20["SRA","Estimate"], coef_50["SRA","Estimate"]), mu = 0)
print(t_test_10_20_SRA)
print(t_test_10_50_SRA)
print(t_test_20_50_SRA)

#DBH2
t_test_10_20_DBH2 <- t.test(c(coef_10["DBH2","Estimate"], coef_20["DBH2","Estimate"]), mu = 0)
t_test_10_50_DBH2 <- t.test(c(coef_10["DBH2","Estimate"], coef_50["DBH2","Estimate"]), mu = 0)
t_test_20_50_DBH2 <- t.test(c(coef_20["DBH2","Estimate"], coef_50["DBH2","Estimate"]), mu = 0)
print(t_test_10_20_DBH2)
print(t_test_10_50_DBH2)
print(t_test_20_50_DBH2)

#CBD
t_test_10_20_CBD <- t.test(c(coef_10["CBD_10","Estimate"], coef_20["CBD_20","Estimate"]), mu = 0)
t_test_10_50_CBD <- t.test(c(coef_10["CBD_10","Estimate"], coef_50["CBD_50","Estimate"]), mu = 0)
t_test_20_50_CBD <- t.test(c(coef_20["CBD_20","Estimate"], coef_50["CBD_50","Estimate"]), mu = 0)
print(t_test_10_20_CBD)
print(t_test_10_50_CBD)
print(t_test_20_50_CBD)


#shannon_div
t_test_10_20_shannon_div <- t.test(c(coef_10["shannon_div_10","Estimate"], coef_20["shannon_div_20","Estimate"]), mu = 0)
t_test_10_50_shannon_div <- t.test(c(coef_10["shannon_div_10","Estimate"], coef_50["shannon_div_50","Estimate"]), mu = 0)
t_test_20_50_shannon_div <- t.test(c(coef_20["shannon_div_20","Estimate"], coef_50["shannon_div_50","Estimate"]), mu = 0)
print(t_test_10_20_shannon_div)
print(t_test_10_50_shannon_div)
print(t_test_20_50_shannon_div)

#RDi
t_test_10_20_RDi <- t.test(c(coef_10["RDi","Estimate"], coef_20["RDi","Estimate"]), mu = 0)
t_test_10_50_RDi <- t.test(c(coef_10["RDi","Estimate"], coef_50["RDi","Estimate"]), mu = 0)
t_test_20_50_RDi <- t.test(c(coef_20["RDi","Estimate"], coef_50["RDi","Estimate"]), mu = 0)
print(t_test_10_20_RDi)
print(t_test_10_50_RDi)
print(t_test_20_50_RDi)

#shannon_div:RDi
t_test_10_20_shannon_divRDi <- t.test(c(coef_10["shannon_div_10:RDi","Estimate"], coef_20["shannon_div_20:RDi","Estimate"]), mu = 0)
t_test_10_50_shannon_divRDi <- t.test(c(coef_10["shannon_div_10:RDi","Estimate"], coef_50["shannon_div_50:RDi","Estimate"]), mu = 0)
t_test_20_50_shannon_divRDi <- t.test(c(coef_20["shannon_div_20:RDi","Estimate"], coef_50["shannon_div_50:RDi","Estimate"]), mu = 0)
print(t_test_10_20_shannon_divRDi)
print(t_test_10_50_shannon_divRDi)
print(t_test_20_50_shannon_divRDi)

