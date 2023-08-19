#手动计算α多样性
library(dplyr)

shannon_20 <- c();sinpson_20 <- c();berger_20 <- c()

for (i in 1:dim(sp_loc)[1]) {
  
  #这样就筛选出了一个群落
  # 20m circle
  hsd_sub <- hsd_alive_singlebr %>% 
    filter(sqrt((GX - sp_loc$GX[i])^2 + (GY - sp_loc$GY[i])^2) <= 20)
  
  #计算这个群落的阿尔法多样性
  #Shannon
  # 计算每个scientific.name的数量占总行数的比例以及新列ca_shannon
  hsd_sub_ratios <- hsd_sub %>%
    count(scientific.name, name = "count") %>%
    mutate(Ratio = count / nrow(hsd_sub),
           ca_shannon = -Ratio * log(Ratio))
  
  shannon_20[i] <- sum(hsd_sub_ratios$ca_shannon)
  #Sinpson
  hsd_sub_ratios <- hsd_sub_ratios %>%
    mutate(ca_sinpson = Ratio * Ratio)
  
  sinpson_20[i] <- sum(hsd_sub_ratios$ca_sinpson)
  
  #Berger-Parker
  berger_20[i] <- hsd_sub_ratios %>%
    filter(count == max(count)) %>%
    slice(1) %>%
    pull(Ratio)
}

sp_loc <- sp_loc %>% 
  mutate(
    shannon_20 = sinpson_20,
    sinpson_20 = sinpson_20,
    berger_20  = berger_20
  )