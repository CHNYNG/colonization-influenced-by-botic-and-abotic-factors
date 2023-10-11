#####
#计算土壤的resource diversity、resource evenness、resource richness
#然后求RR和RE的平均值

library(dplyr)
library(vegan)

cy_resource <- soil_pred %>% 
  select(soc:ph,moisture) %>% 
  mutate_all(decostand, "max") # scale the data by "max" method, just testing
# other methods are OK

#
Si <- rowSums(pi * cy_resource ^ 2/ncol(cy_resource))

Li_part1 <- apply(cy_resource, 1, function(x) 2 * (max(x) - min(x)))
Li_part2 <- rowSums(2 * pi * cy_resource)/ncol(cy_resource)

Li <- Li_part1 + Li_part2

#
RRi <- Si / max(Si)

REi <- 4 * pi * Si / (Li ^ 2) 

RDi <- sqrt(RRi * REi)



soil_pred <- cbind(soil_pred, RRi, REi, RDi)

# 假设您有一个DataFrame soil_pred，并且要处理其中的RDi列

# 划分为三组并分配标签
# 假设您有一个DataFrame soil_pred，并且要处理其中的RDi列

# 利用ifelse将数据分为三组，并分配标签
soil_pred$resource <- ifelse(soil_pred$RDi <= quantile(soil_pred$RDi, 1/3), "low",
                             ifelse(soil_pred$RDi <= quantile(soil_pred$RDi, 2/3), "medium", "high"))



rm(Si, RRi, REi, RDi,Li_part1, Li_part2, Li)
