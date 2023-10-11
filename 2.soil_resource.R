#####
#计算土壤的resource diversity、resource evenness、resource richness
#然后求RR和RE的平均值

library(dplyr)
library(vegan)

cy_resource <- soil_pred %>% 
  select(soc:ph) %>% 
  mutate_all(decostand, "max") # scale the data by "max" method, just testing
# other methods are OK

#
Si <- rowSums(pi * cy_resource ^ 2/ncol(cy_resource))

Li_part1 <- apply(cy_resource, 1, function(x) 2 * (max(x) - min(x))jiu)
Li_part2 <- rowSums(2 * pi * cy_resource)/ncol(cy_resource)

Li <- Li_part1 + Li_part2

#
RRi <- Si / max(Si)

REi <- 4 * pi * Si / (Li ^ 2) 

RDi <- sqrt(RRi * REi)

soil_pred <- cbind(soil_pred, RRi, REi, RDi)

rm(Si, RRi, REi, RDi,Li_part1, Li_part2, Li)
