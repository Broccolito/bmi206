library(dplyr)
library(data.table)

cms1 = fread("cms_components1.csv")
cms2 = fread("cms_components2.csv")

cms = rbind.data.frame(cms1, cms2)

cms = filter(cms, cms_score >= 0)

save(cms, file = "cms.rda")