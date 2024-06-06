library(dplyr)
library(data.table)

cms = fread("tibetan_cms_components.csv")

cms = filter(cms, cms_score >= 0)

save(cms, file = "cms_tibetan.rda")