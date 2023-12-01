library(dplyr)
library(data.table)

# ov = fread("overall_network_property.csv")

pt = fread("parent_PPI node class.csv")
ms = fread("MS Significant.csv")

ms = ms %>%
  filter(`Node Class` != "")

pt = pt %>%
  filter(`Node Class` != "")

table(ms$`Node Class`)/dim(ms)[1]
table(pt$`Node Class`)/dim(pt)[1]