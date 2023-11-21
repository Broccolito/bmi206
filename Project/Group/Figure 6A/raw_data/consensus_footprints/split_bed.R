library(dplyr)
library(data.table)
library(purrr)

d = fread("consensus_footprints_and_collapsed_motifs_hg38.bed")

d %>%
  split(.$V1) %>%
  map(function(x){
    cat(paste0("Saving ", x$V1[1], "...\n"))
    fwrite(x, sep = "\t", col.names = FALSE, file = paste0("consensus_footprints_", x$V1[1], ".bed"))
  })
