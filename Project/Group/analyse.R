library(dplyr)
library(data.table)
library(purrr)
library(regioneR)

load("dhs.rda")
load("footprints.rda")
load("gwas.rda")

footprints0.05 = footprints %>%
  filter(conf < 0.05)

footprints0.01 = footprints %>%
  filter(conf < 0.01)

footprints0.001 = footprints %>%
  filter(conf < 0.001)

footprints0.0001 = footprints %>%
  filter(conf < 0.0001)

gwas = gwas %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(!is.na(pos))

dhs = GRanges(seqnames = dhs$chromosome, 
              ranges = IRanges(start = dhs$start, 
                               end = dhs$end))

footprints0.05 = GRanges(seqnames = footprints0.05$chromosome, 
                         ranges = IRanges(start = footprints0.05$start, 
                                          end = footprints0.05$end))

footprints0.01 = GRanges(seqnames = footprints0.01$chromosome, 
                         ranges = IRanges(start = footprints0.01$start, 
                                          end = footprints0.01$end))

footprints0.001 = GRanges(seqnames = footprints0.001$chromosome, 
                          ranges = IRanges(start = footprints0.001$start, 
                                           end = footprints0.001$end))

footprints0.0001 = GRanges(seqnames = footprints0.0001$chromosome, 
                           ranges = IRanges(start = footprints0.0001$start, 
                                            end = footprints0.0001$end))

gwas = GRanges(seqnames = gwas$chromosome, 
               ranges = IRanges(start = gwas$pos, 
                                end = gwas$pos))

count_overlap = function(gr1, gr2){
  overlaps = as.data.frame(findOverlaps(gr1, gr2))
  return(dim(overlaps)[1])
}


calc_enrichment = function(gr, gwas){
  n_overlap = count_overlap(gwas, gr)
  fold_enrichment_list = vector()
  cat("Generating random regions...\n")
  for(i in 1:100){
    gwas_random = randomizeRegions(gwas, genome = "hg38")
    random_overlap = count_overlap(gwas_random, gr)
    fold_enrichment = n_overlap/random_overlap
    fold_enrichment_list = c(fold_enrichment_list, fold_enrichment)
  }
  return(fold_enrichment_list)
}


dhs_enrichment = calc_enrichment(gr = dhs, gwas = gwas)
footprints0.05_enrichment = calc_enrichment(gr = footprints0.05, gwas = gwas)
footprints0.01_enrichment = calc_enrichment(gr = footprints0.01, gwas = gwas)
footprints0.001_enrichment = calc_enrichment(gr = footprints0.001, gwas = gwas)
footprints0.0001_enrichment = calc_enrichment(gr = footprints0.0001, gwas = gwas)


enrichment_df = tibble(
  dhs_enrichment = dhs_enrichment, 
  footprints0.05_enrichment = footprints0.05_enrichment,
  footprints0.01_enrichment = footprints0.01_enrichment, 
  footprints0.001_enrichment = footprints0.001_enrichment, 
  footprints0.0001_enrichment = footprints0.01_enrichment
)

write.csv(enrichment_df, file = "enrichment_df.csv", row.names = FALSE)



