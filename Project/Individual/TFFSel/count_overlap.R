library(dplyr)
library(data.table)
library(purrr)
library(GenomicRanges)
library(regioneR)
library(liftOver)
library(tidyr)
library(writexl)

library(ggplot2)
library(ggpubr)

ch = import.chain("hg38ToHg19.over.chain")
load("cms.rda")
load("dhs.rda")
load("footprints.rda")

footprints0.01 = footprints %>%
  filter(conf < 0.01)

dhs = GRanges(seqnames = dhs$chromosome, 
              ranges = IRanges(start = dhs$start, 
                               end = dhs$end)) %>%
  liftOver(chain = ch) %>% 
  unlist()

footprints0.01 = GRanges(seqnames = footprints0.01$chromosome, 
                         ranges = IRanges(start = footprints0.01$start, 
                                          end = footprints0.01$end)) %>%
  liftOver(chain = ch) %>% 
  unlist()

cms = cms %>%
  mutate(bin100kb = paste0(CHROM, "_", floor(BP/1e5))) %>%
  mutate(bin200kb = paste0(CHROM, "_", floor(BP/2e5))) %>%
  mutate(bin50kb = paste0(CHROM, "_", floor(BP/5e4))) %>%
  mutate(bin10kb = paste0(CHROM, "_", floor(BP/1e4)))

cms_1b = cms %>%
  mutate(chromosome = paste0("chr", CHROM)) %>%
  mutate(start_pos = BP) %>%
  mutate(end_pos = BP) %>%
  dplyr::select(chromosome, start_pos, end_pos, cms_score)

count_overlap = function(gr1, gr2){
  overlaps = as.data.frame(findOverlaps(gr1, gr2))
  return(dim(overlaps)[1])
}

calc_enrichment = function(gr, target){
  n_overlap = count_overlap(target, gr)
  fold_enrichment_list = vector()
  cat("Generating random regions...\n")
  for(i in 1:30){
    gtarget_random = randomizeRegions(target, genome = "hg19")
    random_overlap = count_overlap(gtarget_random, gr)
    fold_enrichment = n_overlap/random_overlap
    fold_enrichment_list = c(fold_enrichment_list, fold_enrichment)
  }
  return(fold_enrichment_list)
}

overlap_count_df = vector()
for(cms_threshold in 0:20){
  cms_1b_subset = filter(cms_1b, cms_score >= cms_threshold)
  cms_1b_subset = GRanges(seqnames = cms_1b_subset$chromosome, 
                          ranges = IRanges(start = cms_1b_subset$start_pos, 
                                           end = cms_1b_subset$end_pos),
                          cms_score = cms_1b_subset$cms_score)
  overlap_count_dhs = count_overlap(dhs, cms_1b_subset)
  overlap_count_footprints = count_overlap(footprints0.01, cms_1b_subset)
  overlap_count_df = rbind.data.frame(overlap_count_df,
                                      data.frame(
                                        cms_threshold = cms_threshold,
                                        dhs_count = overlap_count_dhs,
                                        footprints_count = overlap_count_footprints
                                      )
  )
}

overlap_count_df = overlap_count_df %>%
  mutate(dhs_count_percent = dhs_count/max(dhs_count)) %>%
  mutate(footprints_count_percent = footprints_count/max(footprints_count))

write_xlsx(overlap_count_df, path = "overlap_count.xlsx")
