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
load("cms_tibetan.rda")
load("dhs.rda")
load("footprints.rda")

# footprints0.05 = footprints %>%
#   filter(conf < 0.05)

footprints0.01 = footprints %>%
  filter(conf < 0.01)

# footprints0.001 = footprints %>%
#   filter(conf < 0.001)
# 
# footprints0.0001 = footprints %>%
#   filter(conf < 0.0001)

dhs = GRanges(seqnames = dhs$chromosome, 
              ranges = IRanges(start = dhs$start, 
                               end = dhs$end)) %>%
  liftOver(chain = ch) %>% 
  unlist()

# footprints0.05 = GRanges(seqnames = footprints0.05$chromosome, 
#                          ranges = IRanges(start = footprints0.05$start, 
#                                           end = footprints0.05$end)) %>%
#   liftOver(chain = ch) %>% 
#   unlist()

footprints0.01 = GRanges(seqnames = footprints0.01$chromosome, 
                         ranges = IRanges(start = footprints0.01$start, 
                                          end = footprints0.01$end)) %>%
  liftOver(chain = ch) %>% 
  unlist()

# footprints0.001 = GRanges(seqnames = footprints0.001$chromosome, 
#                           ranges = IRanges(start = footprints0.001$start, 
#                                            end = footprints0.001$end)) %>%
#   liftOver(chain = ch) %>% 
#   unlist()
# 
# footprints0.0001 = GRanges(seqnames = footprints0.0001$chromosome, 
#                            ranges = IRanges(start = footprints0.0001$start, 
#                                             end = footprints0.0001$end)) %>%
#   liftOver(chain = ch) %>% 
#   unlist()

names(cms)[1] = "CHROM"
names(cms)[2] = "BP"

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

# cms_100kb = cms %>%
#   split(.$bin100kb) %>%
#   map(function(x){
#     head(arrange(x, desc(cms_score)),1)
#   }) %>%
#   purrr::reduce(rbind.data.frame) %>%
#   mutate(chromosome = paste0("chr", CHROM)) %>%
#   mutate(start_pos = floor(BP/1e5)*1e5) %>%
#   mutate(end_pos = start_pos + 1e5) %>%
#   dplyr::select(chromosome, start_pos, end_pos, cms_score) %>%
#   arrange(chromosome, start_pos)

# cms_200kb = cms %>%
#   split(.$bin200kb) %>%
#   map(function(x){
#     head(arrange(x, desc(cms_score)),1)
#   }) %>%
#   purrr::reduce(rbind.data.frame) %>%
#   mutate(chromosome = paste0("chr", CHROM)) %>%
#   mutate(start_pos = floor(BP/2e5)*2e5) %>%
#   mutate(end_pos = start_pos + 2e5) %>%
#   dplyr::select(chromosome, start_pos, end_pos, cms_score) %>%
#   arrange(chromosome, start_pos)

# cms_50kb = cms %>%
#   split(.$bin50kb) %>%
#   map(function(x){
#     head(arrange(x, desc(cms_score)),1)
#   }) %>%
#   purrr::reduce(rbind.data.frame) %>%
#   mutate(chromosome = paste0("chr", CHROM)) %>%
#   mutate(start_pos = floor(BP/5e4)*5e4) %>%
#   mutate(end_pos = start_pos + 5e4) %>%
#   dplyr::select(chromosome, start_pos, end_pos, cms_score) %>%
#   arrange(chromosome, start_pos)

# cms_10kb = cms %>%
#   split(.$bin10kb) %>%
#   map(function(x){
#     head(arrange(x, desc(cms_score)),1)
#   }) %>%
#   purrr::reduce(rbind.data.frame) %>%
#   mutate(chromosome = paste0("chr", CHROM)) %>%
#   mutate(start_pos = floor(BP/1e4)*1e4) %>%
#   mutate(end_pos = start_pos + 1e4) %>%
#   dplyr::select(chromosome, start_pos, end_pos, cms_score) %>%
#   arrange(chromosome, start_pos)

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


enrichment_df = vector()
for(cms_threshold in 0:16){
  cat(paste0("Extracting markers with CMS Score >= ", cms_threshold, "...\n"))
  cms_1b_subset = filter(cms_1b, cms_score >= cms_threshold)
  cms_1b_subset = GRanges(seqnames = cms_1b_subset$chromosome, 
                          ranges = IRanges(start = cms_1b_subset$start_pos, 
                                           end = cms_1b_subset$end_pos),
                          cms_score = cms_1b_subset$cms_score)
  
  
  
  dhs_enrichment = calc_enrichment(gr = dhs, target = cms_1b_subset)
  footprints0.01_enrichment = calc_enrichment(gr = footprints0.01, target = cms_1b_subset)
  enrichment_df = rbind.data.frame(enrichment_df,
                                   data.frame(
                                     cms_threshold, 
                                     dhs_enrichment,
                                     footprints0.01_enrichment
                                   )
  )
}

# ttest = enrichment_df %>%
#   split(.$cms_threshold) %>%
#   map(function(x){
#     ttest = t.test(x$footprints0.01_enrichment, x$dhs_enrichment)
#     data.frame(
#       cms_threshold = x$cms_threshold[1],
#       beta = diff(ttest$estimate),
#       pvalue = ttest$p.value
#     )
#   }) %>%
#   purrr::reduce(rbind.data.frame)

# enrichment_stats = enrichment_df %>%
#   split(.$cms_threshold) %>%
#   map(function(x){
#     data.frame(
#       cms_threshold = x$cms_threshold[1],
#       dhs_mean = mean(x$dhs_enrichment),
#       dhs_se = sd(x$dhs_enrichment),
#       footprint_mean = mean(x$footprints0.01_enrichment),
#       footprint_se = sd(x$footprints0.01_enrichment)
#     )
#   }) %>%
#   purrr::reduce(rbind.data.frame)

names(enrichment_df) = c("CMS Score", "DHS", "TF Footprints")

enrichment_gathered = gather(enrichment_df, "group", "fold_enrichment", -`CMS Score`)
enrichment_gathered$`CMS Score` = factor(enrichment_gathered$`CMS Score`,
                                         levels = unique(enrichment_gathered$`CMS Score`))

# Normalize to NULL selection
enrichment_offset = mean(filter(enrichment_gathered, `CMS Score` == 0, group == "DHS")$fold_enrichment)
enrichment_gathered$fold_enrichment = enrichment_gathered$fold_enrichment / enrichment_offset

footprint_offset = mean(filter(enrichment_gathered, `CMS Score` == 0, group == "TF Footprints")$fold_enrichment) /
  mean(filter(enrichment_gathered, `CMS Score` == 0, group == "DHS")$fold_enrichment)
enrichment_gathered = enrichment_gathered %>%
  mutate(fold_enrichment = ifelse(group == "TF Footprints", fold_enrichment/footprint_offset, fold_enrichment))

# Save calcuated enrichemnt statistics
write_xlsx(enrichment_gathered, path = "enrichment_tibetan.xlsx")

