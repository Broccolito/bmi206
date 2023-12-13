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
glist = fread("glist_hg19.txt")

glist = glist %>%
  mutate(chromosome = paste0("chr", chromosome)) %>%
  mutate(start = start - 100000) %>%
  mutate(end = end + 100000)

glist = GRanges(seqnames = glist$chromosome, 
                ranges = IRanges(start = glist$start, 
                                 end = glist$end),
                gene_symbol = glist$gene_symbol)

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

cms_threshold = 12

cms_1b_subset = filter(cms_1b, cms_score >= cms_threshold)
cms_1b_subset = GRanges(seqnames = cms_1b_subset$chromosome, 
                        ranges = IRanges(start = cms_1b_subset$start_pos, 
                                         end = cms_1b_subset$end_pos),
                        cms_score = cms_1b_subset$cms_score)

cms_footprint_overlap = subsetByOverlaps(cms_1b_subset, footprints0.01)

overlap_hits = as.data.frame(findOverlaps(cms_footprint_overlap, glist))

glist = as.data.frame(glist)
cms_footprint_overlap = as.data.frame(cms_footprint_overlap)

cms_footprint_overlap$queryHits = 1:dim(cms_footprint_overlap)[1]
glist$subjectHits = 1:dim(glist)[1]
glist = glist %>%
  dplyr::select(gene_symbol, subjectHits)

annotated = cms_footprint_overlap %>%
  full_join(overlap_hits, by = "queryHits") %>%
  full_join(glist, by = "subjectHits")

annotated = annotated %>%
  split(.$gene_symbol) %>%
  map(function(x){
    head(arrange(x, desc(cms_score)),1)
  }) %>%
  purrr::reduce(rbind.data.frame) %>%
  arrange(desc(cms_score)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

annotated = annotated %>%
  dplyr::select(seqnames, start, cms_score, gene_symbol) %>%
  dplyr::rename("chromosome" = "seqnames") %>%
  dplyr::rename("pos" = "start") %>%
  filter(!is.na(cms_score))

annotated = annotated %>%
  mutate(chr_pos = paste0(chromosome, "_", pos)) %>%
  split(.$chr_pos) %>%
  map(function(x){
    genes = paste(x$gene_symbol, collapse = ", ")
    data.frame(
      chromosome = x$chromosome[1],
      pos = x$pos[1],
      cms_score = x$cms_score[1],
      genes = genes
    )
  }) %>%
  purrr::reduce(rbind.data.frame) %>%
  arrange(desc(cms_score))

write_xlsx(annotated, path = "annotated_cms.xlsx")

cytoscape_attribute = annotated %>%
  dplyr::select(gene_symbol, cms_score)
names(cytoscape_attribute) = c("Gene", "cms_score")
fwrite(cytoscape_attribute, sep = "\t", file = "cytoscape_attribute.txt")
