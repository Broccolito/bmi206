library(dplyr)
library(data.table)
library(purrr)

dhs = read.csv(gzfile("wgEncodeAwgDnaseMasterSites.txt.gz"), 
               sep = "\t", header = FALSE)

dhs = dhs %>%
  select(V2, V3, V4)
names(dhs) = c("chromosome", "start", "end")

gwas = read.csv(gzfile("gwas_catalog.tsv.gz"), 
                sep = "\t", header = TRUE)

gwas = gwas %>%
  select(CHR_ID, CHR_POS, DISEASE.TRAIT, REPORTED.GENE.S.)
names(gwas) = c("chromosome", "pos", "trait", "gene_symbol")
gwas$chromosome = paste0("chr", gwas$chromosome)

file_list = file.path("consensus_footprints",
                      list.files(path = "consensus_footprints/",
                                 pattern = "bed.gz"))

bed = vector()
for(f in file_list){
  cat(paste0("Reading ", f, "...\n"))
  chr_file = read.csv(gzfile(f), sep = "\t", header = FALSE)
  bed = rbind.data.frame(bed, chr_file)
}
footprints = bed
names(footprints) = c('chromosome', 'start', 'end', 'identifier',
                      'mean_signal', 'num_samples', 'num_fps',
                      'width', 'summit_pos', 'core_start', 
                      'core_end', 'motif_clusters')

footprints = footprints %>%
  select(-motif_clusters, -(num_samples:width)) %>%
  mutate(p = 1 - (2^(-mean_signal))) %>%
  mutate(conf = 1 - p) %>%
  mutate(chr_num = gsub(pattern = "chr", replacement = "", chromosome)) %>%
  mutate(chr_num = as.numeric(chr_num)) %>%
  arrange(chr_num, chromosome) %>%
  select(chromosome, start, end, conf)

# Save files
save(dhs, file = "dhs.rda")
save(gwas, file = "gwas.rda")
save(footprints, file = "footprints.rda")

