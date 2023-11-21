library(dplyr)
library(data.table)
library(purrr)
library(regioneR)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)

load("dhs.rda")
load("footprints.rda")
load("gwas.rda")

trait_cat = read_excel("10cat_traits.xlsx")
trait_cat = trait_cat %>%
  split(.$Category) %>%
  map(function(x){
    tibble(
      trait = unlist(strsplit(x$Traits, split = ", "))
    ) %>%
      mutate(category = x$Category[1])
  }) %>%
  purrr::reduce(rbind.data.frame)

gwas = gwas %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(!is.na(pos)) %>%
  left_join(trait_cat, relationship = "many-to-many", by = "trait")

category = table(gwas$category)
category = tibble(
  name = names(category),
  count = as.numeric(category)
) %>%
  arrange(desc(count))

category$name = factor(category$name, levels = category$name)

plt = ggplot(data = category, aes(x = name, y = count)) +
  geom_bar(stat = "identity", fill = "gray", color = "black") +
  geom_text(aes(label = count), vjust = -0.5, size = 3) +
  xlab("") +
  ylab("SNP Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l = 0.2, unit = "in"))
plt

ggsave(filename = "SNP count by trait category.png", dpi = 1200, plot = plt,
       width = 3.5, height = 5, bg = "white")

# Subsetting GWAS findings by category
gwas %>%
  filter(!is.na(category)) %>%
  split(.$category) %>%
  map(function(x){
    variable_name = paste0("gwas_", gsub(pattern = " ", replacement = "_", tolower(x$category[1])))
    assign(variable_name, x, envir = globalenv())
  })

get_enrichment_df = function(gwas){
  
  load("dhs.rda")
  load("footprints.rda")
  
  footprints0.05 = footprints %>%
    filter(conf < 0.05)
  
  footprints0.01 = footprints %>%
    filter(conf < 0.01)
  
  footprints0.001 = footprints %>%
    filter(conf < 0.001)
  
  footprints0.0001 = footprints %>%
    filter(conf < 0.0001)
  
  gwas_category = gwas$category[1]
  
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
    for(i in 1:10){
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
    category = gwas_category,
    dhs_enrichment = dhs_enrichment, 
    footprints0.05_enrichment = footprints0.05_enrichment,
    footprints0.01_enrichment = footprints0.01_enrichment, 
    footprints0.001_enrichment = footprints0.001_enrichment, 
    footprints0.0001_enrichment = footprints0.01_enrichment
  )
  
  return(enrichment_df)
  
}

rm(gwas)
load("gwas.rda")
gwas = gwas %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(!is.na(pos)) %>%
  left_join(trait_cat, relationship = "many-to-many", by = "trait")
enrichment_overall = get_enrichment_df(gwas = gwas)
enrichment_overall$category = "Overall"
enrichment_behavioral_and_psychological = get_enrichment_df(gwas = gwas_behavioral_and_psychological)
enrichment_blood_and_immune_system = get_enrichment_df(gwas = gwas_blood_and_immune_system)
enrichment_cardiovascular_and_circulatory = get_enrichment_df(gwas = gwas_cardiovascular_and_circulatory)
enrichment_endocrine_and_hormonal = get_enrichment_df(gwas = gwas_endocrine_and_hormonal)
enrichment_genetic_and_developmental = get_enrichment_df(gwas = gwas_genetic_and_developmental)
enrichment_metabolic_and_lipidomic = get_enrichment_df(gwas = gwas_metabolic_and_lipidomic)
enrichment_nutritional_and_dietary = get_enrichment_df(gwas = gwas_nutritional_and_dietary)
enrichment_physical_anthropometrics = get_enrichment_df(gwas = gwas_physical_anthropometrics)
enrichment_respiratory_and_pulmonary = get_enrichment_df(gwas = gwas_respiratory_and_pulmonary)

enrichment_df = rbind.data.frame(
  enrichment_overall,
  enrichment_behavioral_and_psychological,
  enrichment_blood_and_immune_system,
  enrichment_cardiovascular_and_circulatory,
  enrichment_endocrine_and_hormonal,
  enrichment_genetic_and_developmental,
  enrichment_metabolic_and_lipidomic,
  enrichment_nutritional_and_dietary,
  enrichment_physical_anthropometrics,
  enrichment_respiratory_and_pulmonary
)

fold_diff = enrichment_df %>%
  split(.$category) %>%
  map(function(x){
    fold_diff_df = tibble(
      category = x$category[1],
      fold_diff = x$footprints0.01_enrichment / x$dhs_enrichment
    ) %>%
      filter(is.finite(fold_diff))
    
    upper_limit = mean(fold_diff_df$fold_diff) + 3*sd(fold_diff_df$fold_diff)
    
    fold_diff_df = fold_diff_df %>%
      filter(fold_diff <= upper_limit)
    
    return(fold_diff_df)
  }) %>%
  purrr::reduce(rbind.data.frame)

fold_diff_ranked = fold_diff %>%
  split(.$category) %>%
  map(function(x){
    mean(x$fold_diff[is.finite(x$fold_diff)], na.rm = TRUE)
  })

fold_diff_ranked = tibble(
  category = names(fold_diff_ranked),
  fold = unlist(fold_diff_ranked)
) %>%
  arrange(desc(fold))

fold_diff$category = factor(fold_diff$category, 
                            levels = c("Overall", fold_diff_ranked$category[fold_diff_ranked$category != "Overall"]) )

plt = ggplot(data = fold_diff, aes(x = category, y = fold_diff)) +
  # geom_point() + 
  geom_boxplot(fill = "gray") +
  xlab("") + 
  ylab("TF Footprint-DHS \nEnrichment Difference") + 
  stat_compare_means(ref.group = "Overall", 
                     label = "p.signif", method = 't.test',
                     symnum.args=list(
                       cutpoints = c(0, 0.05, 1), 
                       symbols = c("*", "NS")),
                     family = "mono") + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(l = 0.2, unit = "in"),
        legend.position = "none")

ggsave(filename = "Differential enrichment.png", dpi = 1200, plot = plt,
       device = "png", width = 4, height = 5)

aov_result = aov(formula = fold_diff ~ category, data = fold_diff)
TukeyHSD(aov_result)

write.csv(enrichment_df, file = "enrichment_df.csv", row.names = FALSE)



