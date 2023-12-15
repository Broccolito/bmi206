library(dplyr)
library(data.table)
library(purrr)
library(GenomicRanges)
library(regioneR)
library(liftOver)
library(tidyr)
library(writexl)
library(readxl)
library(ggplot2)
library(ggpubr)

enrichment = read_xlsx(path = "enrichment.xlsx")
enrichment$`CMS Score` = factor(enrichment$`CMS Score`, 
                                levels = unique(enrichment$`CMS Score`))


stats = enrichment %>% 
  split(.$`CMS Score`) %>%
  map(function(x){
    ttest_data = data.frame(
      just_one = 1,
      dhs_enrichment = filter(x, group == "DHS")$fold_enrichment,
      footprint_enrichment = filter(x, group == "TF Footprints")$fold_enrichment
    )
    # return(ttest_data)
    dhs2baseline_beta = t.test(ttest_data$dhs_enrichment, ttest_data$just_one)$estimate %>% diff()
    footprint2baseline_beta = t.test(ttest_data$footprint_enrichment, ttest_data$just_one)$estimate %>% diff()
    dhs2footprint_beta = t.test(ttest_data$footprint_enrichment, ttest_data$dhs_enrichment)$estimate %>% diff()
    
    dhs2baseline_p = t.test(ttest_data$dhs_enrichment, ttest_data$just_one)$p.value
    footprint2baseline_p = t.test(ttest_data$footprint_enrichment, ttest_data$just_one)$p.value
    dhs2footprint_p = t.test(ttest_data$footprint_enrichment, ttest_data$dhs_enrichment)$p.value
    
    data.frame(
      cms_threshold = x$`CMS Score`[1],
      dhs2baseline_beta, dhs2baseline_p,
      footprint2baseline_beta, footprint2baseline_p,
      dhs2footprint_beta, dhs2footprint_p
    )
    
    
  }) %>%
  purrr::reduce(rbind.data.frame)

rownames(stats) = NULL

write_xlsx(stats, path = "comparison_stats.xlsx")



