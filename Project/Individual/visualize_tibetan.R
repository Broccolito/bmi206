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

enrichment = read_xlsx(path = "enrichment_tibetan.xlsx")
enrichment$`CMS Score` = factor(enrichment$`CMS Score`, 
                                levels = unique(enrichment$`CMS Score`))

plt = ggplot(data = enrichment, 
             aes(x = `CMS Score`,
                 y = fold_enrichment,
                 fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#e68613", "#8494ff")) +
  labs(fill = "") + 
  geom_hline(aes(yintercept = 1), color = "coral1", linetype = "dashed") + 
  xlab(bquote(CMS[BF])) + 
  ylab("Fold Enrichment") + 
  ylim(c(0.8, 4)) + 
  theme_minimal() + 
  theme(text = element_text(size = 15), legend.position = "top"); plt

ggsave(filename = "Enrichment plot Tibetan.png", device = "png", dpi = 1200,
       width = 6, height = 4, bg = "white")
