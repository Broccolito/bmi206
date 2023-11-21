library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

d = read.csv("enrichment_df.csv")

d = d %>%
  select(-footprints0.0001_enrichment)

offset = mean(d$dhs_enrichment)
d = d %>%
  mutate(dhs_enrichment = dhs_enrichment / offset) %>%
  mutate(footprints0.05_enrichment = footprints0.05_enrichment / offset) %>%
  mutate(footprints0.01_enrichment = footprints0.01_enrichment / offset) %>%
  mutate(footprints0.001_enrichment = footprints0.001_enrichment / offset)

names(d) = c("Within DHS,\noutside footprints", "0.05", "0.01", "0.001")

d = d %>% 
  gather(key = group, value = "fold")

d$group = factor(d$group, levels = c("Within DHS,\noutside footprints", "0.05", "0.01", "0.001"))

plt = ggplot(data = d, aes(x = group, y = fold)) + 
  geom_boxplot(aes(fill = group), outlier.shape = NA) + 
  scale_fill_manual(values = c("gray", "#ffbaba", "#ff7b7b", "#a70000")) + 
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "#F8766D") + 
  stat_compare_means(ref.group = "Within DHS,\noutside footprints", 
                     label = "p.signif", method = 't.test',
                     symnum.args=list(
                       cutpoints = c(0, 0.05, 1), 
                       symbols = c("*", "ns")),
                     family = "mono") + 
  xlab("") + 
  ylab("Fold-enrichment over\nsampled 1KGP SNVs") + 
  ylim(c(0.9, 2.2)) + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none"); plt

ggsave(filename = "Figure 6a.png", device = "png", dpi = 1200, 
       width = 3, height = 4.5)
