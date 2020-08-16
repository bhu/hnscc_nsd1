library(tidyverse)
library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
library(pals)


load('misc/dac.rda')

pannos <- c('WT > KO', 'KO > WT', 'All') %>%
  setNames(., .) %>%
  lapply(function(x) {
    pks <- if (x == 'All') {
      dba.peakset(d, bRetrieve = T)
    } else if (x == 'WT > KO') {
      d.o[d.o$Fold < 0]
    } else if (x == 'KO > WT') {
      d.o[d.o$Fold > 0]
    } 
    annotatePeak(pks, tssRegion = c(-3000, 3000),
                 TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                 annoDb = "org.Hs.eg.db")
  })

anns <- lapply(pannos, function(x) x@annoStat) %>%
  bind_rows(.id = "comp") %>%
  mutate(ttl = 'Distribution of genomic region annotations',
         comp = fct_inorder(comp))

clrs <- c(stepped3(4)[c(1,2,4)], tableau20(12)[c(3:9,11)]) %>%
  setNames(c("Promoter (<=1kb)",
             "Promoter (1-2kb)",
             "Promoter (2-3kb)",
             "5' UTR",
             "3' UTR",
             "1st Exon",
             "Other Exon",
             "1st Intron",
             "Other Intron",
             "Downstream (<=300)",
             "Distal Intergenic"))

brks <- c("Promoter (<=1kb)",
          "Promoter (1-2kb)",
          "Promoter (2-3kb)",
          "5' UTR",
          "3' UTR",
          "Downstream (<=300)",
          "1st Exon",
          "Other Exon",
          "Distal Intergenic",
          "1st Intron",
          "Other Intron")
p <- ggplot(anns, aes(x = comp, y = Frequency)) +
  geom_col(aes(fill = Feature)) +
  coord_flip() +
  scale_fill_manual(values = clrs,
                    breaks = brks,
                    guide = guide_legend(byrow = F)) +
  scale_y_continuous(expand = expansion(0)) +
  #facet_grid(. ~ ttl) +
  ylab('Proportion (%)') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = margin(r = 10, l = 5))

ggsave('figs/3e.pdf', p, height = 2.8, width = 7)

