library(tidyverse)
library(GenomicRanges)
library(msigdbr)
library(fgsea)
library(RobustRankAggreg)
library(forestplot)
library(pals)


load('misc/expd.rda')
load('misc/deg.rda')

expd <- lapply(expd, function(x) x$ID) %>%
  unlist() %>%
  unique()

t1 <- lapply(res, function(x) {
  na.omit(x) %>%
    dplyr::filter(ID %in% expd) %>%
    mutate(fc = log2FoldChangeAdj) %>%
    group_by(gene) %>%
    summarise(fc = mean(fc)) %>%
    arrange(fc) %>%
    pull(gene)
}) %>%
  {inner_join(aggregateRanks(.),
              aggregateRanks(lapply(., rev)),
              by = 'Name')} %>%
  mutate(sc = log10(Score.x) - log10(Score.y)) %>%
  dplyr::select(Name, sc) %>%
  deframe()

pways <- msigdbr(species = "Homo sapiens") %>%
  mutate(db = paste(gs_cat, gs_subcat, sep = '.') %>%
           sub('\\.$', '', .)) %>%
  split(., .$db) %>%
  lapply(function(db) {
    split(db$gene_symbol, db$gs_name)
  })


fres <- pblapply(pways, function(x) {
  fgsea(pathways = x,
        stats = t1,
        minSize = 15,
        maxSize = 500)
})

fresc <- names(pways) %>%
  setNames(., .) %>%
  pblapply(function(x) {
    r <- fres[[x]]
    clps <- collapsePathways(r[padj < .05], pways[[x]], t1)
    r[pathway %in% clps$mainPathways][order(-NES),]
  })

tbltxt <- fresc$H %>%
  as_tibble() %>%
  arrange(-NES) %>%
  dplyr::select(pathway, padj, NES, size) %>%
  mutate(pathway = sub('^HALLMARK_', '', pathway),
         padj = sprintf('%.2g', padj),
         NES = sprintf('%.2g', NES)) %>%
  rbind(c('Pathway', 'FDR', 'NES', 'Size'), .) %>%
  as.matrix()

getWhisks <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, c(0.25, 0.75), na.rm = T)
  data.frame(lower = qs[1], upper = qs[2], middle = median(x, na.rm = T),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))], na.rm = T),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))], na.rm = T)) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}

d <- fresc$H %>%
  as_tibble() %>%
  mutate(nm = paste(pathway, padj, NES, sep = '>')) %>%
  split(., .$nm) %>%
  lapply(function(x) {
    lapply(res, function(y) {
      pways$H[[x$pathway]] %>%
        {y[y$gene %in% ., ]} %>%
        mutate(e = stat) %>%
        group_by(gene) %>%
        summarise(fc = -mean(log2FoldChange),
                  fa = -mean(log2FoldChangeAdj),
                  e = -mean(e))
    }) %>% bind_rows(.id = 'comp')
  }) %>% bind_rows(.id = 'nm') %>%
  separate(nm, c('path', 'padj', 'NES'), '>') %>%
  mutate(padj = as.numeric(padj),
         NES = as.numeric(NES)) %>%
  group_by(comp, path, NES) %>%
  do(getWhisks(.$fc)) %>%
  ungroup() %>%
  arrange(-NES) %>%
  split(., .$comp) %>%
  {list(mean = lapply(., `[[`, 'middle'),
        lower = lapply(., `[[`, 'lower'),
        upper = lapply(., `[[`, 'upper'))} %>%
  lapply(function(x) {
    rbind(c(NA, NA), bind_cols(x))
  })


clrs <- tableau20(3)[c(1,3)]
pdf('figs/5c.pdf', onefile = F, width = 15)
forestplot(tbltxt,
           fn.ci_norm = c(fpDrawCircleCI, fpDrawCircleCI),
           boxsize = 0.2,
           mean = d$mean,
           upper = d$upper,
           lower = d$lower,
           col = fpColors(box = clrs, lines = clrs),
           legend = c('WT vs MT', 'WT vs KO'),
           is.summary = c(TRUE, rep(FALSE, nrow(tbltxt)-1)),
           xlab = 'LFC of genes in pathway (KO or MT - WT)',
           graphwidth = unit(10, 'cm'),
           xticks = -2:3,
           lwd.ci = 1.5,
           ci.vertices = T)
dev.off()
