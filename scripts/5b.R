library(tidyverse)
library(GenomicRanges)
library(RRHO2)


load('misc/deg.rda')
load('misc/expd.rda')

expd <- lapply(expd, function(x) x$ID) %>%
  unlist() %>%
  unique()

res <- lapply(res, function(x) {
  x$es <- -log10(x$pvalue) * sign(x$stat)
  x
})

dv <- lapply(list(adj = 'svalue', fc = 'padj', stat = 'es'), function(sig) {
  okg <- lapply(res, function(x) {
    x$ID[!is.na(x[[sig]])]
  }) %>%
    Reduce(intersect, .) %>%
    unique() %>%
    {.[. %in% expd]}
  lapply(res, function(x) {
    x[x$ID %in% okg,] %>%
      mutate(val = -case_when(sig == 'svalue' ~ log2FoldChangeAdj,
                              sig == 'padj' ~ log2FoldChange,
                              T ~ es)) %>%
      dplyr::select(ID, val)
  })
})

rr2 <- RRHO2(dv$stat$MT, dv$stat$KO, plots = T, outputdir = 'figs',
             labels = c('MT', 'KO'), BY = T, alternative = 'split',
             method = "hyper")
list.files('figs', pattern = 'RRHO_', full.names = T) %>%
  sapply(file.remove)
file.rename('figs/RRHOMap_combined_MT_VS_KO.tiff',
            'figs/5b.tiff')
