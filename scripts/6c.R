library(tidyverse)
library(RRHO2)


load('misc/deg.rda')
r1 <- res$KO

load('tcga/deg.rda')
r2 <- res.tcga %>%
  dplyr::rename(ID = gene, gene = name)

okg <- intersect(na.omit(r1)$gene, na.omit(r2)$gene)
dres <- list(cell = r1, tcga = r2) %>%
  lapply(function(x) {
    x[match(okg, x$gene),]
  })

dv <- lapply(c(fc = 1, adj = 2, stat = 3), function(kind) {
  lapply(dres, function(x) {
    if (kind == 1) {
      dplyr::select(x, ID = gene, val = log2FoldChange)
    } else if (kind == 2) {
      dplyr::select(x, ID = gene, val = log2FoldChangeAdj)
    } else if (kind == 3) {
      mutate(x, es = sign(stat) * -log10(pvalue)) %>%
        dplyr::select(ID = gene, val = es)
    }
  })
})

rr2 <- RRHO2(dv$stat$cell, dv$stat$tcga, plots = T, outputdir = 'figs',
             labels = c('cell', 'TCGA'), BY = T, alternative = 'split',
             method = "hyper")
list.files('figs', pattern = 'RRHO_', full.names = T) %>%
  sapply(file.remove)
file.rename('figs/RRHOMap_combined_cell_VS_TCGA.tiff',
            'figs/6c.tiff')
