library(tidyverse)
library(tximport)
library(DESeq2)

load("misc/tx2gene.rda")
load('misc/g33.rda')
e2g <- g33[g33$type == 'gene'] %>%
  {setNames(.$gene_name, .$ID)}

samps <- list.files('data', pattern = '.sf.txt.gz$', full.names = T) %>%
  setNames(., sub('.*/(.*).sf.*', '\\1', .)) 

mdat <- data.frame(samp = names(samps)) %>%
  mutate(type = case_when(
    grepl('KO',samp) ~ 'KO',
    grepl('^S|^B', samp) ~ 'MT',
    TRUE ~ 'WT'
  )) %>%
  dplyr::filter(!(samp %in% c('Cal27_KO17', 'Det562_KO1', 'Det562_KO4'))) %>%
  mutate(line = factor(sub('_.*', '', samp)),
         type = factor(type)) %>%
  column_to_rownames("samp")

res <- c('MT', 'KO') %>%
  setNames(., .) %>%
  lapply(function(x) {
    m <- mdat[mdat$type %in% c('WT', x),] %>% droplevels()
    s <- samps[names(samps) %in% rownames(m)] %>%
      tximport(type = "salmon", tx2gene = tx2gene)
    s$length[s$length == 0] <- 1
    dds <- if (x == 'KO') {
      DESeqDataSetFromTximport(s, m, ~line + type)
    } else {
      DESeqDataSetFromTximport(s, m, ~type)
    }
    dds <- DESeq(dds)
    res <- results(dds)
    rnm <- resultsNames(dds)[length(resultsNames(dds))]
    print(rnm)
    lfcShrink(dds, rnm, res = res,
              type = "apeglm", lfcThreshold = 0.5) %>%
      data.frame() %>%
      rownames_to_column("gene") %>%
      left_join(rownames_to_column(data.frame(res), "gene"),
                by = "gene") %>%
      dplyr::select(-baseMean.y) %>%
      `colnames<-`(c('ID', 'baseMean', 'log2FoldChangeAdj', 'lfcAdjSE',
                     'svalue', 'log2FoldChange','lfcSE','stat','pvalue','padj')) %>%
      mutate(gene = e2g[ID]) %>%
      dplyr::select(ID, gene, baseMean, log2FoldChange, lfcSE,
                    stat, pvalue, padj, log2FoldChangeAdj, lfcAdjSE, svalue)
  })

save(res, file = 'misc/deg.rda')
