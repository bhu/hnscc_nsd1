library(data.table)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(DESeq2)
library(tximport)
library(limma)
library(patchwork)
library(pals)


load("misc/tx2gene.rda")
load('misc/g33.rda')
e2g <- g33[g33$type == 'gene'] %>%
  {setNames(.$gene_name, .$ID)}

load('tcga/DNA.rda')
load('tcga/RNA.rda')
RNA <- data
load('tcga/DNAme.rda')
DNAme <- data
load('data/mCG.450k.rda')
load('misc/barcodeData.rda')

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

m <- mdat[mdat$type %in% c('WT', 'KO'),] %>% droplevels()
s <- samps[names(samps) %in% rownames(m)] %>%
  tximport(type = "salmon", tx2gene = tx2gene)
s$length[s$length == 0] <- 1
dds <- DESeqDataSetFromTximport(s, m, ~type + line)
dds2 <- DESeqDataSet(RNA, design = ~patient)

d2 <- counts(dds2)
d1 <- counts(dds) %>%
  `rownames<-`(sub('\\..*', '', rownames(.))) %>%
  {.[rownames(.) %in% rownames(d2),]}

mm <- rownames(d1) %>%
  sort() %>%
  {cbind(d1[.,], d2[.,])}

md <- tibble(samp = colnames(mm)) %>%
  mutate(src = case_when(grepl('^TCGA', samp) ~ 'TCGA',
                         T ~ 'cell')) %>%
  column_to_rownames('samp')

dds3 <- DESeqDataSetFromMatrix(mm, md, ~src)
vsd <- vst(dds3, blind = T)
v <- assay(vsd)
cp <- cor(v, method = 'pearson')
cs <-  cor(v, method = 'spearman')

muts <- lapply(DNA, function(x) {
  x[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')] %>%
    dplyr::filter(Hugo_Symbol == 'NSD1')
}) %>% bind_rows(.id = 'Caller') %>%
  dplyr::filter(Variant_Classification != 'Silent') %>%
  pull(Tumor_Sample_Barcode) %>%
  sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', .) %>%
  unique()
clus <- read_xlsx('tcga/papillon2017.xlsx', skip = 1) %>%
  dplyr::select(bc = bcr_patient_barcode, grp = methylation_group) %>%
  dplyr::filter(grp != 'NA') %>%
  {setNames(.$grp, .$bc)}
smps <- mdat[mdat$type %in% c('WT', 'KO'),] %>%
  rownames_to_column('samp') %>%
  mutate(samp = as.character(samp)) %>%
  split(., .$type, drop = T) %>%
  lapply(`[[`, 'samp')
tsamps <- rownames(md)[md$src != 'cell']

cor.rna <- list(pearson = cp, spearman = cs) %>%
  lapply(function(x) {
    lapply(smps, function(y) {
      x[rownames(x) %in% tsamps, y] %>%
        rowMeans() %>%
        data.frame(v = .) %>%
        rownames_to_column('samp')
    }) %>%
      bind_rows(.id = 'cond') %>%
      mutate(samp = sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', samp)) %>%
      group_by(samp, cond) %>%
      summarise(v = mean(v)) %>%
      summarise(v = v[cond == 'WT'] - v[cond == 'KO']) %>%
      deframe()
  })

r <- rowRanges(DNAme) %>%
  as.data.frame() %>%
  dplyr::rename(ref = Composite.Element.REF) %>%
  dplyr::select(chr = seqnames, start, end, ref) %>%
  as.data.table()

d <- as.data.table(assays(data)[[1]])
mat <- mat[, .SD, .SDcols = grep('Cal27|Det|FaDu', names(mat), value = T)]
ok <- lapply(list(mat, d), function(x) {
  which(apply(x, 1, function(y) any(is.finite(y))))
}) %>% Reduce(intersect, .)

cor.5mc <- c('spearman', 'pearson') %>%
  setNames(., .) %>%
  lapply(function(meth) {
    lapply(mat[ok,], function(x) {
      cor(x, d[ok,], use = 'complete.obs', method = meth) %>%
        t() %>%
        data.frame(v = .) %>%
        rownames_to_column('samp')
    }) %>%
      bind_rows(.id = 'ref') %>%
      dplyr::filter(grepl('Cal27|Detroit|FaDu', ref)) %>%
      mutate(samp = sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', samp),
             cond = case_when(grepl('NSD1KO', ref) ~ 'KO',
                              T ~ 'WT')) %>%
      group_by(samp, cond) %>%
      summarise(v = mean(v)) %>%
      summarise(v = v[cond == 'WT'] - v[cond == 'KO']) %>%
      deframe()
  })

ps <- mapply(function(dd, nm) {
  bd <- list('NSD1-' = names(dd) %in% muts) %>%
    {c(., list('NDS1+' = !.[[1]]))} %>%
    lapply(function(x) {
      barcodeData(dd, which(x))$worm %>%
        tibble(y = .) %>%
        mutate(x = 1:n())
    }) %>%
    bind_rows(.id = 'kind')
  
  eclrs <- range(bd$y) %>%
    {setNames(c(seq(.[1], 0, length.out = 51),
                seq(0, .[2], length.out = 51)[2:51]),
              pals::coolwarm(101))}
  xs <- sort(dd) %>%
    names() %>%
    {. %in% muts} %>%
    which()
  
  p1 <- tibble(y = sort(dd),
               nm = nm) %>%
    mutate(x = 1:n()) %>%
    ggplot(aes(x, y)) +
    geom_line() +
    labs(y =  expression(rho('WT') - rho('KO'))) +
    coord_cartesian(clip = 'off', xlim = c(0,525)) +
    theme(panel.background = element_rect(fill = 'grey90'),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey70',
                                            linetype = 'dashed'),
          axis.text.x = element_blank(),
          axis.text = element_text(color = 'black'),
          axis.ticks = element_blank(),
          legend.justification = c(0, 0),
          legend.position = c(0, 0),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.title.x = element_blank(),
          legend.background = element_blank(),
          strip.background = element_rect(fill = 'black'),
          strip.text = element_text(color = 'white'),
          axis.line = element_blank())
  
  
  p2 <- hist(xs, 10, plot = F) %>%
    {tibble(x = .$mids, y = .$counts)} %>%
    ggplot(aes(x, y)) +
    geom_col(width = 50, fill = pals::tableau20(11)[11]) +
    geom_rug(aes(x = x), data = tibble(x = xs),
             inherit.aes = F, sides = 'b',
             size = 0.1,
             length = unit(0.1, "npc"), color = 'black') +
    labs(y = '# of NSD1-') +
    coord_cartesian(xlim = c(0,525), clip = 'off') +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey70',
                                            linetype = 'dashed'),
          axis.text = element_text(color = 'black'),
          axis.ticks.y = element_blank(),
          legend.justification = c(1, 1),
          legend.position = c(1, 1),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.background = element_blank(),
          strip.background = element_rect(fill = 'black'),
          strip.text = element_text(color = 'white'),
          axis.line = element_blank())
  p3 <- ggplot(bd, aes(x, y)) +
    geom_hline(yintercept = 1) +
    geom_line(aes(color = kind), size = 1) +
    scale_color_manual(values = pals::tableau20(3)[c(1,3)]) +
    labs(y = 'Enrichment') +
    scale_x_continuous(breaks = c(1, 500)) +
    coord_cartesian(xlim = c(0,525), clip = 'off') +
    theme(panel.background = element_rect(fill = 'grey90'),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.y = element_line(color = 'grey70',
                                            linetype = 'dashed'),
          axis.text = element_text(color = 'black'),
          axis.ticks.y = element_blank(),
          legend.justification = c(1, 1),
          legend.position = c(1, 1),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.line.x = element_line(color = 'black'),
          axis.title.x = element_blank(),
          legend.background = element_blank(),
          strip.background = element_rect(fill = 'black'),
          strip.text = element_text(color = 'white'),
          plot.caption = element_text(size = 13.1, vjust = 0),
          axis.line = element_blank())
  
  if (nm == 'Gene expression') {
    p1 <- p1 + facet_grid(. ~ nm)
    p3 <- p3 + theme(legend.position = 'none')
  } else {
    p1 <- p1 + facet_grid('Relative similarity' ~ nm) +
      theme(axis.title.y = element_blank())
    p2 <- p2 + facet_grid('Similarity Ranking' ~ .) +
      theme(axis.title.y = element_blank())
    p3 <- p3 + facet_grid('Overrepresentation' ~ .) +
      theme(axis.title.y = element_blank())
  }
  wrap_plots(p1, p2, p3, ncol = 1)
}, list(rna = cor.rna$spearman, dna = cor.5mc$spearman),
c('Gene expression', 'CpG methylation'), SIMPLIFY = F)

wrap_plots(ps, ncol = 2) %>%
  ggsave('figs/6a.pdf', ., height = 4.8, width = 7.8)
