library(tidyverse)
library(readxl)
library(rtracklayer)
library(msigdbr)
library(fgsea)
library(RobustRankAggreg)
library(forestplot)
library(pals)


load('misc/g33.rda')
load('tcga/deg.rda')
load('tcga/dmr.rda')

asamps <- read_delim('tcga/atac_samps.txt', '\t') %>%
  dplyr::filter(grepl('^HNSC', bam_prefix)) %>%
  mutate(barcode = sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*',
                       '\\1', Case_ID),
         bpfix = gsub('-', '_', bam_prefix)) %>%
  group_by(barcode) %>%
  mutate(nm = gsub('-', '\\.', paste0(barcode, '_', 1:n()))) %>%
  ungroup()

aclus <- read_xlsx('tcga/papillon2017.xlsx', skip = 1) %>%
  dplyr::select(bc = bcr_patient_barcode, grp = methylation_group,
                loc = simple_tissue) %>%
  merge(tibble(bc = unique(asamps$barcode))) %>%
  mutate(bc = gsub('-', '\\.', bc))

rnm <- deframe(asamps[,c('bpfix', 'nm')])
atac <- read_delim('tcga/atac_peaks.txt', '\t') %>%
  makeGRangesFromDataFrame() 

dmat <- readRDS('tcga/atac_cts.rds')
dmat <- dmat[, asamps$bpfix]
names(dmat) <- rnm[names(dmat)]

atac$score <- tibble(nm = names(dmat)) %>%
  mutate(bc = sub('_[12]$', '', nm)) %>%
  dplyr::filter(bc %in% aclus$bc[aclus$grp != 'HPV+' &
                             aclus$loc == 'Larynx']) %>%
  split(., .$bc) %>%
  lapply(function(x) {
    rowMeans(dmat[,x$nm,drop=F])
  }) %>%
  bind_cols() %>%
  scale(scale = F) %>%
  as.data.frame() %>%
  mutate(score = TCGA.CN.A49B - TCGA.KU.A66S) %>%
  pull(score)

meth <- mutate(stats,
               start = pos - 1,
               end = pos + 1,
               score = -estimate) %>%
  dplyr::select(chr, start, end, score) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)



ints <- read_delim('misc/GH_interactions1_doubleElite.bed', '\t',
                   col_names = c("chrom","chromStart", "chromEnd",
                                 "name", "score", "value",
                                 "geneAssociationMethods",
                                 "color", "geneHancerChrom",
                                 "geneHancerStart", "geneHancerEnd",
                                 "geneHancerIdentifier",
                                 "geneHancerStrand", "geneChrom",
                                 "geneStart", "geneEnd",
                                 "geneName", "geneStrand"))

gh <- read_delim('misc/GeneHancer.bed', '\t',
                 col_names = c("chrom", "start", "end", "name", "score",
                               "strand", "thickStart", "thickEnd", "reserved",
                               "evidenceSources", "elementType", "eliteness")) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

igr <- import.bed('ensembl/intergenic.bed')
igr_or_not <- setNames(overlapsAny(gh, igr), gh$name)
enh_or_not <- setNames(gh$elementType == 'Enhancer', gh$name)

enhs <- ints %>%
  dplyr::select(chr = geneHancerChrom,
                start = geneHancerStart,
                end = geneHancerEnd,
                strand = geneHancerStrand,
                id = geneHancerIdentifier) %>%
  mutate(start = start + 1,
         igr = igr_or_not[id],
         enh = enh_or_not[id]) %>%
  distinct(id, .keep_all = T) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)



genes <- ints %>%
  dplyr::select(chr = geneChrom,
                start = geneStart,
                end = geneEnd,
                strand = geneStrand,
                gene = geneName,
                id = geneHancerIdentifier) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

gs <- genes[genes$gene %in% g33$gene_name]
g2id <- setNames(g33$ID, g33$gene_name)
gs$ensembl <- g2id[gs$gene]

others <- genes[!(genes$gene %in% g33$gene_name)]
others.ok <- others[overlapsAny(others, g33)]
others <- others[!overlapsAny(others, g33)]

hits <- findOverlaps(others.ok, g33) %>% as("List")
ids <- extractList(g33$ID, hits) %>% as.list()
init <- T
while(T) {
  uniq <- which(sapply(ids, length) == 1)
  others.ok$ensembl <- sapply(ids, `[`, 1)
  if (init) {
    ok <- others.ok
    init <- F
  } else {
    ok <- c(ok, others.ok)
  }
  others.ok <- others.ok[-uniq]
  if (length(others.ok) < 1) {
    break
  }
  hits <- findOverlaps(others.ok, g33) %>% as("List")
  ids <- extractList(g33$ID, hits) %>% as.list
}
gs <- c(gs, ok)

dat <- na.omit(res.tcga) %>%
  dplyr::filter(name %in% gs$gene) %>%
  mutate(es = sign(stat) * -log10(pvalue)) %>%
  dplyr::select(id = gene, name, adj = log2FoldChangeAdj,
                bm = baseMean,
                fc = log2FoldChange, s = svalue, p = pvalue,
                padj, es) %>%
  group_by(name) %>%
  summarise(id = id[1],
            adj = mean(adj),
            fc = mean(fc),
            bm = mean(bm),
            es = mean(es),
            s = sqrt(sum(s^2)),
            p = sqrt(sum(p^2)),
            padj = sqrt(sum(padj^2))) %>%
  arrange(name)

list(atac = atac, meth = meth) %>%
  mapply(function(x, nm) {
    sc <- findOverlaps(enhs, x) %>%
      as("List") %>%
      extractList(score(x), .) %>%
      mean(na.rm = T) %>%
      setNames(enhs$id)
    dat[[nm]] <<- mcols(gs[gs$gene %in% dat$name]) %>% 
      as_tibble() %>%
      mutate(sc = sc[id]) %>%
      group_by(gene) %>%
      summarise(sc = mean(sc, na.rm = T)) %>%
      arrange(gene) %>%
      pull(sc)
    
    NULL
  }, ., names(.))

gls <- c('fc', 'padj', 'es') %>%
  setNames(., .) %>%
  lapply(function(kind) {
    list(dat$atac, dat$meth, dat[[kind]]) %>%
      lapply(function(x) {
        setNames(x, dat$name) %>%
          na.omit() %>%
          sort() %>%
          names()
      }) %>%
      {inner_join(aggregateRanks(.),
                  aggregateRanks(lapply(., rev)),
                  by = 'Name')} %>%
      mutate(sc = log10(Score.x) - log10(Score.y)) %>%
      dplyr::select(Name, sc) %>%
      deframe()
  })

pways <- msigdbr(species = "Homo sapiens") %>%
  mutate(db = paste(gs_cat, gs_subcat, sep = '.') %>%
           sub('\\.$', '', .)) %>%
  split(., .$db) %>%
  lapply(function(db) {
    split(db$gene_symbol, db$gs_name)
  })

fres <- lapply(gls, function(l) {
  res <- pblapply(pways, function(x) {
    fgsea(pathways = x,
          stats = l,
          minSize = 15,
          maxSize = 500,
          eps = 0.0)
  })
  cons <-  names(res) %>%
    setNames(., .) %>%
    pblapply(function(x) {
      r <- res[[x]]
      clps <- collapsePathways(r[padj < .05], pways[[x]], l)
      r[pathway %in% clps$mainPathways][order(-NES),]
    })
  list(res = res, cons = cons)
})

r <- fres$es$cons$H
p <- pways$H

tbltxt <- r %>%
  top_n(10, -padj) %>%
  arrange(-NES) %>%
  dplyr::select(pathway, padj, NES, size) %>%
  mutate(pathway = sub('^HALLMARK_', '', pathway),
         padj = sprintf('%.2g', padj),
         NES = sprintf('%.2g', NES)) %>%
  rbind(list('Pathway', 'FDR', 'NES', 'Size'), .) %>%
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

d <- r %>%
  top_n(10, -padj) %>%
  arrange(-NES) %>%
  mutate(nm = paste(pathway, padj, NES, sep = '>')) %>%
  split(., .$nm) %>%
  lapply(function(x) {
    dat[dat$name %in% p[[x$pathway]], ] %>%
      mutate(e = es) %>%
      group_by(name) %>%
      summarise_at(c('es', 'atac','meth'),
                   function(y) -mean(y, na.rm = T)) %>%
      mutate(meth = meth * 50) %>%
      pivot_longer(-name, names_to = 'kind', values_to = 'val') %>%
      na.omit()
  }) %>% bind_rows(.id = 'nm') %>%
  separate(nm, c('path', 'padj', 'NES'), '>') %>%
  mutate(padj = as.numeric(padj),
         NES = as.numeric(NES)) %>%
  group_by(kind, path, NES) %>%
  do(getWhisks(.$val)) %>%
  ungroup() %>%
  arrange(-NES) %>%
  split(., .$kind) %>%
  {list(mean = lapply(., `[[`, 'middle'),
        lower = lapply(., `[[`, 'lower'),
        upper = lapply(., `[[`, 'upper'))} %>%
  lapply(function(x) {
    rbind(c(NA, NA), bind_cols(x))
  })

clrs <- tableau20(5)[c(1,3,5)]
pdf('figs/6d.pdf', width = 15, onefile = F)
forestplot(tbltxt,
           fn.ci_norm = c(fpDrawCircleCI, fpDrawCircleCI, fpDrawCircleCI),
           boxsize = 0.2,
           mean = d$mean,
           upper = d$upper,
           lower = d$lower,
           col = fpColors(box = clrs, lines = clrs),
           legend = c('ATAC', 'RNA', 'DNAme'),
           is.summary = c(TRUE, rep(FALSE, nrow(tbltxt)-1)),
           xlab = 'log2(NSD1+ / NSD1-)',
           graphwidth = unit(5, 'cm'),
           xticks = -4:1,
           lwd.ci = 1.5,
           ci.vertices = T)
dev.off()
