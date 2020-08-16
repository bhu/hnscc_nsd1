library(tidyverse)
library(rtracklayer)
library(msigdbr)
library(fgsea)
library(RobustRankAggreg)


load('misc/cons_ko.rda')
tads <- import.bed('misc/tad.bed')
load('misc/g33.rda')
load('misc/deg.rda')
load('misc/dac.rda')

g33 <- g33[g33$type == 'gene']
dbins <- cons$B
tads <- tads[overlapsAny(tads, g33) & overlapsAny(tads, d.all)]
g33 <- g33[overlapsAny(g33, tads)]
d.all <- d.all[overlapsAny(d.all, tads)]

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


pways <- msigdbr(species = "Homo sapiens") %>%
  mutate(db = paste(gs_cat, gs_subcat, sep = '.') %>%
           sub('\\.$', '', .)) %>%
  split(., .$db) %>%
  lapply(function(db) {
    split(db$gene_symbol, db$gs_name)
  })

enhgs <- mcols(gs[overlapsAny(gs, tads[overlapsAny(tads, dbins)])]) %>%
  as_tibble() %>%
  dplyr::rename(ID = ensembl, gname = gene) %>%
  full_join(res$KO, by = "ID") %>%
  na.omit() %>%
  distinct()
enhs.ok <- enhs[match(enhgs$id, enhs$id)]
k <- findOverlaps(enhs.ok, d.all)
rdat <- tibble(Fold = d.all$Fold,
               pval = d.all$`p-value`,
               fdr = d.all$FDR)
enhgs2 <- cbind(enhgs[k@from,], rdat[k@to,]) %>%
  mutate(sr = ifelse(log2FoldChange > 0, 1, -1),
         sa = ifelse(Fold < 0, 1, -1),
         er = -log10(pvalue) * sr,
         ea = -log10(pval) * sa,
         tmp = 0) %>%
  group_by(gname, id, ID, gene) %>%
  mutate(er = mean(er),
         ea = mean(ea)) %>%
  ungroup() %>%
  distinct(gname, id, ID, gene, .keep_all = T) %>%
  mutate(ID = sub('\\..*', '', ID)) %>%
  mutate(ID = gene) %>%
  group_by(ID) %>%
  summarise(er = mean(er),
            ea = mean(ea)) %>%
  dplyr::filter(sign(er) == sign(ea))


l3 <- list(arrange(enhgs2, er) %>% pull(ID),
           arrange(enhgs2, ea) %>% pull(ID)) %>%
  {inner_join(aggregateRanks(.),
              aggregateRanks(lapply(., rev)),
              by = 'Name')} %>%
  mutate(sc = log10(Score.x) - log10(Score.y)) %>%
  dplyr::select(Name, sc) %>%
  deframe()

fres <- lapply(pways, function(x) {
  fgsea(pathways = x,
        stats = l3,
        minSize = 15,
        maxSize = 500)
})

fresc <- names(pways) %>%
  setNames(., .) %>%
  lapply(function(x) {
    r <- fres[[x]]
    clps <- collapsePathways(r[padj < .05], pways[[x]], l3)
    r[pathway %in% clps$mainPathways][order(-NES),]
  })

pe <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
}

xthres <- which(l3 > 0)[1]
pp <- fresc$H$pathway %>%
  setNames(., .) %>%
  lapply(function(x) {pe(pways$H[[x]],l3)}) %>%
  bind_rows(.id = 'pathway') %>%
  mutate(pathway = sub('^HALLMARK_', '', pathway),
         pathway = fct_inorder(pathway)) %>%
  ggplot(aes(x, y, color = pathway)) +
  geom_hline(yintercept = 0, color = 'black') +
  geom_vline(xintercept = xthres, color = 'black') +
  annotate(geom = 'text', x = xthres, y = Inf, label = 'KO > WT',
           hjust = -0.1, vjust = 1) +
  annotate(geom = 'text', x = xthres, y = Inf, label = 'WT > KO',
           hjust = 1.1, vjust = 1) +
  geom_line() +
  scale_color_manual(values = pals::tableau20(20)[seq(1,20,2)],
                     guide = guide_legend(nrow = 3)) +
  labs(x = "Rank (from WT >> KO to KO >> WT)",
       y = "Enrichment score",
       title = 'MSigDB hallmark enriched gene sets') +
  coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.line.x = element_line(color = 'black'),
        panel.grid.major.y = element_line(color = 'grey50', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.background = element_blank())

ggsave('figs/5a.pdf', pp, height = 3.4, width = 6.3)
