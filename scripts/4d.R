library(tidyverse)
library(rtracklayer)
library(ggrastr)
library(patchwork)
library(furrr)
plan(multiprocess(workers = 10))


load('misc/cons_ko.rda')
tads <- import.bed('misc/tad.bed')
load('misc/g33.rda')
load('misc/deg.rda')
load('misc/dac.rda')

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols = c("***", "**", "*", "˙", ""))
}

getWhisks <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, c(0.25, 0.75), na.rm = T)
  data.frame(lower = qs[1], upper = qs[2], middle = median(x, na.rm = T),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))], na.rm = T),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))], na.rm = T)) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}

scientific_10 <- function(x) {
  xout <- gsub("1e", "10^{", format(x),fixed=TRUE)
  xout <- gsub("{-0", "{-", xout,fixed=TRUE)
  xout <- gsub("{+", "{", xout,fixed=TRUE)
  xout <- gsub("{0", "{", xout,fixed=TRUE)
  xout <- paste(xout,"}",sep="")
  parse(text=xout)
}

scale_x_log10nice <- function(name=NULL,omag=seq(-10,20),...) {
  breaks10 <- 10^omag
  scale_x_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

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

res <- res$KO
res$dtad <- res$ID %in%
  g33[overlapsAny(g33, tads[overlapsAny(tads, dbins)])]$ID
res$tany <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.all[
    d.all$FDR < .05])]$id]$ensembl %>% unique())
res$tpos <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.all[
    d.all$FDR < .05 & d.all$Fold > 0])]$id]$ensembl %>% unique())
res$tneg <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.all[
    d.all$FDR < .05 & d.all$Fold < 0])]$id]$ensembl %>% unique())
res$targ <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, dbins)]$id]$ensembl %>% unique())
res$dpos <- !is.na(res$svalue) & res$svalue < .05 & res$log2FoldChangeAdj < 0
res$dneg <- !is.na(res$svalue) & res$svalue < .05 & res$log2FoldChangeAdj > 0

res <- res[res$ID %in% g33$ID,]
gg <- g33[match(res$ID, g33$ID)]

res$str <- distanceToNearest(gg, dbins)@elementMetadata$distance %>%
  cut(., breaks = quantile(., probs = seq(0, 1, by = .2)),
      include.lowest = T) %>%
  as.factor()

res$dist <- distanceToNearest(gg, dbins)@elementMetadata$distance

o <- res
o <- o[o$dist > 0,]
ht <- o$dist[o$dneg] %>%
  log10() %>%
  hist(10, plot = F)
o$str2 <- cut(o$dist, breaks = 10^ht$breaks,
              include.lowest = T)
o$int <- as.numeric(o$str2)
probs <- lapply(seq_along(ht$counts), function(x) {
  rep(x, ht$counts[x])
}) %>% unlist()
samp <- function(dat, idxs) {
  lapply(idxs, function(i) {
    dat[[as.character(i)]] %>%
      {.[sample(nrow(.), 1),]}
  }) %>% bind_rows()
}
dat <- split(o, o$int)
nums <- future_map(1:1e4, function(x) {
  samp(dat, probs)$dtad %>%
    sum()
}, .progress = T) %>% unlist()
  
num <- sum(o$dtad[o$dneg])
p <- sum(nums >= num) / length(nums)
ptxt <- tibble(x = -Inf, y = Inf,
               p = sprintf('p = %.3f', p))

p <- ggplot(tibble(x = nums, ttl = 'Permutation test on \'TAD\''), aes(x)) +
  geom_histogram(binwidth = 1, fill = 'dodgerblue4') +
  geom_vline(xintercept = num, color = 'firebrick3') +
  geom_label(aes(x, y, label = p), data = ptxt, vjust = 1.5, hjust = -0.1) +
  labs(x = sprintf("# genes in TAD with \u2193%s", 'H3K36me2 bin'),
       y = "Frequency") +
  facet_grid(ttl ~ .) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(.025, .90),
        legend.justification = c(0, 1),
        legend.title = element_text(family = 'Arial', size = 9),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title = element_text(family = 'Arial'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        axis.line = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        #axis.line.y = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.ticks = element_line(color = 'black'))
ggsave('figs/4d.pdf', p, height = 2.1, width = 3.6,
       device = cairo_pdf)

