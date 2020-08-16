library(tidyverse)
library(rtracklayer)


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

g33 <- g33[g33$type == 'gene']

dbins <- cons$B

tads <- tads[overlapsAny(tads, g33) & overlapsAny(tads, d.o)]
g33 <- g33[overlapsAny(g33, tads)]
d.o <- d.o[overlapsAny(d.o, tads)]


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
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.o[
    d.o$FDR < .05])]$id]$ensembl %>% unique())
res$tpos <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.o[
    d.o$FDR < .05 & d.o$Fold > 0])]$id]$ensembl %>% unique())
res$tneg <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, d.o[
    d.o$FDR < .05 & d.o$Fold < 0])]$id]$ensembl %>% unique())
res$targ <- res$ID %in%
  (gs[gs$id %in% enhs[overlapsAny(enhs, dbins)]$id]$ensembl %>% unique())
res$dpos <- !is.na(res$svalue) & res$svalue < .05 & res$log2FoldChangeAdj < 0
res$dneg <- !is.na(res$svalue) & res$svalue < .05 & res$log2FoldChangeAdj > 0


deg <- res[!is.na(res$svalue) & res$svalue < .05,]
fc <- deg$log2FoldChange

dd <- list(All = fc,
           `Cluster B TAD` = fc[deg$dtad],
          `Cluster B CRE` = fc[deg$targ]) %>%
          #`K27ac (WT != KO)` = fc[deg$tany],
          #`K27ac (WT > KO)` = fc[deg$tneg],
          #`K27ac (KO > WT)` = fc[deg$tpos]) %>%
  rev() %>%
  lapply(function(x) -x)

samps <- names(dd)
bx <- lapply(dd, getWhisks) %>%
  bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, levels = samps),
         x = (as.numeric(samp) - 1) * 10 - 3,
         ttl = 'Differentially expressed genes')
vln <- lapply(dd, function(x) {
  density(x) %>%
    {data.frame(loc = .$x,
                dens = .$y)} %>%
    mutate(dens = dens / max(dens) * 4)
}) %>% bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, samps),
         dens = dens + (as.numeric(samp) - 1) * 10)
pts <- lapply(dd, function(x){tibble(y = x)}) %>%
  bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, levels = samps),
         x = (as.numeric(samp) - 1) * 10 - 1.1)
xttl <- 'LFC (KO - PA)'
clrs <- rev(pals::tableau20(20)[seq(1,6,2)]) %>%
  setNames(levels(bx$samp)) %>%
  c(c('b' = '#00000000', 'a' = '#00000019'))
rects <- bx %>%
  distinct(samp, x, ttl) %>%
  mutate(xmin = x - 2,
         xmax = xmin + 10,
         ymin = -Inf,
         ymax = Inf,
         alt = ifelse(as.numeric(samp) %% 2 == 0, 'a', 'b'))
anns <- list(pos = bx %>%
               distinct(samp) %>%
               mutate(num = sprintf('%d', sapply(dd, function(x) {sum(x>0)})[samp]),
                      x = (as.numeric(samp) - 1) * 10,
                      y = Inf,
                      dir = 'pos'),
             neg = bx %>%
               distinct(samp) %>%
               mutate(num = sprintf('%d', sapply(dd, function(x) {sum(x<0)})[samp]),
                      x = (as.numeric(samp) - 1) * 10,
                      y = -Inf,
                      dir = 'neg'
               ))
plt <- ggplot(bx, aes(x = x, group = samp, y = middle)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin,
                ymax = ymax, fill = alt),
            data = rects, inherit.aes = F) +
  geom_text(aes(x = x, y = y, label = num), hjust = 1, vjust = 0,
            data = anns$pos, inherit.aes = F, nudge_x = 0.5) +
  geom_text(aes(x = x, y = y, label = num), hjust = 0, vjust = 0,
            data = anns$neg, inherit.aes = F, nudge_x = 0.5) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_jitter(aes(x = x, y = y, color = samp), data = pts,
              inherit.aes = F, width = 0.5, size = .05, alpha = 0.5) +
  geom_polygon(aes(x = dens, y = loc, fill = samp),
               data = vln, inherit.aes = F,
               color = NA,
               show.legend = F, alpha = 1) +
  geom_boxplot(aes(ymin = ymin, ymax = ymax, fill = samp,
                   lower = lower, upper = upper,
                   notchupper = notchupper,
                   notchlower = notchlower,
                   middle = middle, color = samp),
               notch = F, stat = "identity", width = 1.5,
               position = position_dodge(5),
               show.legend = F) +
  geom_crossbar(aes(ymin = middle, ymax = middle),
                color = "white", width = 0.7, fatten = 0,
                position = position_dodge(5)) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs) +
  scale_x_continuous(breaks = seq(0,length(levels(bx$samp)) - 1) * 10,
                     labels = levels(bx$samp),
                     expand = expansion(0)) +
  ylab(xttl) +
  facet_grid(. ~ ttl) +
  coord_flip() +
  #coord_cartesian(ylim = c(-3, 3)) +
  theme(panel.background = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.grid.major.y = element_line(color = "grey30", linetype = "solid"),
        panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
        axis.line.y.right = element_line(color = "black"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black"))

ggsave('figs/4a.pdf', plt, height = 2.3, width = 3)
