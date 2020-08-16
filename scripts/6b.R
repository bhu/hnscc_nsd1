library(tidyverse)
library(GenomicRanges)
library(pals)
library(ggsignif)


load('tcga/dmr.rda')
load('misc/cons_ko.rda')

dbins <- cons$B
dmrss <- dmrs[dmrs$p.adjust < .05,]
rs <- makeGRangesFromDataFrame(dmrss)

dd <- list('All' = -dmrss$estimate,
           'Cluster B' = -dmrss$estimate[overlapsAny(rs, dbins)])

cgs <- dplyr::select(stats, chr, start = pos, end = pos) %>%
  makeGRangesFromDataFrame()

dd <- list('All' = stats$estimate,
           'Cluster B' = stats$estimate[overlapsAny(cgs, dbins)])



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

samps <- names(dd)
dta <- 7.5
bx <- lapply(dd, getWhisks) %>%
  bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, levels = samps),
         x = (as.numeric(samp) - 1) * dta - 1,
         ttl = 'mCG in TCGA vs')
vln <- lapply(dd, function(x) {
  density(x) %>%
    {data.frame(loc = .$x,
                dens = .$y)} %>%
    mutate(dens = dens / max(dens) * 2)
}) %>% bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, samps),
         dens = dens + (as.numeric(samp) - 1) * dta + 1)
pts <- lapply(dd, function(x){tibble(y = x)}) %>%
  bind_rows(.id = "samp") %>%
  mutate(samp = factor(samp, levels = samps),
         x = (as.numeric(samp) - 1) * dta)
rects <- bx %>%
  distinct(x, ttl, samp) %>%
  mutate(xmax = x +1 + dta / 2,
         xmin = xmax - dta,
         ymin = -Inf,
         ymax = Inf)
xttl <- 'MT - WT'

anns <- list(pos = bx %>%
               distinct(samp) %>%
               mutate(num = sprintf('%d', sapply(dd, function(x) {sum(x>0)})[samp]),
                      x = (as.numeric(samp) - 1) * dta,
                      y = Inf,
                      dir = 'pos'),
             neg = bx %>%
               distinct(samp) %>%
               mutate(num = sprintf('%d', sapply(dd, function(x) {sum(x<0)})[samp]),
                      x = (as.numeric(samp) - 1) * dta,
                      y = -Inf,
                      dir = 'neg'
               ))
plt <- ggplot(bx, aes(x = x, group = samp, y = middle)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin,
                ymax = ymax, fill = samp), alpha = 0.1,
            data = rects, inherit.aes = F) +
  geom_hline(yintercept = 0) +
  geom_polygon(aes(x = dens, y = loc, fill = samp),
               data = vln, inherit.aes = F,
               color = NA,
               show.legend = F) +
  geom_signif(y_position=.15, xmin=0, xmax=dta,
              annotation="", tip_length = .02) +
  annotate("text", x = dta / 2, y = .165, label = "***", angle = 90) +
  geom_boxplot(aes(ymin = ymin, ymax = ymax, fill = samp,
                   lower = lower, upper = upper,
                   notchupper = notchupper,
                   notchlower = notchlower,
                   middle = middle, color = samp),
               notch = T, stat = "identity", width = 1.5,
               position = position_dodge(5),
               show.legend = F) +
  geom_crossbar(aes(ymin = middle, ymax = middle),
                color = "white", width = 0.7, fatten = 0,
                position = position_dodge(5)) +
  scale_fill_manual(values = pals::tableau20(20)[seq(1,8,2)]) +
  scale_color_manual(values = pals::tableau20(20)[seq(1,8,2)]) +
  scale_x_continuous(breaks = seq(0,length(levels(bx$samp)) - 1) * dta,
                     labels = levels(bx$samp), expand = expansion(0)) +
  ylab(xttl) +
  facet_grid(. ~ ttl,
             labeller = label_bquote(cols = Delta ~ .('mCG between TCGA samples in NSD1-KO affected regions'))) +
  coord_flip(ylim = c(-.25,.15)) +
  theme(panel.background = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black"))

ggsave('figs/6b.pdf', plt, height = 2.15, width = 4.35)

