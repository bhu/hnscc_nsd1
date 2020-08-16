library(tidyverse)
library(rtracklayer)


igr <- import.bed('misc/igr.bed')
load('misc/expd.rda')
load('misc/cons_ko.rda')
load('misc/g33.rda')
pcg <- g33[g33$gene_type == 'protein_coding' & g33$type == 'gene']

comm <- function(a, b) { granges(a)[overlapsAny(a, b)] }
getWhisks <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, c(0.25, 0.75), na.rm = T)
  data.frame(lower = qs[1], upper = qs[2], middle = median(x, na.rm = T),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))], na.rm = T),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))], na.rm = T)) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}

s <- list.files('data', pattern = '.10kb.beta') %>%
  grep('diff|WT|MT', ., invert = T, value = T) %>%
  tibble(f = .) %>%
  separate(f, c('samp', 'mark', NA, NA, NA), '\\.', F) %>%
  dplyr::filter(!(samp %in% c('Det562_KO1', 'Det562_KO4', 'Cal27_KO17'))) %>%
  mutate(cond = case_when(grepl('KO', samp) ~ 'KO',
                          grepl('^S|^B', samp) ~ 'MT',
                          T ~ 'WT'),
         line = sub('_.*', '', samp),
         f = file.path('data', f),
         mark = sub('WGBS', 'DNAme', mark)) %>%
  arrange(cond) %>%
  dplyr::filter(cond != 'MT')

d <- by(s, s$line, function(x) {
  d <- deframe(x[,c('samp', 'f')]) %>%
    lapply(import.bw)
  r <- Reduce(comm, d) 
  r$score <- lapply(d, function(y) {
    findOverlaps(r, y) %>%
      to() %>%
      {y[.]} %>%
      score()
  }) %>%
    {.[[1]] - .[[2]]} * 100
  d <- comm(expd[[x$samp[1]]], expd[[x$samp[2]]]) %>%
    {r$score[overlapsAny(r, .)]} %>%
    {list(act = .,
          Genic = r$score[overlapsAny(r, pcg)],
          Intergenic = r$score[overlapsAny(r, igr)],
          k36 = r$score[overlapsAny(r, cons$B)])}
  samps <- names(d)
  bx <- lapply(d, getWhisks) %>%
    bind_rows(.id = "samp") %>%
    mutate(samp = factor(samp, levels = samps),
           ttl = x$samp[2],
           x = (as.numeric(samp) - 1) * 10 - 3)
  vln <- lapply(d, function(x) {
    density(x) %>%
      {data.frame(loc = .$x,
                  dens = .$y)} %>%
      mutate(dens = dens / max(dens) * 4)
  }) %>% bind_rows(.id = "samp") %>%
    mutate(samp = factor(samp, samps),
           dens = dens + (as.numeric(samp) - 1) * 10,
           ttl = x$samp[2])
  list(bx = bx, vln = vln)
})

bx <- lapply(d, `[[`, 'bx') %>% bind_rows() %>% mutate(m = round(middle, 1))
vln <- lapply(d, `[[`, 'vln') %>% bind_rows()
xlabs <- bx %>% distinct(x, samp) %>%
  mutate(x = x - min(x),
         samp = sub('k36', sprintf('H3K36me2%s', '\u2193'), samp) %>%
           sub('act', 'Transcr\'d', .))
  

rects <- bx %>%
  distinct(x, samp, ttl) %>%
  mutate(xmin = x - min(x) - 5,
         xmax = xmin + 10,
         alt = rep(c('a', 'b'), n()/2))
clrs <- viridis::viridis(4, begin = 0.85, end = 0) %>%
  setNames(levels(bx$samp)) %>%
  {c(., a = '#00000000', b = '#00000019')}

p <- ggplot(bx, aes(x = x, group = samp, y = middle)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                fill = alt), color = NA, data = rects, inherit.aes = F) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(aes(x = x + 3, y = Inf, label = m), vjust = 3,
            size = 3) +
  geom_polygon(aes(x = dens, y = loc, fill = samp),
               data = vln, inherit.aes = F,
               color = NA,
               show.legend = F, alpha = 1) +
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
  facet_grid(. ~ ttl) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs) +
  labs(x = "Region", y = sprintf('\u0394%s\n(KO - WT)', 'DNAme%')) +
  scale_x_continuous(breaks = xlabs$x,
                     labels = xlabs$samp,
                     limits = c(min(xlabs$x) - 5, max(xlabs$x) + 5),
                     expand = c(0,0)) +
  scale_y_continuous(breaks = c(-20, 0, 20),
                     labels = c(-20, 0, 20),
                     limits = c(-25, 25)) +
  theme(plot.background = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        axis.text.y = element_text(color = 'black'),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, family='Arial'),
        axis.title.x = element_blank(),
        axis.title = element_text(family = 'Arial'),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"))
ggsave('figs/2d.pdf', p, height = 2.2, width = 4.6,
       device = cairo_pdf, bg = "transparent")
