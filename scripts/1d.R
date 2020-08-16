library(data.table)
library(tidyverse)
library(ggsignif)
library(scales)


getWhisks <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  qs <- quantile(x, c(0.25, 0.75))
  data.frame(lower = qs[1], upper = qs[2], middle = median(x),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))]),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))])) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}

signif.num <- function(x) {
  symnum(x, corr = F, na = F, legend = F,
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "**", "*", "˙", ""))
}


d <- list.files('data', 'WTvsMT') %>%
  grep('mat.gz', ., value = T) %>%
  lapply(function(x) {
    f <- gzcon(file(file.path('data', x), open = 'rb'))
    s <- readLines(f, 1) %>%
      gsub('"', '', .) %>%
      sub('.*sample_labels:\\[', '', .) %>%
      sub(']}', '', . ) %>%
      strsplit('],sample_boundaries:[0,', fixed = T) %>%
      {tibble(samp = .[[1]][1], j = .[[1]][2])} %>%
      separate_rows(samp, j) %>%
      mutate(j = as.numeric(j),
             i = j - diff(j)[1] + 1,
             cond = case_when(grepl('^S|^B', samp) ~ 'MT',
                              T ~ 'WT'))
    close(f)
    
    m <- fread(file.path('data', x), skip = 1,
               header = F, drop = 1:6) %>%
      as.matrix()
    d <- by(s, s$cond, function(y) {
      by(y, y$samp, function(z) {
        m[,z$i:z$j] %>%
          {apply(.[,26:35], 1, median, na.rm = T) /
              apply(.[,c(6:15, 46:55)], 1, median, na.rm = T)} %>%
          log2()
      }, simplify = F) %>% unlist()
    }, simplify = F)
    
    bx <- lapply(d, getWhisks) %>%
      bind_rows(.id = "samp") %>%
      mutate(mark = strsplit(x, '\\.')[[1]][2],
             samp = factor(samp, levels = c('WT', 'MT')),
             x = (as.numeric(samp) - 1) * 10 - 2)
    vln <- lapply(d, function(x) {
      na.omit(x) %>%
        density() %>%
        {data.frame(loc = .$x,
                    dens = .$y)} %>%
        mutate(dens = dens / max(dens) * 4)
    }) %>% bind_rows(.id = "samp") %>%
      mutate(samp = factor(samp, c('WT', 'MT')),
             dens = dens + (as.numeric(samp) - 1) * 10,
             mark = bx$mark[1])
    sig <- wilcox.test(d$WT, d$MT)$p.value
    list(bx = bx, vln = vln, sig = sig)
  })

s <- lapply(c(bx = 'bx', vln = 'vln'), function(x) {
  lapply(d, `[[`, x) %>%
    bind_rows() %>%
    dplyr::rename(cond = samp) %>%
    mutate(mark = case_when(
             grepl('^K', mark) ~ paste0('H3', mark),
             T ~ 'DNAme'),
           mark = factor(mark, c('H3K36me2', 'H3K27me3', 'DNAme')),
           grp = interaction(mark, cond))
})

sigs <- lapply(d, `[[`, 'sig') %>%
  unlist() %>%
  signif.num()

clrs <- c('H3K27me3.WT' = '#9467bd',
          'H3K27me3.MT' = '#ff7f0e',
          'H3K36me2.WT' = '#1f77b4',
          'H3K36me2.MT' = '#d62728',
          'DNAme.WT' = '#17becf',
          'DNAme.MT' = '#e377c2')

rects <- s$bx %>%
  distinct(cond, mark, grp) %>%
  group_by(mark) %>%
  arrange(cond) %>%
  mutate(x = (1:n() - 1) * 10,
         xmin = x - 5,
         xmax = x + 5) %>%
  ungroup()

brks <- distinct(rects, x, cond)
ylims <- by(s$bx, s$bx$mark, function(x)
  c(min(x$ymin), max(x$ymax)), simplify = F)
yrang <- lapply(ylims, expand_range, c(0.01, 0.1))
vln2 <- by(s$vln, s$vln$mark, function(x) {
  x[x$loc %between% yrang[[as.character(x$mark[1])]],]
}, simplify = F) %>% do.call(rbind, .)

anns <- rects %>%
  distinct(mark) %>%
  mutate(ann = as.character(sigs),
         xmin = 0,
         xmax = 10) %>%
  rowwise() %>%
  mutate(y = ylims[[as.character(mark)]][2] +
           diff(ylims[[as.character(mark)]]) * .05) %>%
  ungroup()

p <- ggplot(s$bx, aes(x = x, group = cond, y = middle)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                ymax = Inf, fill = grp), alpha = .1,
            color = NA, inherit.aes = F, data = rects) +
  geom_polygon(aes(x = dens, y = loc, fill = grp),
               data = vln2, inherit.aes = F,
               color = "white",
               show.legend = F) +
  geom_boxplot(aes(ymin = ymin, ymax = ymax, fill = grp,
                   lower = lower, upper = upper,
                   notchupper = notchupper,
                   notchlower = notchlower,
                   middle = middle, color = grp),
               notch = T, stat = "identity", width = 1.5,
               position = position_dodge(5),
               show.legend = F) +
  geom_crossbar(aes(ymin = middle, ymax = middle),
                color = "white", width = 0.7, fatten = 0,
                position = position_dodge(5)) +
  geom_signif(aes(xmin = xmin, xmax = xmax,
                  annotations = ann, y_position = y),
              vjust = 0.5,
              data = anns, manual = T, inherit.aes = F) +
  facet_wrap(. ~ mark, scales = "free_y") +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_x_continuous(expand = expansion(0),
                     limits = c(min(rects$xmin), max(rects$xmax)),
                     breaks = brks$x,
                     labels = brks$cond) +
  scale_y_continuous(name = 'log2(intergenic/genic)') +
  theme(legend.position = "none",
        axis.line.x = element_line(color = "black", size = .5),
        legend.background = element_blank(),
        legend.title = element_blank(),
        strip.text = element_text(color = 'white'),
        strip.background = element_rect(fill = "black"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major.y = element_line(color = 'grey80',
                                          linetype = 'dashed',
                                          size = 0.5),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave('figs/1d.pdf', p, height = 2.1, width = 4.5)
