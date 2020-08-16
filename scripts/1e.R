library(purrr)
library(tidyverse)
library(rtracklayer)
library(ggnewscale)
library(patchwork)
library(ggforce)


getWhisks <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, c(0.25, 0.75), na.rm = T)
  data.frame(lower = qs[1], upper = qs[2], middle = median(x, na.rm = T),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))], na.rm = T),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))], na.rm = T)) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols = c("***", "**", "*", "˙", ""))
}

igr <- import.bed('misc/igr.bed')
load('misc/expd.rda')

d <- list.files('data', pattern = '.WGBS.10kb.beta.bw') %>%
  grep('diff|WT|MT|KO', ., invert = T, value = T) %>%
  setNames(sub('.WGBS.*', '', .)) %>%
  mapply(function(x, nm) {
    d <- import.bw(file.path('data', x)) %>%
      {list(Genic = .$score[overlapsAny(., expd[[nm]])],
            IGR = .$score[overlapsAny(., igr)])}
    bx <- lapply(d, getWhisks) %>%
      bind_rows(.id = 'reg') %>%
      mutate(samp = nm)
    vln <- lapply(d, function(y) tibble(y)) %>%
      bind_rows(.id = 'reg') %>%
      drop_na() %>%
      group_by(reg) %>%
      do(density(.$y) %>%
           {data.frame(loc = .$x,
                       dens = .$y)}) %>%
      mutate(samp = nm)
    list(bx = bx, vln = vln)
  }, ., names(.), SIMPLIFY = F)

ttl <- "Actively transcribed (left) vs Intergenic (right)"

bx <- lapply(d, `[[`, 'bx') %>%
  bind_rows() %>%
  mutate(cond = ifelse(grepl('^S|^B', samp), 'MT', 'WT'),
         samp = factor(samp, c('Cal27', 'Det562', 'FaDu',
                               'BICR78', 'SCC4', 'SKN3')),
         grp = interaction(cond, reg),
         ttl = ttl) %>%
  arrange(samp) %>%
  mutate(x = (round((1:n())/2 + .25) - 1) * 10)

shif <- distinct(bx, samp, x) %>% deframe()

vln <- lapply(d, `[[` , 'vln') %>%
  bind_rows() %>%
  mutate(cond = ifelse(grepl('^S|^B', samp), 'MT', 'WT'),
         samp = factor(samp, c('Cal27', 'Det562', 'FaDu',
                               'BICR78', 'SCC4', 'SKN3')),
         grp = interaction(cond, reg),
         ttl = ttl) %>%
  group_by(samp) %>%
  mutate(dens = dens / max(abs(dens)) * 2.25) %>%
  ungroup() %>%
  {rbind(., map_df(mutate(., dens = -dens), rev))} %>%
  mutate(dens = ifelse(reg == "Genic", dens - 2.25,
                       dens + 2.25) + shif[samp])

clrs <- c(WT.Genic = '#9edae5',
          WT.IGR = '#17becf',
          MT.Genic = '#f7b6d2',
          MT.IGR = '#e377c2',
          b = '#00000019',
          a = '#00000000')

lns <- bx %>%
  mutate(x = x + ifelse(reg == 'Genic', -2.25, 2.25))

cns <- lns %>%
  group_by(samp) %>%
  mutate(reg = factor(reg, c('Genic', 'IGR'))) %>%
  arrange(reg) %>%
  mutate(y = middle[1],
         yend = middle[2],
         xend = x[2],
         x = x[1]) %>%
  ungroup() %>%
  select(x, xend, y, yend, samp,
                reg, grp, cond, ttl) %>%
  distinct(samp, .keep_all = T) %>%
  split(., .$cond)
  
rects <- bx %>%
  distinct(samp, x, ttl) %>%
  mutate(xmin = x - 5,
         xmax = x + 5,
         alt = ifelse(1:n() %% 2, 'a', 'b'))

d2 <- by(bx, bx$samp, function(x) {
  x$middle[x$reg == 'IGR'] - x$middle[x$reg == 'Genic']
}, simplify = F) %>%
  unlist() %>%
  data.frame(d = .) %>%
  rownames_to_column("samp") %>%
  left_join(distinct(bx,samp, .keep_all = T)) %>%
  mutate(d = sprintf('%+.2f', d),
         y = 1.05,
         xmin = seq(-2, 48, 10),
         xmax = seq(2, 52, 10),
         ym = 1.13)

p1 <- ggplot(bx, aes(x = x)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                fill = alt), inherit.aes = F, data = rects)  +
  scale_fill_manual(values = clrs) +
  geom_link(aes(x = x, y = y, xend = xend, yend = yend,
                color = stat(index)), size = .4, data = cns$WT) +
  scale_color_gradient(low = clrs['WT.Genic'], high = clrs['WT.IGR']) +
  new_scale_color() +
  geom_link(aes(x = x, y = y, xend = xend, yend = yend,
                color = stat(index)), size = .4, data = cns$MT) +
  scale_color_gradient(low = clrs['MT.Genic'], high = clrs['MT.IGR']) +
  new_scale_color() +
  geom_polygon(aes(x = dens, y = loc, fill = grp,
                   group = interaction(samp, reg)),
               data = vln, inherit.aes = F,
               color = NA, size = 0.25,
               show.legend = F, alpha = 1) +
  geom_point(aes(y = middle),
             size = 1, color = 'white',
             data = lns) +
  geom_segment(aes(x = x, xend = x, y = lower, yend = upper,
                   group = interaction(samp, reg)),
               size = 0.4,
               color = 'white',
               lineend = "round",
               data = lns) +
  geom_segment(aes(x = xmin, y = y, yend = y, xend = xmax),
               data = d2, inherit.aes = F) +
  geom_text(aes(x = x, y = ym, label = d), data = d2,
            inherit.aes = F, size = 3.5) +
  scale_color_manual(values = clrs) +
  labs(y = 'Beta') +
  coord_cartesian(xlim = c(-5, length(shif) * 10 - 5)) +
  scale_y_continuous(limits = c(0,1.15),
                     breaks = c(0,.5,1)) +
  scale_x_continuous(breaks = shif,
                     labels = names(shif),
                     limits = c(min(rects$xmin), max(rects$xmax)),
                     expand = expansion(0)) +
  facet_wrap(. ~ ttl) +
  theme(panel.background = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        axis.line.y.right = element_line(color = "black"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y.right = element_blank(),
        axis.title.y.left = element_text(color = "black"))


dd <- list.files('data', pattern = '.WGBS.10kb.beta.bw') %>%
  grep('diff|WT|MT|KO', ., invert = T, value = T) %>%
  setNames(sub('.WGBS.*', '', .)) %>%
  mapply(function(x, nm) {
    d <- import.bw(file.path('data', x)) %>%
      {list(Genic = .$score[overlapsAny(., expd[[nm]])],
            IGR = .$score[overlapsAny(., igr)])} %>%
      lapply(function(b) tibble(b)) %>%
      bind_rows(.id = 'g') %>%
      mutate(samp = nm)
  }, ., names(.), SIMPLIFY = F) %>%
  bind_rows() %>%
  mutate(cond = factor(ifelse(grepl('^S|^B', samp), 'MT', 'WT')),
         samp = factor(samp))

s <- lm(b ~ g * cond, data = dd) %>%
  summary() %>%
  coefficients() %>%
  {.["gIGR:condWT", "Pr(>|t|)"]} %>%
  signif.num() %>%
  as.character()

dd <- split(dd, dd$g) %>%
  lapply(function(x) {
    split(x, x$cond) %>%
      lapply(`[[`, 'b')
  }) %>%
  unlist(recursive = F)
names(dd) <- sapply(names(dd), function(x) {
  strsplit(x, '\\.')[[1]] %>%
    rev() %>%
    paste(collapse = '.')
})

bx2 <- lapply(dd, getWhisks) %>%
  bind_rows(.id = "samp") %>%
  separate(samp, c('cond', 'reg'), '\\.', remove = F) %>%
  mutate(reg = factor(reg, c('Genic', 'IGR')),
         cond = factor(cond, c('WT', 'MT')),
         x = (as.numeric(reg) - 1) * 10 + ifelse(cond == 'WT', 2, -2),
         ttl = 'Summary')

rcts <- tibble(xmin = c(-4, 5),
                xmax = c(5, 14),
                cond = c('a', 'b'))

pts <- tibble(x = .5,
              y = .5,
              cond = c('WT', 'MT'))

cns2 <- bx2 %>%
  select(reg, cond, y = middle, x) %>%
  group_by(cond) %>%
  summarise(xmin = x[1], xmax = x[2],
            ymin = y[1], ymax = y[2]) %>%
  mutate(m = (ymax - ymin) / (xmax - xmin),
         b = ymin - xmin * m,
         x = 5,
         y = m * x + b)

ds <- split(bx2, bx2$samp) %>%
  lapply(`[[`, 'middle') %>%
  {c(.$MT.Genic - .$WT.Genic, .$MT.IGR - .$WT.IGR)} %>%
  sprintf('%+.2f', .) %>%
  tibble(d = ., y = 1.05, yend = y,
         x = c(-2, 8), xend = c(2, 12)) %>%
  mutate(xm = rowMeans(.[,c('x', 'xend')]),
         ym = 1.13)

p2 <- ggplot(bx2, aes(x = x, group = samp, y = middle)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                ymax = Inf, fill = cond),
            color = NA, inherit.aes = F, data = rcts,
            show.legend = F) +
  geom_link(aes(x = xmin, y = ymin, xend = xmax, yend = ymax,
                color = stat(index)), size = .4, inherit.aes = F,
            data = cns2[cns2$cond == 'WT',]) +
  scale_color_gradient(low = clrs['WT.Genic'], high = clrs['WT.IGR']) +
  new_scale_color() +
  geom_link(aes(x = xmin, y = ymin, xend = xmax, yend = ymax,
                color = stat(index)), size = .4, inherit.aes = F,
            data = cns2[cns2$cond == 'MT',]) +
  scale_color_gradient(low = clrs['MT.Genic'], high = clrs['MT.IGR']) +
  new_scale_color() +
  geom_segment(aes(x = 5, xend = 5, y = cns2$y[1], yend = cns2$y[2]),
                inherit.aes = F, arrow = arrow(length = unit(0.03, "npc"),
                                               ends = 'both'),
               size = .4) +
  geom_text(aes(x = 5, y = mean(cns2$y), label = s),
            vjust = 1.2, angle = 90) +
  geom_boxplot(aes(ymin = ymin, ymax = ymax, fill = samp,
                   lower = lower, upper = upper,
                   notchupper = notchupper,
                   notchlower = notchlower,
                   middle = middle, color = samp),
               notch = T, stat = "identity", width = 2.5,
               show.legend = F) +
  geom_crossbar(aes(ymin = middle, ymax = middle),
                color = "white", width = 1.5, fatten = 0) +
  guides(color = guide_legend(override.aes= list(size = 5))) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values =clrs) +
  scale_x_continuous(expand = expansion(0),
                     breaks = c(0, 10),
                     labels = c('Transcr\'d',
                                'IGR')) +
  scale_y_continuous(limits = c(0,1.15),
                     breaks = c(0,.5,1),
                     name = 'Beta') +
  facet_grid(. ~ ttl) +
  geom_segment(aes(x = x, y = y, yend = yend, xend = xend),
               data = ds, inherit.aes = F) +
  geom_text(aes(x = xm, y = ym, label = d), data = ds,
            inherit.aes = F, size = 3.5) +
  theme(legend.position = "none",
        axis.line.x = element_line(color = "black", size = .5),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        strip.text = element_text(color = 'white'),
        strip.background = element_rect(fill = "black"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.margin=margin(0,0,0,0),
        panel.grid.major.y = element_line(color = 'grey80',
                                          linetype = 'dashed',
                                          size = 0.5),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(color = "black"))
  
wrap_plots(p1, p2, widths = c(3, 1)) %>%
  ggsave('figs/1e.pdf', ., height = 2.1, width = 4.5)
