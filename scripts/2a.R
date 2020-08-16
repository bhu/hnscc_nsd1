library(tidyverse)
library(reshape2)
library(readxl)
library(cowplot)
library(patchwork)
library(ggsignif)


clrs <- c('H3K27me3.WT' = '#9467bd',
          'H3K27me3.KO' = '#cd8162',
          'H3K36me2.WT' = '#1f77b4',
          'H3K36me2.KO' = '#CD9B1D')

shps <-  c('Cal27' = 15, 'Det562' = 16, 'FaDu' = 17,
            'Cal27_KO1' = 15, 'Det562_KO2' = 16, 'FaDu_KO1' = 17)

d <- read_excel('data/ms.xlsx') %>%
  pivot_longer(-c(sample, condition),
               names_to = 'ptm', values_to = 'y') %>%
  filter(sample %in% names(shps)) %>%
  mutate(condition = factor(condition, c('WT', 'KO', 'MT')),
         grp = interaction(ptm, condition),
         ptm = factor(ptm, c('H3K36me2', 'H3K27me3')),
         shape = shps[sample]) %>%
  na.omit()

bg <- d %>%
  distinct(condition, ptm) %>%
  mutate(xmin = dense_rank(as.numeric(condition)) - 0.5,
         xmax = xmin + 1,
         grp = interaction(ptm, condition))

ext <- .1
anns <- d %>%
  group_by(ptm) %>%
  mutate(rang = diff(range(y))) %>%
  group_by(ptm, condition) %>%
  summarise(mu = mean(y),
            ma = max(y),
            rang = rang[1]) %>%
  arrange(condition) %>%
  summarise(y = max(ma[1], ma[2]) + rang[1] * ext,
            ann = sprintf("%+.3g", mu[2] - mu[1])) %>%
  mutate(xmin = 1, xmax = 2)

lns <- d %>%
  mutate(line = sub('_.*', '', sample),
         x = c('a', 'b')[as.numeric(condition)]) %>%
  distinct(x, line, y, ptm) %>%
  pivot_wider(names_from = x, values_from = y)

p <- ggplot(d, aes(x = condition, y = y, color = grp)) +
  scale_x_discrete(expand = expansion(add = .5)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                ymax = Inf, fill = grp), alpha = .1,
            color = NA, inherit.aes = F, data = bg) +
  geom_segment(aes(x = 1, xend = 2, y = a, yend = b), data = lns,
               inherit.aes = F,
               color = "black", alpha = 0.3) +
  geom_point(aes(shape = shape), size = 1.5, stroke = 1.5) +
  scale_shape_identity() +
  stat_summary(aes(color = grp), geom = "point",
               fun = mean, shape = 23, size = 2, stroke = 1) +
  facet_wrap(. ~ ptm, nrow = 1, scales = "free") +
  geom_signif(aes(xmin = xmin, xmax = xmax,
                  annotations = ann, y_position = y),
              vjust = -.8, textsize = 3.5, 
              data = anns, manual = T, inherit.aes = F) +
  labs(y = "% of peptides with modification") +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_y_continuous(expand = expansion(c(.05,.15))) +
  guides(shape = guide_legend(nrow = 2, byrow = T)) +
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

lgd <- d %>%
  filter(condition == 'WT') %>%
  distinct(shape, sample) %>%
  mutate(sample = factor(sample, names(shps))) %>%
  ggplot(aes(x = sample, y = shape, shape = sample)) +
  geom_point() +
  scale_shape_manual(values = shps, name = NULL) +
  guides(shape = guide_legend(byrow = T)) +
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())

wrap_plots(p, get_legend(lgd), nrow = 2, heights = c(5,1)) %>%
  ggsave('figs/2a.pdf', .,
         width =  3 , height = 3.7)
