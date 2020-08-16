library(tidyverse)
library(data.table)
library(ggformula)


f <- gzcon(file('data/b_cre.mat.gz', open = 'rb'))
s <- readLines(f, 1) %>%
  gsub('"', '', .) %>%
  sub('.*sample_labels:\\[', '', .) %>%
  sub(']}', '', . ) %>%
  strsplit('],sample_boundaries:[0,', fixed = T) %>%
  {tibble(mark = .[[1]][1], j = .[[1]][2])} %>%
  separate_rows(mark, j) %>%
  mutate(j = as.numeric(j),
         i = j - diff(j)[1] + 1)
close(f)

m <- fread('data/b_cre.mat.gz', skip = 1,
           header = F, drop = 1:6) %>%
  as.matrix()

d <- by(s, s$mark, function(x) {
    m[,x$i:x$j] %>%
      colMeans(na.rm = T) %>%
      tibble(y = .) %>%
      mutate(mark = x$mark[1],
             x = 1:n())
  }, simplify = F) %>%
  do.call(rbind, .) %>%
  mutate(mark = factor(mark, s$mark))

rects <- distinct(d, mark) %>% arrange(mark) %>%
  mutate(alt = rep(c('a', 'b'), 2))

p <- d %>%
  mutate(grp = 'KO - WT') %>%
  ggplot(aes(x, y)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                fill = alt), data = rects, inherit.aes = F) +
  geom_spline(nknots = 100) +
  scale_x_continuous(breaks = c(100, 500, 900),
                     limits = c(1, 1000),
                     labels = c('-40kb', 'Cluster B CRE', '+40kb'),
                     expand = expansion(0)) +
  scale_fill_manual(values = c(a = '#00000000', b = '#00000019')) +
  labs(y = "KO - WT") +
  facet_wrap(. ~ mark, scales = "free_y", nrow = 1) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black"),
        legend.position = "none",
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white"))

ggsave('figs/3c.pdf', p, height = 1.4, width = 9)
