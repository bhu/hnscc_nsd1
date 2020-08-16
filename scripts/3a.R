library(tidyverse)
library(rtracklayer)
library(isoband)
library(sf)
library(MASS)
library(lwgeom)
library(ggrepel)
library(hexbin)
library(ggrastr)
library(viridis)
library(pals)
library(patchwork)


load('misc/cons_ko.rda')
gene <- import.bed('ensembl/gene.bed')
igr <- import.bed('ensembl/intergenic.bed')

s <- list.files('data', pattern = '36me2.10kb.') %>%
  grep('diff|WT|MT', ., invert = T, value = T) %>%
  tibble(f = .) %>%
  separate(f, c('samp', 'mark', NA, NA, NA), '\\.', F) %>%
  dplyr::filter(!(samp %in% c('Det562_KO1', 'Det562_KO4', 'Cal27_KO17'))) %>%
  mutate(cond = case_when(grepl('KO', samp) ~ 'KO',
                          grepl('^S|^B', samp) ~ 'MT',
                          T ~ 'WT'),
         line = sub('_.*', '', samp),
         f = file.path('data', f)) %>%
  arrange(cond) %>%
  dplyr::filter(cond != 'MT')

d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw)

r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges()

clus <- case_when(
  overlapsAny(r, cons$A) ~ 'A',
  overlapsAny(r, cons$B) ~ 'B',
  overlapsAny(r, cons$C) ~ 'C',
  T ~ 'NA'
)
clus[clus == 'NA'] <- NA

olap <- tibble(gene = overlapsAny(r, gene),
               igr = overlapsAny(r, igr)) %>%
  mutate(out = case_when(
    gene & !igr ~ 1,
    !gene & igr ~ -1,
    TRUE ~ 0
  )) %>%
  pull(out)

lineclr <- "black"
horz <- F
gradientn1 <- brewer.rdylbu(50)
cramp <- colorRampPalette(c("#000000ff","#ffffff00"), alpha = T)(5)
leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
leg.labs <- c(sprintf('Genic\u25bc'), rep('', 3), '50%',
              rep('', 3), sprintf('Genic\u25b2'))
len <- 9
pal <- brewer.rdylbu(len)
cmat <- seq(0, 255, length.out = len + 1) %>%
  {.[-1]} %>%
  round() %>%
  as.hexmode() %>%
  format(width = 2, upper.case = T) %>%
  lapply(function(x) {
    paste0(pal, x)
  }) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  `colnames<-`(1:len) %>%
  mutate(clr = 1:dplyr::n()) %>%
  reshape2::melt(id.vars = "clr", variable.name = "opa") %>%
  mutate(opa = as.integer(opa))
leg <- ggplot() +
  geom_tile(aes(x = opa, y = clr, fill = value),
            data = cmat) +
  scale_fill_identity() +
  labs(x = "# of bins \u25ba",
       y = "% genic \u25ba") +
  coord_fixed(expand = F) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "white", fill = NA, size = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Arial',
                                  color = "white",
                                  size = 7))

ps <- split(d, s$line) %>%
  lapply(function(x) {
    pdat <- lapply(x[2:1], function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y')) %>%
      mutate(r = olap) %>%
      dplyr::filter(x > quantile(x, .01),
                    x < quantile(x, .99),
                    y > quantile(y, .01),
                    y < quantile(y, .99))
    
    
    hex <- hexbin(pdat$x, pdat$y, xbins = 75, IDs = T)
    pdat$cell <- hex@cID
    hex <- data.frame(hcell2xy(hex),
                      cell = hex@cell,
                      count = hex@count)
    
    t1 <- 'Intergenic vs genic ratio'
    pdat2 <- pdat %>%
      group_by(cell) %>%
      summarise(prop = mean(r, na.rm = T)) %>%
      ungroup %>%
      right_join(hex, by = "cell") %>%
      mutate(logcount = log10(count),
             ttl = t1)
    lim <- data.frame(x = c(min(pdat2$x[pdat2$count > 10]),
                            max(pdat2$x[pdat2$count > 10])) %>%
                        scales::expand_range(mul = .05),
                      y = c(min(pdat2$y[pdat2$count > 10]),
                            max(pdat2$y[pdat2$count > 10])) %>%
                        scales::expand_range(mul = .05))
    yrang <- diff(lim$y)
    tdat <- data.frame(x = rep(median(pdat$x), 2),
                       y = rep(median(pdat$y), 2),
                       c = c(0,1),
                       ttl = t1)
    pow <- 1.25
    pm <- ggplot() +
      geom_point(aes(x = x, y = y, color = c), alpha = 0, data = tdat,
                 show.legend = T) +
      geom_hex(aes(x = x, y = y, fill = prop, alpha = count, color = prop),
               stat = "identity", color = NA, data = pdat2, size = 5) +
      scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                           breaks = leg.brks, labels = leg.labs,
                           limits = c(-1, 1)) +
      scale_color_gradientn(colors = gradientn1, name = "% Genic",
                            breaks = leg.brks, labels = leg.labs,
                            limits = c(-1, 1)) +
      scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F,
                  trans = scales::trans_new("square",
                                            function(x) {
                                              x^(pow)
                                            },
                                            function(x) x^(1/pow))) +
      coord_cartesian(expand = F,
                      xlim = lim$x,
                      ylim = lim$y) +
      facet_grid(.~ttl) +
      labs(x = sprintf('%s WT \u25ba', names(x)[2]),
           y = sprintf('%s KO \u25ba', names(x)[2])) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      theme(panel.background = element_rect(fill = "black"),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            # strip.background = element_rect(fill = "#726866"),
            strip.background = element_rect(fill = 'black'),
            axis.title = element_text(family = "Arial", color = "black"))
    
    pdat <- kde2d(x = pdat$x, y = pdat$y, n = 100)
    brks <- pretty(c(pdat$z), 30)
    nbrks <- length(brks)
    fac <- diff(brks[1:2]) * nrow(d)/sum(pdat$z)
    b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
    bands <- iso_to_sfg(b)
    bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
    if(!all(st_is_valid(bdat))) {
      bdat <- lwgeom::st_make_valid(bdat)
    }
    bnds <- st_bbox(bdat %>% filter(level > 1))
    
    bdat$ttl <- 'H3K36me2 enrichment'
    
    pl <- ggplot() +
      geom_sf(data = bdat, aes(fill = level), color = NA) +
      scale_fill_gradientn(colors = pals::magma(10),
                           name = "# of bins \u25ba") +
      coord_sf(expand = F,
               xlim = lim$x,
               ylim = lim$y) +
      labs(x = sprintf('%s WT \u25ba', names(x)[2]),
           y = sprintf('%s KO \u25ba', names(x)[2])) +
      facet_grid(. ~ ttl) +
      guides(fill = guide_colorbar(barwidth = 3, barheight = 0.5,
                                   title.position = "top", title.hjust = 0.5,
                                   frame.colour = "white", frame.linewidth = 1,
                                   ticks = F)) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      theme(panel.background = element_rect(fill = "black"),
            panel.grid = element_blank(),
            plot.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(color = "black"),
            legend.title = element_text(color = "white", family = "Arial",
                                        size = 7),
            legend.justification = c(0, 1),
            legend.background = element_blank(),
            legend.position = c(0.05, 0.95),
            legend.direction = "horizontal",
            legend.text = element_blank(),
            #axis.ticks = element_line(color = "black"),
            #axis.text = element_text(color = "black"),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            #strip.background = element_rect(fill = "#726866"),
            strip.background = element_rect(fill = 'black'),
            axis.title = element_text(family = "Arial"))
    pl$coordinates$aspect <- function(f) { NULL }
    
    
    d <- lapply(x[2:1], function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y')) %>%
      mutate(clu = clus) %>%
      dplyr::filter(x > quantile(x, .01),
                    x < quantile(x, .99),
                    y > quantile(y, .01),
                    y < quantile(y, .99))
    
    hull <- na.omit(d) %>%
      group_by(clu) %>%
      dplyr::slice(chull(x, y))
    lab <- na.omit(d) %>%
      group_by(clu) %>%
      summarise(x = mean(x),
                y = mean(y))
    cclr <- setNames(tableau20(6)[seq(1,6,2)],
                     c('C', 'B', 'A'))
    
    bdat2 <- bdat
    bdat2$ttl <- 'Density-based clusters'
    d$ttl <- bdat2$ttl[1]
    gradientn <- paste0('#FFFFFF', as.hexmode(round(seq(0,255, length.out = nbrks))))
    pr <- ggplot() +
      geom_sf(data = bdat2, aes(fill = level), color = NA) +
      geom_point_rast(data = d, aes(x = x, y = y, color = clu),
                      alpha = .1, size = .02, raster.dpi = 300,
                      raster.height = 5, raster.width = 5) +
      geom_label_repel(aes(label = clu, x = x, y = y,
                           color = clu), data = lab,
                       show.legend = F,
                       point.padding = NA, box.padding = 0,
                       direction = "x", inherit.aes = F) +
      scale_fill_gradientn(colors = gradientn) +
      scale_color_manual(values = cclr) +
      coord_sf(expand = F,
               xlim = lim$x,
               ylim = lim$y) +
      facet_grid(.~ttl) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      labs(x = sprintf('%s WT \u25ba', names(x)[2]),
           y = sprintf('%s KO \u25ba', names(x)[2])) +
      theme(panel.background = element_rect(fill = "black"),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(color = "black"),
            legend.background = element_blank(),
            #axis.ticks = element_line(color = "black"),
            #axis.text = element_text(color = "black"),
            plot.title = element_blank(),
            strip.text = element_text(color = "white"),
            #strip.background = element_rect(fill = '#726866'),
            strip.background = element_rect(fill = 'black'),
            legend.margin = margin(0.015, 0, 0, 0, unit="npc"),
            axis.title = element_text(family = "Arial"))
    pr$coordinates$aspect <- function(f) { NULL }
    lo <- c(area(t = 1, l = 1, b = 20, r = 20),
            area(t = 2, l = 2, b = 8, r = 8))
    pm_leg <- pm + leg + plot_layout(design = lo)
    wrap_plots(pl, pm_leg, pr, nrow = 1)
  })

ggsave('figs/3a.pdf', ps$Cal27, height = 2.4,
       width = 6.2, device = cairo_pdf)
