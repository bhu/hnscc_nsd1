library(tidyverse)
library(rtracklayer)
library(isoband)
library(sf)
library(MASS)
library(lwgeom)
library(viridis)
library(patchwork)


s <- list.files('data', pattern = '.10kb.') %>%
  grep('diff|WT|MT', ., invert = T, value = T) %>%
  tibble(f = .) %>%
  separate(f, c('samp', 'mark', NA, NA, NA), '\\.', F) %>%
  dplyr::filter(!(samp %in% c('Det562_KO1', 'Det562_KO4', 'Cal27_KO17')) &
                  mark != 'H3K27ac') %>%
  mutate(cond = case_when(grepl('KO', samp) ~ 'KO',
                          grepl('^S|^B', samp) ~ 'MT',
                          T ~ 'WT'),
         line = sub('_.*', '', samp),
         f = file.path('data', f),
         mark = sub('WGBS', 'DNAme', mark)) %>%
  arrange(cond)

ps <- by(s, s$mark, function(x) {
    d <- deframe(x[,c('samp', 'f')]) %>%
      lapply(import.bw)
    if (x$mark[1] != 'DNAme') {
      r <- lapply(d, function(y) y[y$score != 0]) %>%
        Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
        granges()
    } else {
      r <- Reduce(function(a, b) a[overlapsAny(a, b)], d) %>%
        granges()
    }
    d1 <- list(MT = x$cond == 'MT', WT = x$cond == 'WT') %>%
      lapply(function(y) {
        lapply(d[y], function(z) {
          findOverlaps(r, z) %>%
            to() %>%
            {z[.]} %>%
            score()
        }) %>%
          bind_cols() %>%
          rowMeans(na.rm = T)
      }) %>%
      {.[[1]] - .[[2]]}
    d2 <- split(d[x$cond != 'MT'], x$line[x$cond != 'MT']) %>%
      lapply(function(y) {
        lapply(y, function(z) {
          findOverlaps(r, z) %>%
            to() %>%
            {z[.]} %>%
            score()
        }) %>%
          {.[[1]] - .[[2]]}
      }) %>%
      bind_cols() %>%
      rowMeans(na.rm = T)
    
    pdat <- kde2d(x = d1, y = d2, n = 100)
    brks <- pretty(c(pdat$z), 30)
    nbrks <- length(brks)
    b <- isobands(x = pdat$x, y = pdat$y, z = t(pdat$z), brks[1:(nbrks - 1)], brks[2:nbrks])
    bands <- iso_to_sfg(b)
    bdat <- st_sf(level = 1:length(bands), geometry = st_sfc(bands))
    if(!all(st_is_valid(bdat))) {
      bdat <- lwgeom::st_make_valid(bdat)
    }
    bnds <- st_bbox(bdat %>% dplyr::filter(level > 1))
    
    xlims <- bnds[c(1,3)]
    ylims <- bnds[c(2,4)]
    xrang <- diff(xlims)
    yrang <- diff(ylims)
    xlims <- c(xlims[1] - xrang * 0.05, xlims[2] + xrang * 0.05)
    ylims <- c(ylims[1] - yrang * 0.05, ylims[2] + yrang * 0.05)
    xrang <- diff(xlims)
    yrang <- diff(ylims)
    gradientn <- magma(10)
    gradientn[1] <- sub('..$', '00', gradientn[1])
    
    bdat$ttl <- x$mark[1]
    cc <- cor(d1, d2, method = 'spearman', use = 'complete.obs') %>%
      sprintf("rho == %0.2f", .)
    
    p <- ggplot() +
      geom_sf(data = bdat, aes(fill = level), color = NA) +
      scale_fill_gradientn(colors = gradientn) +
      annotate('text', x = -Inf, y = Inf, hjust = 0, vjust = 1,
               label = cc, color = 'white', parse = T) +
      coord_sf(expand = F,
               xlim = xlims,
               ylim = ylims) +
      facet_grid(. ~ ttl) +
      labs(y = 'KO - WT',
           x = 'MT - WT') +
      theme(panel.background = element_rect(fill = "black"),
            panel.grid = element_blank(),
            legend.position = 'none',
            plot.background = element_blank(),
            text = element_text(color = 'black'),
            axis.ticks = element_line(color = 'black'),
            plot.title = element_blank(),
            axis.text = element_text(color = 'black'),
            strip.text = element_text(color = "white"),
            strip.background = element_rect(fill = 'black'))
    p$coordinates$aspect <- function(f) { NULL }
    if (x$mark[1] == 'H3K36me2') {
      p +
        theme(axis.title.x = element_blank())
    } else if (x$mark[1] == 'H3K27me3') {
      p +
        theme(axis.title.y = element_blank())
    } else {
      p +
        theme(axis.title = element_blank())
    }
  }, simplify = F)

wrap_plots(ps[c('H3K36me2', 'H3K27me3', 'DNAme')]) %>%
  ggsave('figs/2e.pdf', ., height = 2.1, width = 4.5)
