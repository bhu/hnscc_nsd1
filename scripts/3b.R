library(reshape2)
library(tidyverse)
library(rtracklayer)
library(ggrepel)


load('misc/cons_ko.rda')

uni <- c('Cal27', 'Det562', 'FaDu') %>%
  sprintf('misc/%s.bed', .) %>%
  lapply(import.bed) %>%
  GRangesList() %>%
  unlist() %>%
  granges() %>%
  unique() 

signif.num <- function(x) {
  as.character(
    symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "**", "*", "˙", "NS"))
  )
}

g <- import.bed('ensembl/gene.bed')
ig <- import.bed('ensembl/intergenic.bed')

rs <- list.files('ensembl', full.names = T) %>%
  setNames(sub('ensembl/(.*).bed', '\\1', .)) %>%
  lapply(import.bed) %>%
  GRangesList()

ol.min <- 1

ps <- c('Genome-wide', 'Genic', 'Intergenic') %>%
  setNames(., .) %>%
  lapply(function(x) {
    uSet <- uni
    mapply(function(qSet, clu) {
      if (x == "Intergenic") {
        qSet <- qSet[overlapsAny(qSet, ig) & !overlapsAny(qSet, g)]
        uSet <- uSet[overlapsAny(uSet, ig) & !overlapsAny(uSet, g)]
      } else if (x == "Genic") {
        qSet <- qSet[overlapsAny(qSet, g) & !overlapsAny(qSet, ig)]
        uSet <- uSet[overlapsAny(uSet, g) & !overlapsAny(uSet, ig)]
      }
      
      ql <- length(qSet)
      ul <- length(uSet)
      ol.u <- countOverlaps(rs, uSet, minoverlap = ol.min)
      
      countOverlaps(rs, qSet, minoverlap = ol.min) %>%
        data.frame(support = .) %>%
        rownames_to_column('region') %>%
        mutate(b = ol.u[region] - support,
               c = ql - support,
               d = ul - support - b - c,
               userSet = clu) %>%
        cbind(., apply(., 1, function(y) {
          fisher.test(matrix(as.numeric(y[2:5]), 2, 2),
                      alternative = "two.sided")[c("estimate", "conf.int")] %>%
            unlist() %>%
            c(fisher.test(matrix(as.numeric(y[2:5]), 2, 2),
                          alternative = "greater")$p.value)
        }) %>%
          t() %>%
          `colnames<-`(c('or','clo', 'chi', 'p'))
        ) %>%
        left_join(sapply(rs, length) %>%
                    data.frame(size = .) %>%
                    rownames_to_column("region"),
                  by = "region") %>%
        mutate(q = p.adjust(p, method = 'fdr')) %>%
        dplyr::filter(q < 0.25 & support > 100) %>%
        dplyr::select(clus = userSet, sup = support, reg = region,
                      qval = q, clo, chi, or) %>%
        top_n(5, -qval) %>%
        arrange(or) %>%
        mutate(sig = signif.num(qval),
               reg = fct_inorder(reg),
               ttl = sprintf('Cluster %s (%s)', clu, x)) %>%
        ggplot(aes(x = or, y = reg)) +
        geom_segment(aes(x = chi, y = reg, xend = clo, yend = reg,
                         color = sup)) +
        geom_point(aes(size = sup, color = sup), stat = "identity") +
        geom_text(aes(label = sig), vjust = -0.4) +
        scale_size(name = "# overlaps") +
        labs(y = "Region", x = "Odds ratio") +
        scale_color_viridis_c(name = "# overlaps",
                              begin = 0.2, end = 0.8,
                              option = "A", guide = "legend",
                              direction = 1) +
        facet_wrap(. ~ ttl) +
        coord_cartesian(clip = "off") +
        theme(panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank(),
              axis.line.x = element_line(size = 0.5, color = "black"),
              axis.text = element_text(color = "black"),
              panel.grid.major.y = element_line(color = "grey", linetype = "solid"),
              axis.ticks.y = element_line(color = "grey", linetype = "solid"),
              panel.grid.major.x = element_line(color = "grey", linetype = "dashed"),
              legend.key = element_rect(fill = "transparent", color = "transparent"),
              legend.background = element_blank(),
              strip.text.x = element_text(color = "white"),
              strip.background.x = element_rect(fill = "black"),
              legend.position = "none")
    }, cons, names(cons), SIMPLIFY = F)
  })

ggsave('figs/3b.pdf', ps$Intergenic$B, height = 2.2, width = 3,
       device = cairo_pdf, bg = "transparent")
