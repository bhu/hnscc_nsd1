library(tidyverse)
library(rtracklayer)

s <- list.files('data', pattern = '.10kb.') %>%
  grep('diff|WT|MT', ., invert = T, value = T) %>%
  tibble(f = .) %>%
  separate(f, c('samp', 'mark', NA, NA, NA), '\\.', F) %>%
  filter(!(samp %in% c('Det562_KO1', 'Det562_KO4', 'Cal27_KO17')) &
           mark == 'H3K36me2') %>%
  mutate(cond = case_when(grepl('KO', samp) ~ 'KO',
                          grepl('^S|^B', samp) ~ 'MT',
                          T ~ 'WT'),
         line = sub('_.*', '', samp),
         f = file.path('data', f),
         mark = sub('WGBS', 'DNAme', mark)) %>%
  arrange(cond) %>%
  filter(cond != 'MT')

d <- deframe(s[,c('samp', 'f')]) %>%
  lapply(import.bw)

r <- lapply(d, function(y) y[y$score != 0]) %>%
  Reduce(function(a, b) a[overlapsAny(a, b)], .) %>%
  granges()

split(d, s$line) %>%
  lapply(function(x) {
    o <- lapply(x, function(y) {
      findOverlaps(r, y) %>%
        to() %>%
        {y[.]} %>%
        score()
    }) %>%
      bind_cols() %>%
      `names<-`(c('x', 'y'))
    ok <- o$x > quantile(o$x, .01) &
      o$x < quantile(o$x, .99) &
      o$y > quantile(o$y, .01) &
      o$y < quantile(o$y, .99)
    write_csv(o[ok,], sprintf('misc/%s.csv', names(x[2])), col_names = F)
    export.bed(r[ok], sprintf('misc/%s.bed', names(x[2])))
  })
