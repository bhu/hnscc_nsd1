library(tidyverse)
library(rtracklayer)

fs <- list.files('data', pattern = '.bw$') %>%
  tibble(f = .) %>%
  separate(f, c('samp', 'mark', 'res', 'norm' , NA), '\\.', F) %>%
  mutate(cond = case_when(grepl('KO', samp) ~ 'KO',
                          grepl('Cal|Det|FaDu', samp) ~ 'WT',
                          T ~ 'MT'),
         line = sub('_.*', '', samp)) %>%
  filter(cond != 'KO' | samp %in% c('Cal27_KO1',
                                    'Det562_KO2',
                                    'FaDu_KO1'))

fs[fs$cond != 'KO',] %>%
  by(., .$mark, function(x) {
    by(x, x$cond, function(y) {
      d <- lapply(file.path('data', y$f), import.bw)
      r <- lapply(d, granges) %>%
        GRangesList() %>%
        unlist() %>%
        unique() %>%
        sort()
      r$score <- lapply(d, function(z) {
        findOverlaps(r, z) %>%
          as("List") %>%
          extractList(z$score, .) %>%
          mean()
      }) %>%
        bind_cols() %>%
        rowMeans(na.rm = T)
      export.bw(r, sprintf('data/%s.%s.%s.%s.bw', y$cond,
                           y$mark, y$res, y$norm)[1])
    })
  })

fs[fs$cond != 'MT',] %>%
  by(., .$mark, function(x) {
    d <- by(x, x$line, function(y) {
      d <- arrange(y, cond) %>%
        pull(f) %>%
        file.path('data', .) %>%
        lapply(import.bw)
      r <- d[[1]][overlapsAny(d[[1]], d[[2]])] %>%
        granges()
      r$score <- lapply(d, function(z) {
        findOverlaps(r, z) %>%
          as("List") %>%
          extractList(z$score, .) %>%
          mean()
      }) %>%
        {.[[1]] - .[[2]]}
      export.bw(r, sprintf('data/%s_diff.%s.%s.%s.bw', y$line,
                           y$mark, y$res, y$norm)[1])
      r
    })
    r <- lapply(d, granges) %>%
      GRangesList() %>%
      unlist() %>%
      unique() %>%
      sort()
    r$score <- lapply(d, function(z) {
      findOverlaps(r, z) %>%
        as("List") %>%
        extractList(z$score, .) %>%
        mean()
    }) %>%
      bind_cols() %>%
      rowMeans(na.rm = T)
    export.bw(r, sprintf('data/diff.%s.%s.%s.bw',
                         x$mark, x$res, x$norm)[1])
  })
