library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

r <- import.bed('misc/Cal27.bed')
d <- read_csv('misc/Cal27.csv', col_names = F)
d$clu <- read_csv('misc/clus.Cal27.5000.5000.txt', col_names = F)$X1

ctr <- d %>%
  dplyr::filter(clu != -1) %>%
  group_by(clu) %>%
  summarise_all(mean) %>%
  mutate(mu = 0.5 * (X1 + X2)) %>%
  arrange(mu)

a1 <- r[d$clu == ctr$clu[1]]
b1 <- r[d$clu == ctr$clu[2]]
c1 <- r[d$clu == ctr$clu[3]]


r <- import.bed('misc/Det562.bed')
d <- read_csv('misc/Det562.csv', col_names = F)
d$clu <- read_csv('misc/clus.Det562.5000.5000.txt', col_names = F)$X1

ctr <- d %>%
  dplyr::filter(clu != -1) %>%
  group_by(clu) %>%
  summarise_all(mean) %>%
  mutate(mu = 0.5 * (X1 + X2)) %>%
  arrange(mu)

ab2 <- r[d$clu == ctr$clu[1]]
c2 <- r[d$clu == ctr$clu[2]]

r <- import.bed('misc/FaDu.bed')
d <- read_csv('misc/FaDu.csv', col_names = F)
d$clu <- read_csv('misc/clus.FaDu.5000.5000.txt', col_names = F)$X1

ctr <- d %>%
  dplyr::filter(clu != -1) %>%
  group_by(clu) %>%
  summarise_all(mean) %>%
  mutate(mu = 0.5 * (X1 + X2)) %>%
  arrange(mu)

b3 <- r[d$clu == ctr$clu[1]]
c3 <- r[d$clu == ctr$clu[2]]

d$clu <- read_csv('misc/clus.FaDu.1000.1000.txt', col_names = F)$X1

ctr <- d %>%
  dplyr::filter(clu != -1) %>%
  group_by(clu) %>%
  summarise_all(mean) %>%
  mutate(mu = 0.5 * (X1 + X2)) %>%
  arrange(mu)

a3 <- r[d$clu == ctr$clu[1]]

cons <- list(A = list(a1, ab2, a3),
             B = list(b1, ab2, b3),
             C = list(c1, c2, c3)) %>%
  lapply(function(x) {
    Reduce(function(a,b){a[overlapsAny(a,b)]}, x) %>%
      granges()
  })

save(cons, file = 'misc/cons.rda')
load('misc/cons_ko.rda')

read_delim('misc/GeneHancer_double_elite.bed', '\t', col_names = F) %>%
  dplyr::filter(X11 == 'Enhancer') %>%
  dplyr::select(chr = X1, start = X2, end = X3) %>%
  makeGRangesFromDataFrame() %>%
  {.[overlapsAny(., cons$B)]} %>%
  as_tibble() %>%
  dplyr::select(1:3) %>%
  mutate(start = as.integer(start),
         end = as.integer(end)) %>%
  write_delim('misc/b_cre.bed', '\t', col_names = F)
