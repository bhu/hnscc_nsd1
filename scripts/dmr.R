library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(sva)
library(limma)
library(dmrff)

load('tcga/DNA.rda')
load('tcga/DNAme.rda')

r <- rowRanges(data) %>%
  as.data.frame() %>%
  dplyr::rename(ref = Composite.Element.REF) %>%
  dplyr::select(chr = seqnames, start, end, ref) %>%
  as.data.table()

d <- as.data.table(assays(data)[[1]])

muts <- lapply(DNA, function(x) {
  x[,c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification')] %>%
    dplyr::filter(Hugo_Symbol == 'NSD1')
}) %>% bind_rows(.id = 'Caller') %>%
  dplyr::filter(Variant_Classification != 'Silent') %>%
  pull(Tumor_Sample_Barcode) %>%
  sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', .) %>%
  unique()

clu <- read_xlsx('tcga/papillon2017.xlsx', skip = 1) %>%
  dplyr::select(bc = bcr_patient_barcode, grp = methylation_group,
                loc = simple_tissue, K36M,
                cig = tobacco_smoking_history,
                age = days_to_birth,
                sex = gender) %>%
  dplyr::filter(!(grp %in% c('NA', 'HPV+')) & !K36M &
                  age != '[Not Available]' &
                  !(cig %in% c('[Unknown]', '[Not available]'))) %>%
  dplyr::select(-K36M) %>%
  mutate(age = as.numeric(age),
         cig = paste0('t', cig),
         nsd1 = case_when(grp == 'H3K36' ~ 'neg', T ~ 'pos') %>%
           factor(c('pos', 'neg')))

ok <- as.character(seqnames(rowRanges(data))) != '*' &
  start(rowRanges(data)) > 0

methylation <- colnames(data) %>%
  sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', .) %>%
  {. %in% clu$bc} %>%
  {assays(data[, .])[[1]]} %>%
  {.[ok,]} %>%
  {.[apply(., 1, function(x) {any(is.finite(x))}),]}

md <- tibble(samp = colnames(methylation)) %>%
  mutate(bc = sub('^(TCGA-[A-Z0-9]*-[A-Za-z0-9]*).*', '\\1', samp)) %>%
  left_join(clu, by = 'bc')

mod <- model.matrix(~ nsd1 + loc + age + sex + cig, md)
mod0 <- mod[,1]
set.seed(42)
random.idx <- sample(1:nrow(methylation), 10000)
methylation.sva <- methylation[random.idx,]
methylation.mean <- rowMeans(methylation.sva, na.rm=T)
idx <- which(is.na(methylation.sva), arr.ind=T)
if (nrow(idx) > 0)
  methylation.sva[idx] <- methylation.mean[idx[,"row"]]
sva.fit <- sva(methylation.sva, mod=mod, mod0=mod0)
design <- cbind(mod, sva.fit$sv)
fit <- lmFit(methylation, design)
fit <- eBayes(fit)

stats <- tibble(estimate = fit$coefficients[,"nsd1neg"],
                se = sqrt(fit$s2.post) * fit$stdev.unscaled[,"nsd1neg"],
                p.value = fit$p.value[,"nsd1neg"])

locs <- rowRanges(data) %>%
  as.data.frame() %>%
  dplyr::select(chr = seqnames, start, end,
                ref = Composite.Element.REF) %>%
  {.[match(rownames(methylation), .$ref),]}
stats$chr <- as.character(locs$chr)
stats$pos <- round(0.5 * (locs$start + locs$end))
dmrs <- dmrff(estimate=stats$estimate,
              se=stats$se,
              p.value=stats$p.value,
              methylation=methylation,
              chr=stats$chr,
              pos=stats$pos,
              maxgap=500,
              verbose=T)

save(stats, dmrs, file = 'tcga/dmr.rda')
