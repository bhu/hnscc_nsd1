library(tidyverse)
library(readxl)
library(SummarizedExperiment)
library(DESeq2)
library(RRHO2)


load('misc/g33.rda')
e2g <- g33[g33$type == 'gene'] %>%
  {setNames(.$gene_name, .$ID)}
e2g2 <- setNames(e2g, sub('\\..*', '', names(e2g)))

load('tcga/DNA.rda')
load('tcga/RNA.rda')
RNA <- data

clu <- readxl::read_xlsx('tcga/papillon2017.xlsx', skip = 1) %>%
  dplyr::select(bc = bcr_patient_barcode,
                grp = methylation_group,
                loc = simple_tissue,
                cig = tobacco_smoking_history,
                age = days_to_birth,
                sex = gender, K36M) %>%
  dplyr::filter(!(grp %in% c('NA', 'HPV+')) & !K36M &
                  age != '[Not Available]' &
                  !(cig %in% c('[Unknown]', '[Not available]'))) %>%
  dplyr::select(-K36M) %>%
  mutate(age = as.numeric(age),
         cig = paste0('t', cig),
         nsd1 = case_when(grp == 'H3K36' ~ 'neg', T ~ 'pos') %>%
           factor(c('neg', 'pos')))

rsub <- RNA[,RNA$patient %in% clu$bc & !is.na(RNA$patient)]
hit <- match(rsub$patient, clu$bc)
rsub$grp <- clu$grp[hit]
rsub$loc <- clu$loc[hit]
rsub$nsd1 <- clu$nsd1[hit]
rsub$cig <- clu$cig[hit]
rsub$age <- scale(clu$age[hit])
rsub$sex <- clu$sex[hit]

dds <- DESeqDataSet(rsub, design = ~nsd1 + loc + cig + age + sex)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit = 1000)

res.tcga <- results(dds, name = "nsd1_pos_vs_neg") %>%
  {cbind(as.data.frame(lfcShrink(dds, coef = "nsd1_pos_vs_neg",
                                 type = 'apeglm',
                                 lfcThreshold = 0.5, svalue = T)),
         as.data.frame(.))} %>%
  `colnames<-`(c('baseMean', 'log2FoldChangeAdj', 'lfcSEadj',
                 'svalue', 'bm2', 'log2FoldChange', 'lfcSE',
                 'stat', 'pvalue', 'padj')) %>%
  rownames_to_column('gene') %>%
  mutate(gene = sub('\\..*', '', gene),
         name = e2g2[gene]) %>%
  dplyr::select(-bm2)

save(res.tcga, file = 'tcga/deg.rda')
