library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)
library(tximport)
library(DESeq2)
library(AnnotationDbi)
library(zFPKM)

txdb <- makeTxDbFromGFF('misc/gencode.v33.annotation.gff3.gz')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
save(tx2gene, file = "misc/tx2gene.rda")

g33 <- import.gff3('misc/gencode.v33.annotation.gff3.gz')
save(g33, file = 'misc/g33.rda')

samps <- list.files('data', pattern = '.sf.txt.gz$', full.names = T) %>%
  setNames(., sub('.*/(.*).sf.*', '\\1', .)) 

mdat <- data.frame(samp = names(samps)) %>%
  mutate(type = case_when(
    grepl('KO',samp) ~ 'KO',
    grepl('^S|^B', samp) ~ 'MT',
    TRUE ~ 'WT'
  )) %>%
  mutate(line = factor(sub('_.*', '', samp)),
         type = factor(type)) %>%
  column_to_rownames("samp")

salmon <- tximport(samps[names(samps) %in% rownames(mdat)],
                   type = "salmon", tx2gene = tx2gene)
salmon$length[salmon$length == 0] <- 1

dds <- DESeqDataSetFromTximport(salmon,
                                colData = mdat,
                                design = ~type)

pcg <- g33[g33$type == 'gene' & g33$gene_type == 'protein_coding']
expd <- zFPKM(data.frame(fpkm(dds))) %>%
  `colnames<-`(colnames(dds)) %>%
  {apply(., 2, function(x) {
    rownames(.)[x > 0] %>%
      {pcg[pcg$ID %in% .]}
  })}

save(expd, file = 'misc/expd.rda')
