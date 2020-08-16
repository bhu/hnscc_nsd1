library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create('tcga', showWarnings = FALSE)
setwd('tcga')

query <- GDCquery(project = "TCGA-HNSC", 
                  legacy = FALSE,
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
RNA <- GDCprepare(query = query,
                  save = TRUE, 
                  save.filename = "RNA.rda")

clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
save(clin, file = 'clin.rda')

d <- c('muse', 'varscan2', 'somaticsniper', 'mutect2') 
d <- setNames(d, d)
d <- lapply(d, function(x) {
  GDCquery_Maf('HNSC', pipelines = x)
}) 
DNA <- d
save(DNA, file = 'DNA.rda')

query <- GDCquery(project = "TCGA-HNSC",
                  legacy = FALSE,
                  data.category = "DNA Methylation")
GDCdownload(query)
DNAme <- GDCprepare(query = query,
                    save = TRUE,
                    save.filename = "DNAme.rda",
                    summarizedExperiment = TRUE)

download.file('https://api.gdc.cancer.gov/data/7a3d7067-09d6-4acf-82c8-a1a81febf72c',
              'atac_samps.txt')
download.file('https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec',
              'atac_peaks.txt')
download.file('https://api.gdc.cancer.gov/data/a544ed6f-4a27-4430-8f99-48b657da11fe',
              'atac_cts.rds')
