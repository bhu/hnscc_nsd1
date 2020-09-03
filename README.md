* Processed data will first need to be downloaded from GEO and placed in the directory named `data` , then files should be renamed to remove prefixes (e.g., `for f in GSM*; do mv $f ${f#GSM*_} ; done`). Spreadsheet containing mass spectrometry results should also be added there and renamed `ms.xlsx`.
* The following programs are required:
  * [deepTools](https://github.com/deeptools/deepTools)
  * [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)
  * [UCSC Genome Browser Utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)
    * gtfToGenePred
    * genePredToBed
  * [Python](https://www.python.org/downloads/)
    * [HDBSCAN](https://github.com/scikit-learn-contrib/hdbscan)
    * [NumPy](https://pypi.org/project/matplotlib/)
    * [matplotlib](https://pypi.org/project/matplotlib/)
    * [seaborn](https://pypi.org/project/seaborn/)
  * [R](https://cran.r-project.org/)
    * Bioconductor
      * AnnotationDbi
      * ChIPseeker
      * DESeq2
      * DiffBind
      * GenomicFeatures
      * GenomicRanges
      * SummarizedExperiment
      * TCGAbiolinks
      * TxDb.Hsapiens.UCSC.hg38.knownGene
      * fgsea
      * limma
      * org.Hs.eg.db
      * rtracklayer
      * sva
      * tximport
      * zFPKM
    * CRAN
      * MASS
      * cowplot
      * data.table
      * forestplot
      * furrr
      * ggforce
      * ggformula
      * ggnewscale
      * ggrepel
      * ggsignif
      * hexbin
      * isoband
      * lwgeom
      * msigdbr
      * pals
      * patchwork
      * purrr
      * readxl
      * reshape2
      * RobustRankAggreg
      * scales
      * sf
      * tidyverse
      * viridis
    * GitHub
      * [RRHO2](https://github.com/RRHO2/RRHO2)
      * [dmrff](https://github.com/perishky/dmrff)
      * [ggrastr](https://github.com/VPetukhov/ggrastr)
* The scripts below should be ran to pre-process a number of files 
  * `anns.sh` downloads various resources (e.g., gene annotations)
  * `tracks.R` computes difference and average tracks (e.g., KO-PA, all WT, etc.)
  * `cm.sh` generates the enrichment matrices (e.g., heatmaps in 1c / 2c)
  * `active.R` determines actively transcribed genes in each sample
  * `clust.R` prepares H3K36me2 PA-KO scatter plot data for HDBSCAN
  * `clust.sh` runs HDBSCAN to obtain clusters of concordant epigenetic change
  * `cons.R` derives consensus cluster assignments from all three cell lines
  * `cm2.sh` creates enrichment matrices centered around cluster B CREs
  * `dge.R` performs differential gene expression analysis
  * `tcga.R` obtains TCGA-HNSC datasets
  * `dmr.R` calls differentially methylated regions in TCGA-HNSC based on NSD1 status
  * `dge2.R` identifies differentially expressed genes between NSD1+/- TCGA-HNSC samples
* One option is to simply use `all.sh` to run everything (may take several hours depending on the connection speed to download public data / resources, as the actual code execution should total less than an hour), or individual figures can be produced from specific scripts found within the `scripts` directory. However, some steps depend on data saved during the processing of preceeding plots and so it is suggested to run everything in sequence at least once
* Certain figures as generated here will largely resemble but are not exactly the same as ones included in the manuscript, as those have had additional manual aesthetic adjustments. Additionally, some minor discrepancies may arise due to a degree of stochasticity in particular algorithms, but they do not impact the conclusions presented

