#!/bin/bash

computeMatrix scale-regions -S data/Cal27.K36me2.100.logcpm_ms.bw \
  data/Det562.K36me2.100.logcpm_ms.bw \
  data/FaDu.K36me2.100.logcpm_ms.bw \
  data/BICR78.K36me2.100.logcpm_ms.bw \
  data/SCC4.K36me2.100.logcpm_ms.bw \
  data/SKN3.K36me2.100.logcpm_ms.bw \
  -R ensembl/intergenic.bed \
  -o data/WTvsMT.K36me2.mat.gz \
  --samplesLabel Cal27 Det562 FaDu BICR78 SCC4 SKN3 \
  -m 20000 -a 20000 -b 20000 -bs 1000 -bl misc/bl.bed \
  --missingDataAsZero --skipZeros -p 20

computeMatrix scale-regions -S data/Cal27.K27me3.100.logcpm_ms.bw \
  data/Det562.K27me3.100.logcpm_ms.bw \
  data/FaDu.K27me3.100.logcpm_ms.bw \
  data/BICR78.K27me3.100.logcpm_ms.bw \
  data/SCC4.K27me3.100.logcpm_ms.bw \
  data/SKN3.K27me3.100.logcpm_ms.bw \
  -R ensembl/intergenic.bed \
  -o data/WTvsMT.K27me3.mat.gz \
  --samplesLabel Cal27 Det562 FaDu BICR78 SCC4 SKN3 \
  -m 20000 -a 20000 -b 20000 -bs 1000 -bl misc/bl.bed \
  --missingDataAsZero --skipZeros -p 20

computeMatrix scale-regions -S data/Cal27.WGBS.1kb.beta.bw \
  data/Det562.WGBS.1kb.beta.bw \
  data/FaDu.WGBS.1kb.beta.bw \
  data/BICR78.WGBS.1kb.beta.bw \
  data/SCC4.WGBS.1kb.beta.bw \
  data/SKN3.WGBS.1kb.beta.bw \
  -R ensembl/intergenic.bed \
  -o data/WTvsMT.WGBS.mat.gz \
  --samplesLabel Cal27 Det562 FaDu BICR78 SCC4 SKN3 \
  -m 20000 -a 20000 -b 20000 -bs 1000 -bl misc/bl.bed \
  --missingDataAsZero --skipZeros -p 20


plotHeatmap -m data/WTvsMT.K36me2.mat.gz -o figs/1c.pdf 
