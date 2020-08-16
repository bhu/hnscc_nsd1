#!/bin/bash

computeMatrix reference-point -S data/diff.K36me2.100.logcpm_ms.bw \
  data/diff.WGBS.1kb.beta.bw \
  data/diff.K27ac.100.logcpm_ms.bw \
  data/diff.K27me3.100.logcpm_ms.bw \
  -R misc/b_cre.bed \
  -o data/b_cre.mat.gz \
  --samplesLabel H3K36me2 DNAme H3K27ac H3K27me3 \
  -a 50000 -b 50000 -bs 100 -bl misc/bl.bed \
  --missingDataAsZero --skipZeros -p 20

