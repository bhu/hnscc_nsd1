#!/usr/bin/env bash

bash scripts/anns.sh
Rscript --vanilla scripts/tracks.R
bash scripts/cm.sh
Rscript --vanilla scripts/active.R
Rscript --vanilla scripts/clust.R
bash scripts/clust.sh
Rscript --vanilla scripts/cons.R
bash scripts/cm2.sh
Rscript --vanilla scripts/dge.R
Rscript --vanilla scripts/tcga.R
Rscript --vanilla scripts/dmr.R
Rscript --vanilla scripts/dge2.R

Rscript --vanilla scripts/1a.R
bash scripts/1b.sh
bash scripts/1c.sh
Rscript --vanilla scripts/1d.R
Rscript --vanilla scripts/1e.R
Rscript --vanilla scripts/2a.R
bash scripts/2b.sh
Rscript --vanilla scripts/2d.R
Rscript --vanilla scripts/2e.R
Rscript --vanilla scripts/3a.R
Rscript --vanilla scripts/3b.R
Rscript --vanilla scripts/3c.R
Rscript --vanilla scripts/3d.R
Rscript --vanilla scripts/3e.R
Rscript --vanilla scripts/4a.R
Rscript --vanilla scripts/4b.R
Rscript --vanilla scripts/4c.R
Rscript --vanilla scripts/4d.R
Rscript --vanilla scripts/5a.R
Rscript --vanilla scripts/5b.R
Rscript --vanilla scripts/5c.R
Rscript --vanilla scripts/6a.R
Rscript --vanilla scripts/6b.R
Rscript --vanilla scripts/6c.R
Rscript --vanilla scripts/6d.R
