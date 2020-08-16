#!/bin/bash

export LC_COLLATE=C
mkdir -p ensembl

# Download annotations
curl ftp://ftp.ensembl.org/pub/release-97/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz | gunzip > reg.gff
curl ftp://ftp.ensembl.org/pub/release-97/gff3/homo_sapiens/Homo_sapiens.GRCh38.97.gff3.gz | gunzip > ann.gff3

# Download NCBI assembly report (contains mapping between naming conventions)
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt | recode dos..l1 > report.txt

# Generate dictionary connecting Ensembl to UCSC convention
grep -v '^#' report.txt | grep -v -i 'patch' | awk -v OFS="\t" '{if ($2=="assembled-molecule") {print $3,$NF;} else {print $5,$NF}}' > map.dict

# Get list of all regions in Ensembl annotation
regs=($(cut -f 3 ann.gff3 | grep -vE '^#|^biological_region$|^chromosome$|^supercontig$|^scaffold$' | cut -f 9 | sort | uniq))

# Generate bed file for each region
for reg in "${regs[@]}"; do
  grep -v '^#' ann.gff3 | \
    awk -v OFS='\t' -v REG="$reg" '$3==REG {print $1,$4-1,$5}' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - | \
    awk -v OFS='\t' 'NR==FNR{map[$1]=$2;next} {if($1 in map)$1=map[$1]}1' map.dict - | \
    sort -k1,1V -k2,2n > ensembl/${reg}.bed
done

# Repeat for CpG islands
grep -v '^#' ann.gff3 | \
  grep 'logic_name=cpg' | \
  awk -v OFS='\t' '{print $1,$4-1,$5}' | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i - | \
  awk -v OFS='\t' 'NR==FNR{map[$1]=$2;next} {if($1 in map)$1=map[$1]}1' map.dict - | \
  sort -k1,1V -k2,2n > ensembl/CGI.bed

# Extract chromosome/scaffold/contig sizes
grep '^##sequence-region' ann.gff3 | \
  tr -s " " | \
  cut -d " " -f2,4 | \
  awk -v OFS='\t' 'NR==FNR{map[$1]=$2;next} {if($1 in map)$1=map[$1]}1' map.dict - | \
  sort -k1,1V -k2,2n > hg38.chrom.sizes

# Create intergenic & intron based on genes/exons
bedtools subtract -a ensembl/gene.bed -b ensembl/exon.bed > ensembl/intron.bed
bedtools complement -i ensembl/gene.bed -g hg38.chrom.sizes | sort -k1,1V -k2,2n > ensembl/intergenic.bed

# Get unique region types in regulation build
regs=($(cut -f3 reg.gff | sort | uniq))

# Generate bed file for each type
for reg in "${regs[@]}"; do
  grep -v '^#' reg.gff | \
    awk -v OFS='\t' -v REG="$reg" '$3==REG {print $1,$4-1,$5}' | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - | \
    awk -v OFS='\t' 'NR==FNR{map[$1]=$2;next} {if($1 in map)$1=map[$1]}1' map.dict - | \
    sort -k1,1V -k2,2n > ensembl/${reg}.bed
done

# blacklisted regions
mkdir -p misc
awk -v OFS='\t' '{print $1,0,10"\n"$1,$2-10,$2}' hg38.chrom.sizes > misc/ends.bed
curl https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz | gzip -d | cut -f 1-3 | sort -k1,1 -k2,2n >  misc/blacklist.bed
cat misc/ends.bed misc/blacklist.bed | sort -k1,1 -k2,2n | bedtools merge -i - > misc/bl.bed
bedtools subtract -a ensembl/intergenic.bed -b misc/bl.bed -A | awk '$3-$2 > 20000' > misc/igr.bed

# gene annotations for plotting tracks
curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.basic.annotation.gtf.gz \
  | zcat \
  | grep protein_coding \
  > g33.gtf
gtfToGenePred -geneNameAsName2 -genePredExt g33.gtf g33.genePred
awk -v OFS='\t' '!arr[$12]++{$1=$12; print $0}' g33.genePred > g33.uniq.genePred
genePredToBed g33.uniq.genePred g33.bed
mv g33.bed misc/

curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz > misc/gencode.v33.annotation.gff3.gz

wget https://genecards.weizmann.ac.il/geneloc/gh_hub/hg38/GeneHancer.bb
wget https://genecards.weizmann.ac.il/geneloc/gh_hub/hg38/GH_interactions1_doubleElite.bb
wget https://genecards.weizmann.ac.il/geneloc/gh_hub/hg38/GeneHancer_double_elite.bb
bigBedToBed GeneHancer.bb misc/GeneHancer.bed
bigBedToBed GeneHancer_double_elite.bb misc/GeneHancer_double_elite.bed
bigBedToBed GH_interactions1_doubleElite.bb misc/GH_interactions1_doubleElite.bed

# clean up
rm gencode* g33.[gu]* *gff* report* map* hg38* *bb

