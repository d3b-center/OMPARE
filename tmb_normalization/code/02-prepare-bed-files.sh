# This script is adapted from the https://github.com/AlexsLemonade/OpenPBTA-analysis/
# blob/master/analyses/tcga-capture-kit-investigation

#!/bin/bash
####################################
# This is script that reads the bam_manifest.tsv and download
# all the uniq BED and add chr prefix to it.
#
# BAMs with `|` in the capture kit name and url are those with 
# more than one capture kit which neither the GDC nor its origin 
# data center could retrieve/figure out what the actual capture kit 
# had been applied. We should just generate intersec BED for those 
# samples and used that for our analysis.
#
# Those are all hg19 based coordinates, use UCSC online LiftOVer 
# to convert to GRCh38 based bed. Add chr to the TCGA bed to 
# make sure its format is liftover compatible.
# The second step download the chain file from UCSC, converts the  
# BED files from hg19 to Gh38 ans sorts and mereges (in case there 
# are  any overlaps in the BED regions)
####################################

# root and results directories
root_dir="$(git rev-parse --show-toplevel)"
results_dir="$root_dir/tmb_normalization/results"
ref_dir="$root_dir/tmb_normalization/references"
mkdir -p $results_dir
mkdir -p ref_dir

# tmp directory and some sub-directories
scratch_dir="/tmp"
tcga_not_in_pbta_dir="$scratch_dir/tcga_not_in_pbta"
mkdir -p $tcga_not_in_pbta_dir
tcga_in_pbta_dir="$scratch_dir/tcga_in_pbta"
mkdir -p $tcga_in_pbta_dir

# the following part deals with the bed files from TCGA specimens that are NOT part of PBTA
tcga_tsv="$results_dir/tcga_not_in_pbta_bam_manifest.tsv"
sed 1d $tcga_tsv \
| cut -f9 | tr "\|" "\n" | sort -u \
| while read i
do
    filename=`basename $i`
    curl $i | awk '{print "chr"$0}' > "$tcga_not_in_pbta_dir/$filename"
done

# note that by using 'curl' - not all the bed files can be found -
# files with 0kb were deleted and I went into tcga_not_in_pbta_unique_bed_files.tsv file
# - try to see whether I can manually find the bed file
# for some of the samples, multiple bed files/capture kits are listed since GDC does not know
# exactly which one was used - for those, some of them does not have a bed file available
# hence, manually choose the kit that actually has an associated bed file and the selected bed
# is recorded by adding another column called "bed_selected" in the tcga_not_in_pbta_bed_file.tsv
# and save it as "tbga_not_in_pbta_bed_selected.tsv". If among the multiple kits that are listed 
# as available, the one with the highest coverage is chosen. 
# the bed file name that will be used in the subsequent events

# the following part deals with the bed files from TCGA specimens that ARE part of PBTA
# the code to get the manifest is https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/
# master/analyses/tcga-capture-kit-investigation/scripts/get-tcga-capture_kit.py
# the result was saved as '../results/tcga_in_pbta_bam_manifest.tsv'
# Alternatively, the bed files can be downloaded from here: https://github.com/AlexsLemonade/
# OpenPBTA-analysis/tree/master/analyses/tcga-capture-kit-investigation/results/bedfiles

# TCGA in PBTA
tcga_in_pbta_tsv="$results_dir/tcga_in_pbta_bam_manifest.tsv"
sed 1d $tcga_in_pbta_tsv\
| cut -f3 | tr "\|" "\n" | sort -u \
| while read i
do
    filename=`basename $i`
    curl -s $i | awk '{print "chr"$0}' > "$tcga_in_pbta_dir/$filename"
done


########################################################################################
## Downloading chain.gz file for CrossMap command
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz  -O "$scratch_dir/hg19ToHg38.over.chain.gz"

## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls $tcga_in_pbta_dir/*.bed`; do
  out=$(echo $i | sed 's/.bed/.Gh38.bed/g')
  CrossMap.py bed "$scratch/hg19ToHg38.over.chain.gz" $i $out
done

## This command uses sort and merge from bedtools to merge any overalpping BED regions
tcga_in_pbta_outputdir="$results_dir/bed_files/tcga_in_pbta"
mkdir -p $tcga_in_pbta_outputdir
for i in `ls $tcga_in_pbta_dir/*.Gh38.bed`; do
  bedtools sort -i $i \
  | bedtools merge \
  > "$tcga_in_pbta_outputdir/$(basename $i)"
done


## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
for i in `ls $tcga_not_in_pbta_dir/*.bed`; do
  out=$(echo $i | sed 's/.bed/.Gh38.bed/g')
  CrossMap.py bed "$scratch_dir/hg19ToHg38.over.chain.gz" $i $out
done

## This command uses sort and merge from bedtools to merge any overalpping BED regions
tcga_not_in_pbta_outputdir="$results_dir/bed_files/tcga_not_in_pbta"
mkdir -p $tcga_not_in_pbta_outputdir
for i in `ls $tcga_not_in_pbta_dir/*.Gh38.bed`; do
  bedtools sort -i $i \
  | bedtools merge \
  > "$tcga_not_in_pbta_outputdir/$(basename $i)"
done


# the following deals with bed files used for PNOC008 samples
# the bed files were provided and was put into the scratch file without having to download it anywhere
## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
pnoc008_dir="$scratch_dir/pnoc008"
mkdir -p $pnoc008_dir
for i in `ls $pnoc008_dir/*.bed`; do
  out=$(echo $i | sed 's/.bed/.Gh38.bed/g')
  CrossMap.py bed "$scratch_dir/hg19ToHg38.over.chain.gz" $i $out
done

## This command uses sort and merge from bedtools to merge any overalpping BED regions
for i in `ls $pnoc008_dir/*.Gh38.bed`; do
  bedtools sort -i $i \
  | bedtools merge \
  > "$ref_dir/$(basename $i)"
done

# the following deals with bed files used for PBTA samples
# the bed files were provided and was put into the scratch file without having to download it anywhere
## This command takes every BED files downloaded and uses CrossMap tool to convert hg19 coordinates to Gh38
pbta_dir="$scratch_dir/pbta"
mkdir -p $pbta_dir
for i in `ls $pbta_dir/*.bed`; do
  out=$(echo $i | sed 's/.bed/.Gh38.bed/g')
  CrossMap.py bed "$scratch_dir/hg19ToHg38.over.chain.gz" $i $out
done

## This command uses sort and merge from bedtools to merge any overalpping BED regions
mkdir -p "$ref_dir/pbta_bed_files"
for i in `ls $pbta_dir/*.Gh38.bed`; do
  bedtools sort -i $i \
  | bedtools merge \
  > "$ref_dir/pbta_bed_files/$(basename $i)"
done

