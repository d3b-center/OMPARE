# Author: Komal S. Rathi
# Function to create cancer specific gene list

library(tidyverse)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# from oncokb.org
cancerGenes <- read.delim(file.path(ref_input_file_dir, "CancerGeneList.tsv"), stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  dplyr::select(-c(Entrez.Gene.ID, GRCh37.Isoform, GRCh37.RefSeq, GRCh38.Isoform, GRCh38.RefSeq,
                   X..of.occurrence.within.resources..Column.D.J.)) %>%  
  dplyr::rename(Gene_Symbol = Hugo.Symbol) %>%
  gather(key = "type", value = "file", -Gene_Symbol) %>%
  filter(file != "No") %>%
  mutate(file = type) 

# from annofuse
geneListRef <- read.delim(file.path(ref_input_file_dir, "genelistreference.txt"), stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)

# save output
saveRDS(cancerGenes, file.path(ref_dir, "cancer_gene_list.rds"))
