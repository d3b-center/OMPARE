# Author: Komal S. Rathi
# Function to create cancer specific gene list

library(tidyverse)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# from oncokb.org
cancerGenes <- read.delim(file.path(data_dir, "input_files", "CancerGeneList.tsv"), stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  dplyr::select(-c(Entrez.Gene.ID, GRCh37.Isoform, GRCh37.RefSeq, GRCh38.Isoform, GRCh38.RefSeq,
                   X..of.occurrence.within.resources..Column.D.J.)) %>%  
  dplyr::rename(Gene_Symbol = Hugo.Symbol) %>%
  gather(key = "type", value = "file", -Gene_Symbol) %>%
  filter(file != "No") %>%
  mutate(file = type) 

# from annofuse
geneListRef <- read.delim(file.path(data_dir, "input_files", "genelistreference.txt"), stringsAsFactors = F)
geneListRef <- geneListRef[grep("TumorSuppressorGene|CosmicCensus|Oncogene", geneListRef$type),]
cancerGenes <- rbind(cancerGenes, geneListRef)

# save output
saveRDS(cancerGenes, file.path(data_dir, "cancer_gene_list.rds"))
