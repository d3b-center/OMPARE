# Author: Komal S. Rathi
# Function to create cancer specific gene list

library(tidyverse)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# from Pichai
cancerGenes <- read.delim(file.path(ref_input_file_dir, "CancerGeneList.tsv"), stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  filter(Gene_Symbol != "") %>%
  dplyr::select(-Count) %>%
  gather(key = "file", value = "type", -Gene_Symbol) %>%
  mutate(type = file)

# from annofuse
geneListRef <- read.delim(file.path(ref_input_file_dir, "genelistreference.txt"), stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)

# save output
saveRDS(cancerGenes, file.path(ref_dir, "cancer_gene_list.rds"))
