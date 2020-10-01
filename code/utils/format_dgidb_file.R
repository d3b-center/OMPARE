# Author: Komal S. Rathi
# Function: Format DGIDB file

library(tidyverse)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

dgidb <- read.delim(file.path(ref_input_file_dir, "DGIdb.txt"), stringsAsFactors = F)
dgidb <- dgidb %>%
  mutate(Gene_name = gene_name) %>%
  filter(interaction_types != "" & drug_name != "" & Gene_name != "") %>%
  dplyr::select(Gene_name, drug_name) %>%
  unique() %>%
  group_by(Gene_name) %>%
  dplyr::summarize(Drugs=paste(drug_name, collapse=", ")) %>%
  as.data.frame()

# save output
saveRDS(dgidb, file = file.path(ref_dir, "dgidb_output.rds"))
