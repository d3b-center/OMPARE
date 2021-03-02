# merge annotations from the preceding steps with actionable genes information.
library(tidyverse)
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

patient <- 'PNOC008-28'
output_dir <- file.path(results_dir, patient, 'output')

# read annotated merged output
merged_file <- read.delim(file.path(output_dir, 'oncokb_merged_annotated.txt'))

# read actionable genes
act_genes <- read.delim(file.path(ref_dir, 'oncokb', 'oncokb_biomarker_drug_associations.tsv'), check.names = F)
act_genes <- act_genes %>% 
  filter(`Drugs (for therapeutic implications only)` != '') %>%
  mutate(Alterations = strsplit(as.character(Alterations), ", ")) %>% 
  unnest(Alterations)
act_genes$Alterations <- gsub('.*Fusion$|^Fusions$', 'Fusion', act_genes$Alterations)

# merge actionable genes with merged annotated oncokb output (no result)
final_file <- merged_file %>%
  inner_join(act_genes, by = c('GENE' = 'Gene', 'ALTERATION' = 'Alterations')) %>%
  unique()
write.table(final_file, file = file.path(output_dir, 'oncokb_merged_annotated_actgenes.txt'), quote = F, sep = "\t", row.names = F)

