# merge annotations from the preceding steps with actionable genes information.
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# arguments
option_list <- list(
  make_option(c("--output_dir"), type = "character",
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
output_dir <- opt$output_dir

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# read actionable genes
act_genes <- read.delim(file.path(data_dir, "oncokb", "oncokb_biomarker_drug_associations.tsv"), check.names = F)
act_genes <- act_genes %>% 
  filter(`Drugs (for therapeutic implications only)` != '') %>%
  mutate(Alterations = strsplit(as.character(Alterations), ", ")) %>% 
  unnest(Alterations)
act_genes$Alterations <- gsub('.*Fusion$|^Fusions$', 'Fusion', act_genes$Alterations)

# read annotated merged output
merged_file <- file.path(output_dir, paste0('oncokb_merged_consensus_annotated.txt'))
merged_file <- read.delim(merged_file, stringsAsFactors = F)

# merge actionable genes with merged annotated oncokb output (no result)
output_file <- file.path(output_dir, paste0('oncokb_merged_consensus_annotated_actgenes.txt'))
final_file <- merged_file %>%
  inner_join(act_genes, by = c('GENE' = 'Gene', 'ALTERATION' = 'Alterations')) %>%
  unique()

# rename columns
final_file <- final_file %>%
  rename("LEVEL" = "Level",
         "CANCER_TYPES" = "Cancer Types",
         "DRUGS" = "Drugs (for therapeutic implications only)")
write.table(final_file, file = output_file, quote = F, sep = "\t", row.names = F)

