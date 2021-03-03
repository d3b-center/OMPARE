# merge annotations from the preceding steps with actionable genes information.
library(tidyverse)
library(optparse)

# arguments
option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("-s", "--snv_caller"), type = "character",
              help = "SNV caller pattern: lancet, vardict, consensus, strelka2, mutect2 and all")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
snv_caller <- opt$snv_caller

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
output_dir <- file.path(results_dir, patient, 'output')

# patient <- 'PNOC008-28'

# read annotated merged output
merged_file <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated.txt'))
merged_file <- read.delim(merged_file, stringsAsFactors = F)

# read actionable genes
act_genes <- read.delim(file.path(ref_dir, 'oncokb', 'oncokb_biomarker_drug_associations.tsv'), check.names = F)
act_genes <- act_genes %>% 
  filter(`Drugs (for therapeutic implications only)` != '') %>%
  mutate(Alterations = strsplit(as.character(Alterations), ", ")) %>% 
  unnest(Alterations)
act_genes$Alterations <- gsub('.*Fusion$|^Fusions$', 'Fusion', act_genes$Alterations)

# merge actionable genes with merged annotated oncokb output (no result)
output_file <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated_actgenes.txt'))
final_file <- merged_file %>%
  inner_join(act_genes, by = c('GENE' = 'Gene', 'ALTERATION' = 'Alterations')) %>%
  unique()

# rename columns
final_file <- final_file %>%
  rename("LEVEL" = "Level",
         "CANCER_TYPES" = "Cancer Types",
         "DRUGS" = "Drugs (for therapeutic implications only)")
write.table(final_file, file = output_file, quote = F, sep = "\t", row.names = F)

