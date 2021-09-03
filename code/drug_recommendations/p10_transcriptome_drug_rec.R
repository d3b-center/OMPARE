# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "drug_recommendations")
output_dir <- file.path(patient_dir, "output", "drug_recommendations")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "transcriptome_drug_rec.R"))

# input files
genes_up <- read.delim(file.path(patient_dir, "output", "rnaseq_analysis", paste0(patient, '_summary_DE_Genes_Up.txt')))

# transcriptomic drug recommendations
fname <- file.path(output_dir, "transcriptome_drug_rec.rds")
if(!file.exists(fname)){
  # full output save as tsv
  full_fname <- file.path(output_dir, paste0(patient, "_CHEMBL_drug-gene.tsv"))
  transcriptome_drug_rec_output <- transcriptome_drug_rec(diffexpr_genes = genes_up)
  write.table(x = transcriptome_drug_rec_output, file = full_fname, quote = F, sep = "\t", row.names = F)
  
  # subset columns for report
  transcriptome_drug_rec_output <- transcriptome_drug_rec_output %>%
    dplyr::select(Drug, Gene, Source, Comparison, logFC, MOA) %>%
    unique()
  saveRDS(transcriptome_drug_rec_output, file = fname)
} else {
  transcriptome_drug_rec_output <- readRDS(fname)
}
