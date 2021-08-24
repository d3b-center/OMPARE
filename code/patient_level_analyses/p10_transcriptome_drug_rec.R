# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'transcriptome_drug_rec.R'))

# input files
genes_up <- read.delim(file.path(topDir, "output", paste0(sampleInfo$subjectID, '_summary_DE_Genes_Up.txt')))

# call function
fname <- file.path(topDir, "output", "transcriptome_drug_rec.rds")
transcriptome_drug_rec_output <- transcriptome_drug_rec(diffexpr_genes = genes_up)
saveRDS(transcriptome_drug_rec_output, file = fname)
