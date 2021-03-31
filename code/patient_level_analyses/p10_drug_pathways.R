# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'drug_pathways.R'))

# input files
output_dir <- file.path(topDir, "output")

# call function
fname <- file.path(topDir, "output", "drug_pathways_barplot.rds")
drug_pathways_barplot <- drug_pathways(pnoc008_patient = sampleInfo$subjectID, output_dir = output_dir)
saveRDS(drug_pathways_barplot, file = fname)
