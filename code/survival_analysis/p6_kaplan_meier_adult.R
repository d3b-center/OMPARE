# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "survival_analysis")
output_dir <- file.path(patient_dir, "output", "survival_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "kaplan_meier.R"))

# load inputs
nn_table <- file.path(patient_dir, "output", "transcriptomically_similar_analysis", "adult_nn_table.rds")
nn_table <- readRDS(nn_table)
adult_patient_clinical <-  file.path(patient_dir, "output", "transcriptomically_similar_analysis", "adult_patient_combined_clinical_input.rds")
adult_patient_clinical <- readRDS(adult_patient_clinical)

# survival analysis
fname <- file.path(output_dir, "kaplan_meier_adult.pdf")
kaplan_meier_adult <- kaplan_meier(nn_table = nn_table, 
                                   surv_data = adult_patient_clinical)

# save plot
pdf(fname)
print(kaplan_meier_adult, newpage = FALSE)
dev.off()
