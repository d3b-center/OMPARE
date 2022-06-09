# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "survival_analysis")
output_dir <- file.path(patient_dir, "output", "survival_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "kaplan_meier.R"))

# load inputs
nn_table <- file.path(patient_dir, "output", "transcriptomically_similar_analysis", "pediatric_all_nn_table.rds")
nn_table <- readRDS(nn_table)
pediatric_patient_clinical <-  file.path(patient_dir, "output", "transcriptomically_similar_analysis", "pediatric_all_patient_combined_clinical_input.rds")
pediatric_patient_clinical <- readRDS(pediatric_patient_clinical)

# survival analysis
fname <- file.path(output_dir, "kaplan_meier_pediatric.pdf")
kaplan_meier_pediatric <- kaplan_meier(nn_table = nn_table, 
                                       surv_data = pediatric_patient_clinical)

# save plot
pdf(fname)
print(kaplan_meier_pediatric, newpage = FALSE)
dev.off()

