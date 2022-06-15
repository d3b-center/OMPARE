# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "mutationally_similar_analysis")
output_dir <- file.path(patient_dir, "output", "mutationally_similar_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "kaplan_meier.R"))

# load inputs
nn_table <- file.path(patient_dir, "output", "mut_distance_calc", "pediatric", "nn_table.rds")
nn_table <- readRDS(nn_table)
surv_data <-  file.path(patient_dir, "output", "mut_distance_calc", "pediatric", "surv_data.rds")
surv_data <- readRDS(surv_data)

# survival analysis
fname <- file.path(output_dir, "kaplan_meier_pediatric.pdf")
kaplan_meier_pediatric <- kaplan_meier(nn_table = nn_table, 
                                       surv_data = surv_data)

# save plot
pdf(fname)
print(kaplan_meier_pediatric, newpage = FALSE)
dev.off()

