# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "survival_analysis")
output_dir <- file.path(patient_dir, "output", "survival_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "kaplan_meier.R"))

# survival analysis
fname <- file.path(output_dir, "kaplan_meier_pediatric.rds")
if(!file.exists(fname)){
  kaplan_meier_pediatric <- kaplan_meier(all_cor = pbta_hgat_pnoc008_nn_table, 
                                         surv_data = pbta_survival)
  
  # save output
  saveRDS(kaplan_meier_pediatric, file = file.path(fname))
} else {
  kaplan_meier_pediatric <- readRDS(fname)
}
