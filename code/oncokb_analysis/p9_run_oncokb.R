# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "oncokb_analysis")
output_dir <- file.path(patient_dir, "output", "oncokb_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# master genomics files
master_genomics_dir <- file.path(data_dir, "master_genomics")

# patient id
sample_info <- read.delim(file.path(patient_dir, "output", "sample_info.tsv"))
patient <- sample_info %>%
  pull(cohort_participant_id) %>%
  unique()

# format for oncokb
script <- file.path(module_dir, '01-format_for_oncokb.R')
cmd <- paste('Rscript', script, 
             '--patient_of_interest', patient, 
             '--master_genomics_dir', master_genomics_dir,
             '--output_dir', output_dir)
system(cmd)

# run annotator
script <- file.path(module_dir, '02-run_annotator.R')
cmd <- paste('Rscript', script, 
             '--patient_of_interest', patient, 
             '--output_dir', output_dir)
system(cmd)

# merge oncokb annotated outputs
script <- file.path(module_dir, '03-merge_oncokb_annotated.R')
cmd <- paste('Rscript', script, 
             '--output_dir', output_dir)
system(cmd)

# merge oncokb annotated outputs with actionable genes
script <- file.path(module_dir, '04-merge_oncokb_annotated_actionable_genes.R')
cmd <- paste('Rscript', script, 
             '--output_dir', output_dir)
system(cmd)

# read output generated from running the above scripts
fname <- file.path(output_dir, paste0('oncokb_merged_annotated_actgenes.txt'))
oncokb_output <- read.delim(fname, stringsAsFactors = F)
