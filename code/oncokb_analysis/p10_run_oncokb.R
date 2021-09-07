# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "oncokb_analysis")
output_dir <- file.path(patient_dir, "output", "oncokb_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# format for oncokb
script <- file.path(module_dir, '01-format_for_oncokb.R')
cmd <- paste('Rscript', script, 
             '--patient', patient, 
             '--output_dir', output_dir)
system(cmd)

# run annotator
script <- file.path(module_dir, '02-run_annotator.R')
if(snv_caller != "all"){
  cmd <- paste('Rscript', script, 
               '--patient', patient, 
               '--snv_caller', snv_caller,
               '--output_dir', output_dir)
  system(cmd)
}

# merge oncokb annotated outputs
script <- file.path(module_dir, '03-merge_oncokb_annotated.R')
cmd <- paste('Rscript', script, 
             '--patient', patient, 
             '--snv_caller', snv_caller,
             '--output_dir', output_dir)
system(cmd)

# merge oncokb annotated outputs with actionable genes
script <- file.path(module_dir, '04-merge_oncokb_annotated_actionable_genes.R')
cmd <- paste('Rscript', script, 
             '--patient', patient, 
             '--snv_caller', snv_caller,
             '--output_dir', output_dir)
system(cmd)

# read output generated from running the above scripts
fname <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated_actgenes.txt'))
oncokb_output <- read.delim(fname, stringsAsFactors = F)

