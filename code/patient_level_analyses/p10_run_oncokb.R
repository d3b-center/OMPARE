# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
oncokb_scripts <- file.path(patient_level_analyses_utils, 'oncokb_scripts')

# format for oncokb
script <- file.path(oncokb_scripts, '01-format_for_oncokb.R')
cmd <- paste('Rscript', script, '--patient', sampleInfo$subjectID)
system(cmd)

# run annotator
script <- file.path(oncokb_scripts, '02-run_annotator.R')
if(snv_caller != "all"){
  cmd <- paste('Rscript', script, '--patient', sampleInfo$subjectID, '--snv_caller', snv_caller)
  system(cmd)
}

# merge oncokb annotated outputs
script <- file.path(oncokb_scripts, '03-merge_oncokb_annotated.R')
cmd <- paste('Rscript', script, '--patient', sampleInfo$subjectID, '--snv_caller', snv_caller)
system(cmd)

# merge oncokb annotated outputs with actionable genes
script <- file.path(oncokb_scripts, '04-merge_oncokb_annotated_actionable_genes.R')
cmd <- paste('Rscript', script, '--patient', sampleInfo$subjectID, '--snv_caller', snv_caller)
system(cmd)
