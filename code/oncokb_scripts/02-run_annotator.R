# run annotator to annotate mutations, fusions and cnv
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

patient <- 'PNOC008-28'
output_dir <- file.path(results_dir, patient, 'output')
readRenviron('~/.Renviron')
oncokb_tool <- '~/tools/oncokb-annotator/'
oncokb_token <- Sys.getenv('oncokb_token')
maf_annotator <- file.path(oncokb_tool, 'MafAnnotator.py')
cnv_annotator <- file.path(oncokb_tool, 'CnaAnnotator.py')
fusion_annotator <- file.path(oncokb_tool, 'FusionAnnotator.py')

# mutations
maf_dir <- file.path(results_dir, patient, 'simple-variants')
maf_files <- list.files(path = maf_dir, pattern = "*.maf", full.names = T)
for(i in 1:length(maf_files)){
  print(i)
  output_filename <- gsub('_somatic.vep.maf', '', maf_files[i])
  output_filename <- gsub('.*[.]', '', output_filename)
  output_filename <- file.path(output_dir, paste0('oncokb_', output_filename, '_annotated.txt'))
  command <- paste('python', maf_annotator, '-i', maf_files[i], '-o', output_filename, '-b', oncokb_token, '-q hgvsp_short')
  system(command)
}

# cnv
cnv_file <- file.path(output_dir, 'oncokb_cnv.txt')
output_filename <- file.path(output_dir, 'oncokb_cnv_annotated.txt')
command <- paste('python', cnv_annotator, '-i', cnv_file, '-o', output_filename, '-b', oncokb_token)
system(command)

# fusion 
fusion_file <- file.path(output_dir, 'oncokb_fusion.txt')
output_filename <- file.path(output_dir, 'oncokb_fusion_annotated.txt')
command <- paste('python', fusion_annotator, '-i', fusion_file, '-o', output_filename, '-b', oncokb_token)
system(command)

