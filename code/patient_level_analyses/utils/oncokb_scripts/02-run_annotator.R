# run annotator to annotate mutations, fusions and cnv
library(optparse)

# arguments
option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("-s", "--snv_caller"), type = "character",
              help = "SNV caller pattern: lancet, vardict, consensus, strelka2, mutect2")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
snv_caller <- opt$snv_caller

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# oncokb binary
readRenviron('~/.Renviron')
oncokb_token <- Sys.getenv('oncokb_token')
oncokb_tool <- '~/tools/oncokb-annotator/'
maf_annotator <- file.path(oncokb_tool, 'MafAnnotator.py')
cnv_annotator <- file.path(oncokb_tool, 'CnaAnnotator.py')
fusion_annotator <- file.path(oncokb_tool, 'FusionAnnotator.py')

# patient <- 'PNOC008-28'
# output directory
output_dir <- file.path(results_dir, patient, 'output')

# mutations
maf_dir <- file.path(results_dir, patient, 'simple-variants')
maf_files <- list.files(path = maf_dir, pattern = "*.maf", full.names = T)
maf_file <- grep(snv_caller, maf_files, value = TRUE)
maf_out <- file.path(output_dir, paste0('oncokb_', snv_caller, '_annotated.txt'))
command <- paste('python', maf_annotator, '-i', maf_file, '-o', maf_out, '-b', oncokb_token, '-q hgvsp_short')
system(command)

# cnv
cnv_out <- file.path(output_dir, 'oncokb_cnv_annotated.txt')
if(!file.exists(cnv_out)){
  cnv_file <- file.path(output_dir, 'oncokb_cnv.txt')
  command <- paste('python', cnv_annotator, '-i', cnv_file, '-o', cnv_out, '-b', oncokb_token)
  system(command)
}

# fusion 
fusion_out <- file.path(output_dir, 'oncokb_fusion_annotated.txt')
if(!file.exists(fusion_out)){
  fusion_file <- file.path(output_dir, 'oncokb_fusion.txt')
  command <- paste('python', fusion_annotator, '-i', fusion_file, '-o', fusion_out, '-b', oncokb_token)
  system(command)
}

