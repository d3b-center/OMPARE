# run annotator to annotate mutations, fusions and cnv
suppressPackageStartupMessages({
  library(optparse)
})

# arguments
option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("--snv_caller"), type = "character",
              help = "SNV caller pattern: lancet, vardict, consensus, strelka2, mutect2"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
snv_caller <- opt$snv_caller
output_dir <- opt$output_dir

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
patient_dir <- file.path(root_dir, "results", patient)

# oncokb binary
readRenviron("~/.Renviron")
oncokb_token <- Sys.getenv("oncokb_token")
oncokb_tool <- "~/tools/oncokb-annotator/"
maf_annotator <- file.path(oncokb_tool, "MafAnnotator.py")
cnv_annotator <- file.path(oncokb_tool, "CnaAnnotator.py")
fusion_annotator <- file.path(oncokb_tool, "FusionAnnotator.py")

# mutations
maf_dir <- file.path(patient_dir, "simple-variants")
maf_files <- list.files(path = maf_dir, pattern = "*.maf", full.names = T)
maf_file <- grep(snv_caller, maf_files, value = TRUE)
maf_out <- file.path(output_dir, paste0('oncokb_', snv_caller, '_annotated.txt'))
command <- paste('/usr/bin/python', maf_annotator, '-i', maf_file, '-o', maf_out, '-b', oncokb_token, '-q hgvsp_short')
system(command)

# cnv
cnv_out <- file.path(output_dir, "oncokb_cnv_annotated.txt")
if(!file.exists(cnv_out)){
  cnv_file <- file.path(output_dir, "oncokb_cnv.txt")
  command <- paste('python', cnv_annotator, '-i', cnv_file, '-o', cnv_out, '-b', oncokb_token)
  system(command)
}

# fusion 
fusion_out <- file.path(output_dir, "oncokb_fusion_annotated.txt")
if(!file.exists(fusion_out)){
  fusion_file <- file.path(output_dir, "oncokb_fusion.txt")
  command <- paste('python', fusion_annotator, '-i', fusion_file, '-o', fusion_out, '-b', oncokb_token)
  system(command)
}

