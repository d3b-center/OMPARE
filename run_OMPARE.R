# Author: Komal S. Rathi
# Date: 04/13/2020
# Function: Generate patient report

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rmarkdown))

option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("-s", "--sourcedir"), type = "character", 
              default = NULL,
              help = "Source directory with all files"),
  make_option(c("-c", "--clin_file"), type = "character",
              default = NULL,
              help = "Manifest file (.xlsx)")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
clinical_sheet <- opt$clin_file
sourceDir <- opt$sourcedir

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# set variables
topDir <- file.path(root_dir, 'results', patient)
set_title <- paste0(patient,' Patient Report')
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")

# 1. Create Project directory (if sourcedir param is provided)
if(!is.null(sourceDir)){
  print("Create Project Directory...")
  create.dirs <- file.path(code_dir, 'create_project_dir.R')
  cmd1 <- paste('Rscript', create.dirs, '-s', sourceDir, '-d', topDir)
  print(cmd1)
  system(cmd1)
} else {
  print("Project Directory found...")
}

# 2. Create clinical file  (if clin_file param is provided)
if(!is.null(clinical_sheet)){
  print("Create Clinical file...")
  create.clinfile <- file.path(code_dir, 'create_clinfile.R')
  cmd2 <- paste('Rscript', create.clinfile, '-s', clinical_sheet, '-p', patient, '-d', topDir)
  print(cmd2)
  system(cmd2)
} else {
  print("Clinical file present...")
}

# 3. Update PNOC008 matrices for each new patient
print("Update PNOC008 data matrices...")
pnoc.format <- file.path(patient_level_analyses, 'pnoc_format.R')
cmd3 <- paste('Rscript', pnoc.format, '--clin_file', clinical_sheet)
print(cmd3)
system(cmd3)

# 4. Update GSEA enrichment for each new patient
print("Update PNOC008 GSEA summary...")
gsea.enrichment <- file.path(patient_level_analyses, 'gsea_enrichment.R')
cmd4 <- paste('Rscript', gsea.enrichment, '-p', patient)
print(cmd4)
system(cmd4)

# 5. Generate genes and pathway enrichment output
print("Generate output for RNA-seq enrichment...")
rnaseq_enrichment <- file.path(patient_level_analyses, 'enrichment_output.R')
cmd5 <- paste('Rscript', rnaseq_enrichment, '-i', topDir, '-o', paste0(patient, '_summary'), '-t text')
print(cmd5)
system(cmd5)

# 6. Run html reports
# fusion_method can be either arriba, star, both or not specified
print("Run reports...")
if(dir.exists(topDir)){
  for(i in 1:length(callers)) {
    output_dir <- file.path(topDir, 'Reports')
    output_file <- paste0(patient, '_', callers[i], '.html')
    input_file <- file.path(root_dir, 'OMPARE.Rmd')
    rmarkdown::render(input = input_file,
                      params = list(topDir = topDir,
                                    fusion_method = 'arriba',
                                    set_title = set_title,
                                    snv_caller = callers[i],
                                    tmb = 77.46), 
                      output_dir = output_dir, 
                      intermediates_dir = output_dir,
                      output_file = output_file)
  }
}
