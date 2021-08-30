# Author: Komal S. Rathi
# Date: 04/13/2020
# Function: Generate patient report

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rmarkdown))

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("--sourcedir"), type = "character", 
              default = NULL,
              help = "Source directory with all files"),
  make_option(c("--clin_file"), type = "character",
              default = NULL,
              help = "Manifest file (.xlsx)"),
  make_option(c("--update_pbta"), type = "logical",
              default = NULL,
              help = "Update PBTA adapt file (TRUE or FALSE)"),
  make_option(c("--sync_data"), type = "logical",
              default = NULL,
              help = "Sync reference data to s3 (TRUE or FALSE)"),
  make_option(c("--upload_reports"), type = "logical",
              default = NULL,
              help = "Upload reports to cavatica (TRUE or FALSE)"),
  make_option(c("--study"), type = "logical",
              default = "PNOC008",
              help = "Study ID (PNOC008 or CBTN)")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
clinical_sheet <- opt$clin_file
sourceDir <- opt$sourcedir
update_pbta <- opt$update_pbta
sync_data <- opt$sync_data
upload_reports <- opt$upload_reports
study <- opt$study

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
  cmd1 <- paste('Rscript', create.dirs, '--sourcedir', sourceDir, '--destdir', topDir)
  print(cmd1)
  system(cmd1)
} else {
  print("Project Directory found...")
}

# 2. Create clinical file  (if clin_file param is provided)
if(!is.null(clinical_sheet)){
  print("Create Clinical file...")
  create.clinfile <- file.path(code_dir, 'create_clinfile.R')
  cmd2 <- paste('Rscript', create.clinfile, '--sheet', clinical_sheet, '--patient', patient, '--dir', topDir)
  print(cmd2)
  system(cmd2)
} else {
  print("Clinical file present...")
}

# 3. update histologies file from adapt (needs VPN)
if(update_pbta){
  print("Update PBTA histologies file...")
  update_pbta <- file.path(code_dir, 'update_pbta.R')
  cmd3 <- paste('Rscript', update_pbta)
  print(cmd3)
  system(cmd3)
}

# 4. Update PNOC008 matrices for each new patient
print("Update PNOC008 data matrices...")
pnoc.format <- file.path(patient_level_analyses, 'pnoc_format.R')
cmd4 <- paste('Rscript', pnoc.format, '--clin_file', clinical_sheet)
print(cmd4)
system(cmd4)

# 5. Update GSEA enrichment for each new patient
print("Update PNOC008 GSEA summary...")
gsea.enrichment <- file.path(patient_level_analyses, 'gsea_enrichment.R')
cmd5 <- paste('Rscript', gsea.enrichment, '--patient', patient)
print(cmd5)
system(cmd5)

# 6. Generate genes and pathway enrichment output
print("Generate output for RNA-seq enrichment...")
rnaseq_enrichment <- file.path(patient_level_analyses, 'enrichment_output.R')
cmd6 <- paste('Rscript', rnaseq_enrichment, '--input', topDir, '--output', paste0(patient, '_summary'), '--type text')
print(cmd6)
system(cmd6)

# 7. Run html reports
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
                                    snv_caller = callers[i]), 
                      output_dir = output_dir, 
                      intermediates_dir = output_dir,
                      output_file = output_file)
  }
}

# 8. Sync to s3 (needs VPN)
# exclude chembl folder as it is huge ~20GB and connection breaks before it is uploaded.
if(sync_data){
  print("Sync data back to s3...")
  cmd8 <- paste("aws s3 --profile saml sync", ref_dir, "s3://d3b-bix-dev-data-bucket/PNOC008/reference --exclude 'chembl/*'")
  print(cmd8)
  system(cmd8)
}

# 9. Upload reports to cavatica
if(upload_reports){
  print("Upload reports to cavatica...")
  cmd9 <- paste("Rscript upload_reports.R --patient", patient, "--study", study)
  print(cmd9)
  system(cmd9)
}
