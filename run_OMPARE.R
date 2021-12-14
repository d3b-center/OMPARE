# Author: Komal S. Rathi
# Date: 04/13/2020
# Function: Generate patient report

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier (PNOC008-22, C3342894...)"),
  make_option(c("--source_dir"), type = "character", 
              default = NULL,
              help = "Source directory with all files"),
  make_option(c("--clin_file"), type = "character",
              default = NULL,
              help = "Manifest file (.xlsx)"),
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
source_dir <- opt$source_dir
sync_data <- opt$sync_data
upload_reports <- opt$upload_reports
study <- opt$study

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# set variables
patient_dir <- file.path(root_dir, 'results', patient)
set_title <- paste0(patient,' Patient Report')
callers <- c("lancet", "mutect2", "strelka2", "vardict", "consensus", "all")

# 1. Create Project directory (if source_dir param is provided)
if(!is.null(source_dir)){
  print("Create Project Directory...")
  create.dirs <- file.path(root_dir, "code", 'create_project_dir.R')
  cmd <- paste('Rscript', create.dirs, 
                '--source_dir', source_dir, 
                '--patient_dir', patient_dir)
  print(cmd)
  system(cmd)
}

# 2. Create clinical file  (if clin_file param is provided)
if(!is.null(clinical_sheet)){
  print("Create Clinical file...")
  create.clinfile <- file.path(root_dir, "code", 'create_clinfile.R')
  cmd <- paste('Rscript', create.clinfile, 
                '--sheet', clinical_sheet, 
                '--patient', patient, 
                '--dir', patient_dir)
  print(cmd)
  system(cmd)
}

# 3. Update PNOC008 matrices for each new patient
print("Update PNOC008 data matrices...")
update_pnoc008_matrices <- file.path(root_dir, "code", 'update_pnoc008_matrices.R')
cmd <- paste('Rscript', update_pnoc008_matrices, 
              '--clin_file', clinical_sheet)
print(cmd)
system(cmd)

# 4. Run html reports
print("Run reports...")
for(i in 1:length(callers)) {
  output_dir <- file.path(patient_dir, 'reports')
  output_file <- paste0(patient, '_', callers[i], '.html')
  input_file <- file.path(root_dir, 'OMPARE.Rmd')
  rmarkdown::render(input = input_file,
                    params = list(patient_dir = patient_dir,
                                  set_title = set_title,
                                  snv_caller = callers[i]), 
                    output_dir = output_dir, 
                    intermediates_dir = output_dir,
                    output_file = output_file)
}

# 5. Sync to s3 (needs VPN)
# exclude chembl folder as it is huge ~20GB and connection breaks before it is uploaded.
if(sync_data){
  print("Sync data back to s3...")
  cmd <- paste("aws s3 --profile saml sync", data_dir, "s3://d3b-bix-dev-data-bucket/PNOC008/data --exclude 'chembl/*'")
  print(cmd)
  system(cmd)
}

# 6. Upload reports to cavatica
if(upload_reports){
  print("Upload reports to cavatica...")
  cmd <- paste("Rscript upload_reports.R", 
                "--patient", patient, 
                "--study", study)
  print(cmd)
  system(cmd)
}
