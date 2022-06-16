# Author: Komal S. Rathi
# Function: generalized patient report

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "patient identifier (cohort participant id)"),
  make_option(c("--patient_cancer_type"), type = "character",
              help = "patient cancer type"),
  make_option(c("--mut_only"), type = "logical", default = FALSE,
              help = "DNA only report (TRUE or FALSE)"),
  make_option(c("--rnaseq_only"), type = "logical", default = FALSE,
              help = "RNA only report (TRUE or FALSE)"),
  make_option(c("--source_dir"), type = "character", default = NULL,
              help = "source directory with files downloaded from data delivery project"),
  make_option(c("--sync_data"), type = "logical", default = FALSE,
              help = "sync reference data to s3 (TRUE or FALSE)"),
  make_option(c("--upload_reports"), type = "logical", default = FALSE,
              help = "upload reports to cavatica (TRUE or FALSE)"),
  make_option(c("--study"), type = "character",
              help = "study identifier i.e. PNOC008, CBTN, etc required to upload output to relevant cavatica project")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
patient_cancer_type <- opt$patient_cancer_type
mut_only <- opt$mut_only
rnaseq_only <- opt$rnaseq_only
source_dir <- opt$source_dir
sync_data <- opt$sync_data
upload_reports <- opt$upload_reports
study <- opt$study

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# set variables
patient_dir <- file.path(root_dir, 'results', patient)
set_title <- paste0(patient, ' Patient Report')

# set driver script and markdown file
if(mut_only){
  driver_script <- file.path(root_dir, "code", "driver_mut_only.R")
  input_file <- file.path(root_dir, 'patient_report_mut_only.Rmd')
} else if(rnaseq_only){
  input_file <- file.path(root_dir, 'patient_report_rnaseq_only.Rmd')
  driver_script <- file.path(root_dir, "code", "driver_rnaseq_only.R")
} else {
  input_file <- file.path(root_dir, 'patient_report.Rmd')
  driver_script <- file.path(root_dir, "code", "driver.R")
}

# organize downloaded files into patient directory
if(!is.null(source_dir)){
  print("Create Project Directory...")
  create.dirs <- file.path(root_dir, 'create_project_dir.R')
  cmd <- paste('Rscript', create.dirs, 
               '--source_dir', source_dir, 
               '--patient_dir', patient_dir)
  print(cmd)
  system(cmd)
}

# output file
if(rnaseq_only){
  output_file <- paste0(patient, ".html")
} else {
  output_file <- paste0(patient, "_consensus.html")
}

# call driver script to generate all outputs
source(driver_script)
run_driver(patient = patient, 
           patient_cancer_type = patient_cancer_type, 
           patient_dir = patient_dir)

# generate html reports
print("Run reports...")
output_dir <- file.path(patient_dir, 'reports')
rmarkdown::render(input = input_file,
                  params = list(patient = patient,
                                patient_dir = patient_dir,
                                set_title = set_title), 
                  output_dir = output_dir, 
                  output_file = output_file)

# sync to s3 (needs VPN)
# exclude OpenPedCan-analysis and chembl folder as it is huge ~20GB and connection breaks before it is uploaded.
if(sync_data){
  print("Sync data back to s3...")
  cmd <- paste("aws s3 --profile Mgmt-Console-Dev-chopd3bprod@684194535433 sync", data_dir, 
               "s3://d3b-bix-dev-data-bucket/d3b-patient-report-analysis/data --exclude 'chembl/*' --exclude 'OpenPedCan-analysis/*'")
  print(cmd)
  system(cmd)
}

# upload reports to cavatica
if(upload_reports){
  print("Upload reports to cavatica...")
  cmd <- paste("Rscript upload_reports.R", 
               "--patient_dir", patient_dir, 
               "--study", study)
  print(cmd)
  system(cmd)
}
