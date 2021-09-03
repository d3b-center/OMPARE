# Author: Komal S. Rathi
# Function: Upload patient output and reports to data delivery project

suppressPackageStartupMessages({
	library(optparse)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient Identifier (PNOC008-22, etc...)"),
  make_option(c("--study"), type = "character", 
              default = "PNOC008",
              help = "PNOC008 or CBTN")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
study <- opt$study

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
patient_dir <- file.path(root_dir, 'results', patient)

# set variables for upload commands
readRenviron("~/.Renviron")
cav <- Sys.getenv('CAV') # path to cavatica-uploader.sh
auth <- Sys.getenv('AUTH_TOKEN') # authentication token
if(study == "PNOC008"){
  project <- file.path('cavatica', 'sd-8y99qzjj') # pnoc008 project id
} else {
  project <- file.path('d3b-bixu', 'sd-bhjxbdqk-realtime') # cbtn project id
}

# destination folder
dest_folder <- file.path(patient) 

# source folders
output <- file.path(patient_dir, "output") # all output
reports <- file.path(patient_dir, "reports") # all reports

# upload output
cmd <- paste(cav, 
	'-t', auth, 
	'-p', project, '-f', dest_folder, 
	'-pf', output)
print(cmd)
system(cmd)

# upload reports
cmd <- paste(cav, 
	'-t', auth, 
	'-p', project, 
	'-f', dest_folder, 
	'-pf', reports)
print(cmd)
system(cmd)

