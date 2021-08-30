# Author: Komal S. Rathi
# Function: Upload patient output and reports to data delivery project

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient Identifier (PNOC008-22, etc...)"),
  make_option(c("--workdir"), type = "character",
              help = "OMPARE working directory"),
  make_option(c("--study"), type = "character", 
              default = "PNOC008",
              help = "PNOC008 or CBTN")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
patient <- opt$patient
workdir <- opt$workdir
study <- opt$study

# set variables
setwd(workdir) # This should be the OMPARE directory
print(paste0("Working directory:", getwd()))
topDir <- file.path(getwd(), 'results', patient)

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
dest.folder <- file.path(patient) 

# source folders
output <- file.path(topDir, "output") # all output
reports <- file.path(topDir, "reports") # all reports
cemitools_reports <- file.path(topDir, "CEMITools") # cemitools output

# upload output
cmd1 <- paste(cav, '-t', auth, '-p', project, '-f', dest.folder, '-pf', output, sep = " ")
print(cmd1)
system(cmd1)

# upload reports
cmd2 <- paste(cav, '-t', auth, '-p', project, '-f', dest.folder, '-pf', reports, sep = " ")
print(cmd2)
system(cmd2)

# upload cemitools reports and output
cmd3 <- paste(cav, '-t', auth, '-p', project, '-f', dest.folder, '-pf', cemitools_reports, sep = " ")
print(cmd3)
system(cmd3)

