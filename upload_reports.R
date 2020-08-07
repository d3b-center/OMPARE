# Author: Komal S. Rathi
# Date: 07/12/2020
# Function: Upload patient reports to data delivery project

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient Number (1, 2...)"),
  make_option(c("-w", "--workdir"), type = "character",
              help = "OMPARE working directory")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
p <- opt$patient
workdir <- opt$workdir

# set variables
setwd(workdir) # This should be the OMPARE directory
print(paste0("Working directory:", getwd()))
patient <- paste0('PNOC008-', p)
topDir <- file.path(getwd(), 'data', patient)

# set variables for upload commands
readRenviron("~/.Renviron")
cav <- Sys.getenv('CAV') # path to cavatica-uploader.sh
auth <- Sys.getenv('AUTH_TOKEN') # authentication token
project <- 'cavatica/sd-8y99qzjj' # project id
report.folder <- paste0(patient,'/reports') # destination folder
summary <- file.path(topDir, 'Summary', paste0(patient, '_summary.xlsx')) # summary excel file
reports <- file.path(topDir, 'Reports', paste0(patient, '_*.html')) # html reports
umap <- file.path(topDir, 'Summary', paste0('*_umap_output.rds')) # umap output
cmd1 <- paste(cav, '-t', auth, '-p', project, '-f', report.folder, summary, sep = " ")
print(cmd1)
system(cmd1)
cmd2 <- paste(cav, '-t', auth, '-p', project, '-f', report.folder, reports, sep = " ")
print(cmd2)
system(cmd2)
cmd3 <- paste(cav, '-t', auth, '-p', project, '-f', report.folder, umap, sep = " ")
print(cmd3)
system(cmd3)
