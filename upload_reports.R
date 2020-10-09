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
topDir <- file.path(getwd(), 'results', patient)

# set variables for upload commands
readRenviron("~/.Renviron")
cav <- Sys.getenv('CAV') # path to cavatica-uploader.sh
auth <- Sys.getenv('AUTH_TOKEN') # authentication token
project <- file.path('cavatica', 'sd-8y99qzjj') # project id

# destination folder
dest.folder <- file.path(patient) 

# source folders
output <- file.path(topDir, "output") # all output
reports <- file.path(topDir, "reports") # all reports

cmd1 <- paste(cav, '-t', auth, '-p', project, '-f', dest.folder, '-pf', output, sep = " ")
print(cmd1)
# system(cmd1)
cmd2 <- paste(cav, '-t', auth, '-p', project, '-f', dest.folder, '-pf', reports, sep = " ")
print(cmd2)
# system(cmd2)
