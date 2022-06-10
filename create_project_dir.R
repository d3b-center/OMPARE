# Author: Komal S. Rathi
# Date: 02/28/2020
# Function: Code to create and organize directory structure
# This will be called from within run_OMPARE.R
# E.g.: Rscript create_project <path to directory with all files> <path to project directory>

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--source_dir"), type = "character",
              help = "Source directory with all files"),
  make_option(c("--patient_dir"), type = "character",
              help = "Destination directory. Should be /path/to/OMPARE/results/PNOC008-22/ for Patient 22")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
source_dir <- opt$source_dir
patient_dir <- opt$patient_dir

# specify destination directories
clinicaldir <- file.path(patient_dir, 'clinical')
cnvdir <- file.path(patient_dir, 'copy-number-variations')
exprdir <- file.path(patient_dir, 'gene-expressions')
fusionsdir <- file.path(patient_dir, 'gene-fusions')
mutdir <- file.path(patient_dir, 'simple-variants')
reports <- file.path(patient_dir, 'reports')
output <- file.path(patient_dir, 'output')

# create directories to move input files
dir.create(clinicaldir, showWarnings = F, recursive = T)
dir.create(cnvdir, showWarnings = F, recursive = T)
dir.create(exprdir, showWarnings = F, recursive = T)
dir.create(fusionsdir, showWarnings = F, recursive = T)
dir.create(mutdir, showWarnings = F, recursive = T)

# create directories to save output
dir.create(output, showWarnings = F, recursive = T)
dir.create(reports, showWarnings = F, recursive = T)

# organize data
# copy number
cmd <- file.path(source_dir, '*.{gainloss.txt,info.txt,diagram.pdf}')
cmd <- paste('mv', cmd, cnvdir)
system(cmd)

# clinical file is to be obtained from KF data tracker

# expression
cmd <- file.path(source_dir, '*rsem*')
cmd <- paste('mv', cmd, exprdir)
system(cmd)

# fusions
cmd <- file.path(source_dir, '*fusion*')
cmd <- paste('mv', cmd, fusionsdir)
system(cmd)

# mutations - somatic/germline
cmd <- file.path(source_dir, '*.{maf,hg38_multianno.txt.gz}')
cmd <- paste('mv', cmd, mutdir)
system(cmd)



