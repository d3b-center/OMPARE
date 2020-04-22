# Author: Komal S. Rathi
# Date: 02/28/2020
# Function: Code to create and organize directory structure
# This will be called from within run_OMPARE.R
# E.g.: Rscript create_project <path to directory with all files> <path to project directory>

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-s", "--sourcedir"), type = "character",
              help = "Source directory with all files"),
  make_option(c("-d", "--destdir"), type = "character",
              help = "Destination directory. Should be /path/to/OMPARE/data/PNOC008-13/ for Patient 13")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
sourceDir <- opt$sourcedir
topDir <- opt$destdir

# specify destination directories
cnvdir <- file.path(topDir, 'CNV')
clinicaldir <- file.path(topDir, 'Clinical')
exprdir <- file.path(topDir, 'ExpressionGene')
fusionsdir <- file.path(topDir, 'Fusions')
immunescores <- file.path(topDir, 'ImmuneScores')
gsvascores  <-  file.path(topDir, 'GSVA')
mutdir <- file.path(topDir, 'MutationsMAF')
reports <- file.path(topDir, 'Reports')
summary <- file.path(topDir, 'Summary')

# create destination directories
system(paste0('mkdir -p ', cnvdir))
system(paste0('mkdir -p ', clinicaldir))
system(paste0('mkdir -p ', exprdir))
system(paste0('mkdir -p ', fusionsdir))
system(paste0('mkdir -p ', immunescores))
system(paste0('mkdir -p ', gsvascores))
system(paste0('mkdir -p ', mutdir))
system(paste0('mkdir -p ', reports))
system(paste0('mkdir -p ', summary))

# organize data
# copy number
cmd <- file.path(sourceDir, '*.{CNVs.p.value.txt,controlfreec.ratio.txt}')
cmd <- paste0('mv ', cmd, ' ', cnvdir)
system(cmd)
# clinical file
cmd <- file.path(sourceDir, 'patient_report.txt')
cmd <- paste0('mv ', cmd, ' ', clinicaldir)
system(cmd)
# expression
cmd <- file.path(sourceDir, '*rsem*')
cmd <- paste0('mv ', cmd, ' ', exprdir)
system(cmd)
# fusions
cmd <- file.path(sourceDir, '*fusion*')
cmd <- paste0('mv ', cmd, ' ', fusionsdir)
system(cmd)
# mutations - somatic/germline
cmd <- file.path(sourceDir, '*.{maf,hg38_multianno.txt.gz}')
cmd <- paste0('mv ', cmd, ' ', mutdir)
system(cmd)



