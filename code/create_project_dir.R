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
clinicaldir <- file.path(topDir, 'Clinical')
cnvdir <- file.path(topDir, 'CNV')
exprdir <- file.path(topDir, 'ExpressionGene')
fusionsdir <- file.path(topDir, 'Fusions')
mutdir <- file.path(topDir, 'MutationsMAF')
immunescores <- file.path(topDir, 'ImmuneScores')
gsvascores  <-  file.path(topDir, 'GSVA')
reports <- file.path(topDir, 'Reports')
summary <- file.path(topDir, 'Summary')

# create directories to move input files
dir.create(clinicaldir, showWarnings = F, recursive = T)
dir.create(cnvdir, showWarnings = F, recursive = T)
dir.create(exprdir, showWarnings = F, recursive = T)
dir.create(fusionsdir, showWarnings = F, recursive = T)
dir.create(mutdir, showWarnings = F, recursive = T)

# create directories to save output
dir.create(immunescores, showWarnings = F, recursive = T)
dir.create(gsvascores, showWarnings = F, recursive = T)
dir.create(reports, showWarnings = F, recursive = T)
dir.create(summary, showWarnings = F, recursive = T)

# organize data
# copy number
cmd <- file.path(sourceDir, '*.{CNVs.p.value.txt,controlfreec.ratio.txt}')
cmd <- paste('mv', cmd, cnvdir)
system(cmd)

# clinical file
cmd <- file.path(sourceDir, 'patient_report.txt')
cmd <- paste('mv', cmd, clinicaldir)
system(cmd)

# expression
cmd <- file.path(sourceDir, '*rsem*')
cmd <- paste('mv', cmd, exprdir)
system(cmd)

# fusions
cmd <- file.path(sourceDir, '*fusion*')
cmd <- paste('mv', cmd, fusionsdir)
system(cmd)

# mutations - somatic/germline
cmd <- file.path(sourceDir, '*.{maf,hg38_multianno.txt.gz}')
cmd <- paste('mv', cmd, mutdir)
system(cmd)



