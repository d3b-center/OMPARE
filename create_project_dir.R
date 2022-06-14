# Author: Komal S. Rathi
# Function: code to create and organize directory structure

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--source_dir"), type = "character",
              help = "source directory with files downloaded from data delivery project"),
  make_option(c("--patient_dir"), type = "character",
              help = "patient-specific output directory")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
source_dir <- opt$source_dir
patient_dir <- opt$patient_dir

# specify destination directories
cnvdir <- file.path(patient_dir, 'copy-number-variations')
mutdir <- file.path(patient_dir, 'simple-variants')
reports <- file.path(patient_dir, 'reports')
output <- file.path(patient_dir, 'output')

# create directories to move input files
dir.create(cnvdir, showWarnings = F, recursive = T)
dir.create(mutdir, showWarnings = F, recursive = T)

# create directories to save output
dir.create(output, showWarnings = F, recursive = T)
dir.create(reports, showWarnings = F, recursive = T)

# organize patient specific files
# copy number
# cnvkit diagram.pdf file
cmd <- file.path(source_dir, '*.diagram.pdf')
cmd <- paste('mv', cmd, cnvdir)
system(cmd)

# mutations
# somatic and germline variant files
cmd <- file.path(source_dir, '*.{*_somatic.norm.annot.protected.maf,consensus_somatic.protected.maf,gatk.PASS.vcf.gz.hg38_multianno.txt.gz}')
cmd <- paste('mv', cmd, mutdir)
system(cmd)