################################################################################
# This script conducts principal component analysis (PCA) and generates a diagnostic plot for optimization of tSNE, using prcomp function.
#
#
# Written by Adam Kraya for PNOC008, 2020
#
#
#
# ####### USAGE, assumed to be run from top-level of project:
# Rscript --vanilla 'code/runPCAdiagnostic.R --input <expression input file> --output <output files for writing image pdf files>
#     --input_file: The name of the input rds expression data file to use for calculating PCA.
#     --output_file: The name of the pdf output file for PCA proportion of variance explained (pve)

knitr::opts_chunk$set(echo = TRUE, error = TRUE, results = "hide", fig.width=8, fig.height=8)
if(!require('pacman')) {
  install.packages('pacman')
}

if (!("BiocManager" %in% installed.packages())) {
  install.packages("BiocManager")
}
library(BiocManager, quietly = TRUE)
if (!("pcaMethods" %in% installed.packages())) {
  BiocManager::install("pcaMethods")
}
library(pcaMethods)
pacman::p_load(glmnet, tidyverse, data.table, dplyr, gmodels, IDPmisc)

#### Set Up optparse --------------------------------------------------------------------

## Define arguments

option_list <- list(
  optparse::make_option(
    c("--input_file"),
    type = "character",
    default = NA,
    help = "The input file of expression data from which scores will be calculated."
  ),
  optparse::make_option(
    c("--output_file"),
    type = "character",
    default = NA,
    help = "The output file for writing a PCA diagnostic plot in /data/pca/OpenPBTA"
  )
)


## Read in arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
if (is.na(opt$input_file)) stop("\n\nERROR: You must provide an input file with TPM expression data with the flag --input")
if (is.na(opt$output_file)) stop("\n\nERROR: You must provide an output file for saving a PCA diagnostic plot with the flag --output, assumed to be placed in the `results/` directory of this analysis.")


#### Set Up paths and file names --------------------------------------------------------------------

# If the output directory does not exist, create it
output_dir <- dirname(opt$output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Ensure the input file exists in `data/` and specify input/output files
expression_data_file <- opt$input_file
if (!file.exists(expression_data_file)) stop("\n\nERROR: Provided input file does not exist.")
image_output_file <- file.path(output_dir, basename(opt$output_file))

#### Load input file --------------------------------------------------------------------
message('Reading in and re-scaling input expression data')
expression_data <- as.data.frame( readr::read_rds(expression_data_file) )

## Reformat and scale expression data
rownames(expression_data) = expression_data[,'gene_id']
expression_data=select(expression_data, select=-c(gene_id))
expression_data.t=t(expression_data)
expression_data.t.scale=prep(expression_data.t, scale='uv', center=T)

## Perform PCA and compute proportion of variance explained (pve) for plotting
message('Computing principle components and generating pve plot')
pc <- pca(expression_data.t.scale, method="svd", center=FALSE, nPcs=100)
pve <- 100* (sDev(pc))/sum ((sDev(pc))^2)

pdf(image_output_file)

  plot(pve, pch=16, 
       main = 'PCA Analysis of OpenPBTA TPM Stranded Data',
       xlab="Principal Components",
       ylab="Prop. of variance explained"
       )
  
dev.off()



