# function to create histology file and rnaseq matrices
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
})

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# protein coding genes and remove HISTONE genes
# read gencode from OpenPedCan
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, 'OpenPedCan-analysis/data/gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_pc <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding",
         !grepl("^HIST", gene_name)) %>%
  unique()

create_rnaseq_matrices <- function(hist_file, tpm_mat, counts_mat, collapse = FALSE, output_dir, prefix){
  
  # create output directory
  dir.create(output_dir, showWarnings = F, recursive = T)
  
  # subset to rnaseq and keep only those samples overlapping with expression matrix
  hist_file <- hist_file %>%
    filter(experimental_strategy == "RNA-Seq",
           Kids_First_Biospecimen_ID %in% colnames(tpm_mat))
  
  # tpm
  if(collapse == TRUE){
    tpm_mat <- tpm_mat %>%
      dplyr::select(c('gene_id',hist_file$Kids_First_Biospecimen_ID))
    tpm_mat <- collapse_rnaseq(expr.mat = tpm_mat)
  } else {
    tpm_mat <- tpm_mat %>%
      dplyr::select(hist_file$Kids_First_Biospecimen_ID)
  }
  tpm_mat <- tpm_mat %>%
    rownames_to_column("gene") %>%
    filter(gene %in% gencode_pc$gene_name) %>%
    column_to_rownames("gene")
  tpm_output <- file.path(output_dir, paste(prefix, "tpm.rds", sep = "_"))
  saveRDS(tpm_mat, file = tpm_output)
  
  # counts
  if(collapse == TRUE){
    counts_mat <- counts_mat %>%
      dplyr::select(c('gene_id',hist_file$Kids_First_Biospecimen_ID))
    counts_mat <- collapse_rnaseq(expr.mat = counts_mat)
  } else {
    counts_mat <- counts_mat %>%
      dplyr::select(hist_file$Kids_First_Biospecimen_ID)
  }
  counts_mat <- counts_mat %>%
    rownames_to_column("gene") %>%
    filter(gene %in% gencode_pc$gene_name) %>%
    column_to_rownames("gene")
  counts_output <- file.path(output_dir, paste(prefix, "counts.rds", sep = "_"))
  saveRDS(counts_mat, file = counts_output)
}
