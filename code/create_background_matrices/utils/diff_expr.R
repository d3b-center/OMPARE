# script to perform differential expression between pediatric and normal samples
# this will be used by Oncogrid 
# TMM normalization takes a very long time especially because the N for normals is large
# so we will use z-score method (just like PNOC003) and revisit this later

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir('.git'))
data_dir <- file.path(root_dir, "data")

# function
source(file.path(root_dir, "code", "rnaseq_analysis", "utils", "ss_diffexpr.R"))

# z-score and return only patient's value (i.e. last column)
get_zscore <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out[length(out)])
}

calc_degs <- function(expData, gtexData, thresh = 2.0) {
  
  sample_of_interest <- unique(expData$sample)
  print(sample_of_interest)
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(gtexData), expData$gene_symbol)
  mergeDF <- gtexData %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(expData %>%
                 dplyr::select(-c(sample)), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  colnames(mergeDF)[ncol(mergeDF)] <- sample_of_interest 
  
  # Filter in Patient: TPM > 10
  mergeDF <- mergeDF[mergeDF[,sample_of_interest] > 10,] 
  
  # z-score
  output <- apply(mergeDF, FUN = get_zscore, MARGIN = 1) 
  
  # full data
  # combine tpm and z-score for sample of interest
  genes_df <- data.frame(z_score = output, tpm = mergeDF[names(output),sample_of_interest], sample = sample_of_interest)
  genes_df <- genes_df %>%
    mutate(diff_expr = ifelse(z_score < (-1*thresh), "down", 
                              ifelse(z_score > thresh, "up", NA))) %>%
    rownames_to_column("gene_symbol")
  return(genes_df)
}

# ssexpr
calc_degs_ssexpr <- function(expData_counts, gtexData_counts) {
  
  sample_of_interest <- unique(expData_counts$sample)
  print(sample_of_interest)
  
  # Merge GTEx and Patient data on common genes
  intGenesTmp <- intersect(rownames(gtexData_counts), expData_counts$gene_symbol)
  mergeDF_counts <- gtexData_counts %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(expData_counts %>%
                 dplyr::select(-c(sample)), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol")
  colnames(mergeDF_counts)[ncol(mergeDF_counts)] <- "sample_of_interest" 
  
  # apply single sample differential expression
  genes_df <- ss_diffexpr(expr = mergeDF_counts, norm_method = "tmm", housekeeping_genes = NULL)
  genes_df$sample <- sample_of_interest
  genes_df <- genes_df %>%
    rownames_to_column("gene_symbol")
  return(genes_df)
}

diff_expr <- function(hist_file, normal_tissue_dir, pediatric_cancer_dir, output_dir, prefix){
  # normal tissue tpm
  normal_tissue_files <- list.files(path = normal_tissue_dir, full.names = T)
  normal_tissue_tpm <- normal_tissue_files[grep('tpm', normal_tissue_files)]
  normal_tissue_tpm <- readRDS(normal_tissue_tpm)
  
  # normal tissue counts
  # normal_tissue_counts <- normal_tissue_files[grep('count', normal_tissue_files)]
  # normal_tissue_counts <- readRDS(normal_tissue_counts)
  
  # pediatric tumor tpm
  pediatric_cancer_files <- list.files(path = pediatric_cancer_dir, full.names = T)
  pediatric_cancer_tpm <- pediatric_cancer_files[grep('tpm', pediatric_cancer_files)]
  pediatric_cancer_tpm <- readRDS(pediatric_cancer_tpm)
  pediatric_cancer_tpm <- pediatric_cancer_tpm %>% 
    rownames_to_column("gene_symbol") %>%  
    gather('sample', "tpm", -gene_symbol) 
  
  # pediatric tumor counts
  # pediatric_cancer_files <- list.files(path = pediatric_cancer_dir, full.names = T)
  # pediatric_cancer_counts <- pediatric_cancer_files[grep('count', pediatric_cancer_files)]
  # pediatric_cancer_counts <- readRDS(pediatric_cancer_counts)
  # pediatric_cancer_counts <- pediatric_cancer_counts %>% 
  #   rownames_to_column("gene_symbol") %>%  
  #   gather('sample', "tpm", -gene_symbol) 

  # do single sample deg
  
  # method 1: z-score
  # to be conservative, we have increased the threshold to 2 from 1.5
  # per sample comparison
  # user  system elapsed 
  # 1.362   1.505   3.798
  res <- plyr::ddply(pediatric_cancer_tpm, 
                     .variables = "sample", 
                     .fun = function(x) calc_degs(expData = x, gtexData = normal_tissue_tpm, thresh = 2))
  
  # method 2: TMM normalization
  # user  system elapsed 
  # 59.091  15.276  77.606
  # res <- plyr::ddply(pediatric_cancer_counts, 
  #                    .variables = "sample", 
  #                    .fun = function(x) calc_degs_ssexpr(expData_counts = x, gtexData_counts = normal_tissue_counts))
  
  # so for the oncogrid, we will go with the z-score because it will take a long time to process hundreds of samples
  res <- res %>%
    filter(!is.na(diff_expr))
  
  # combine with histology
  res <- res %>%
    inner_join(hist_file %>%
    dplyr::select(Kids_First_Biospecimen_ID, sample_id, cohort, cohort_participant_id), by = c("sample" = "Kids_First_Biospecimen_ID"))

  saveRDS(res, file = file.path(output_dir, paste(prefix, 'degs.rds', sep = "_")))
}
