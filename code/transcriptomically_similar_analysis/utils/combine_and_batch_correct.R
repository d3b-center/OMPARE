# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")

source(file.path(root_dir, "code", "utils", 'quiet.R'))
source(file.path(root_dir, "code", "utils", 'batch_correct.R'))

combine_and_batch_correct <- function(ref_expr, ref_clinical, subject_expr, subject_clinical){
  
  common_cols <- intersect(colnames(subject_clinical), colnames(ref_clinical))
  
  # combine both clinical files
  clinical <- ref_clinical %>%
    dplyr::select(common_cols) %>%
    rbind(subject_clinical %>%
            dplyr::select(common_cols)) %>%
    column_to_rownames("Kids_First_Biospecimen_ID")
  
  # combine tpm matrix
  expr_uncorrected <- ref_expr %>%
    rownames_to_column('gene_symbol') %>%
    inner_join(subject_expr %>%
                 rownames_to_column("gene_symbol"), by = 'gene_symbol') %>%
    column_to_rownames('gene_symbol')
  
  # match clinical file to expression data
  clinical <- clinical[colnames(expr_uncorrected),]
  
  # correct for batch effect: study_id + library_name
  clinical <- clinical %>%
    mutate(batch = paste0(cohort, "_", RNA_library))
  expr_corrected <- quiet(batch_correct(mat = expr_uncorrected, clin = clinical))
  
  # return batch corrected matrix and combined clinical
  res <- list(expr_corrected = expr_corrected, expr_uncorrected = expr_uncorrected, clinical = clinical)
  return(res)
}
