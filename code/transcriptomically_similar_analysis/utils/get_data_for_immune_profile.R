# only patient of interest and tcga/pbta (use uncorrected TPM)
get_data_for_immune_profile <- function(expr_uncorrected, ref_expr, patient_of_interest){
  # for immune profile
  # now get only sample of interest + samples from ref_expr
  smps <- c(colnames(ref_expr), patient_of_interest)
  expr_uncorrected_subset <- expr_uncorrected[,colnames(expr_uncorrected) %in% smps]
  return(expr_uncorrected_subset)
}
