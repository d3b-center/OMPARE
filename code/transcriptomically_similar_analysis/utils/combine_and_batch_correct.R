# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "transcriptomically_similar_analysis")

source(file.path(root_dir, "code", "utils", 'quiet.R'))
source(file.path(root_dir, "code", "utils", 'batch_correct.R'))

combine_and_batch_correct <- function(expr1, expr2 = pnoc008_tpm, clinical1, clinical2 = pnoc008_clinical){
  # common columns
  cols <- c('kids_first_biospecimen_id', 'sample_id', 'cohort_participant_id', 'subject_id',
            'gender', 'ethnicity', 'age_at_diagnosis_in_days',
            'short_histology', 'broad_histology', 
            'pathology_diagnosis', 'integrated_diagnosis',
            'study_id', 'library_name')
  
  # subset clinical file
  if(unique(clinical1$study_id) == "PBTA"){
    clinical1 <- clinical1 %>%
      mutate(kids_first_biospecimen_id = Kids_First_Biospecimen_ID,
             subject_id = Kids_First_Biospecimen_ID) %>%
      dplyr::select(cols)
  } else if(unique(clinical1$study_id == "TCGA")){
    clinical1 <- clinical1 %>%
      mutate(kids_first_biospecimen_id = '',
             subject_id = sample_barcode,
             cohort_participant_id = '') %>%
      dplyr::select(cols)
  }
  
  clinical2 <- clinical2 %>%
    mutate(kids_first_biospecimen_id = Kids_First_Biospecimen_ID,
           subject_id = subjectID,
           gender = sex,
           age_at_diagnosis_in_days = age_diagnosis_days,
           short_histology = 'HGAT',
           broad_histology = 'Diffuse astrocytic and oligodendroglial tumor',
           pathology_diagnosis = tumorType, 
           integrated_diagnosis = tumorType) %>%
    dplyr::select(cols)
  
  # combine both clinical files
  clinical <- clinical1 %>%
    rbind(clinical2) %>%
    mutate(tmp = subject_id,
           batch = paste0(study_id, '_', library_name)) %>%
    remove_rownames() %>%
    column_to_rownames('tmp')
  
  # combine tpm matrix
  expr_uncorrected <- expr1 %>%
    rownames_to_column('genes') %>%
    inner_join(expr2 %>%
                 rownames_to_column('genes'), by = 'genes') %>%
    column_to_rownames('genes')
  
  # match clinical file to expression data
  clinical <- clinical[colnames(expr_uncorrected),]
  
  # correct for batch effect: study_id + library_name
  expr_corrected <- quiet(batch_correct(mat = expr_uncorrected, clin = clinical))
  
  # return batch corrected matrix and combined clinical
  res <- list(expr = expr_corrected, expr_uc = expr_uncorrected, clinical = clinical)
  return(res)
}
