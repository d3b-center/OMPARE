# functions to get tpm matrix subsets for various modules
# source functions
source(file.path(patient_level_analyses_utils, 'quiet.R'))
source(file.path(patient_level_analyses_utils, 'batch_correct.R'))

# step 1
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
           short_histology = tumorType,
           broad_histology = tumorType,
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
  expr_corrected <- quiet(batch.correct(mat = expr_uncorrected, clin = clinical))
  
  # return batch corrected matrix and combined clinical
  res <- list(expr = expr_corrected, expr_uc = expr_uncorrected, clinical = clinical)
  return(res)
}

# subset dataset
# only patient of interest and tcga/pbta (use uncorrected TPM)
get_data_for_immune_profile <- function(expr_uncorrected, expr1, sampleInfo){
  # for immune profile
  # now get only sample of interest + samples from expr1
  smps <- c(colnames(expr1), sampleInfo$subjectID)
  expr_uncorrected_subset <- expr_uncorrected[,colnames(expr_uncorrected) %in% smps]
  return(expr_uncorrected_subset)
}

# get 10000 most variable genes for umap plotting & nearest neighbor analysis
get_most_variable_for_umap <-  function(expr_corrected){
  # filter out low expression genes
  maxVals <- apply(expr_corrected, FUN = max, MARGIN = 1)
  expr_filtered <- expr_corrected[maxVals > 20,]
  
  # for dimensionality reduction visualization (dim_reduction_plot.R)
  rv <- matrixStats::rowVars(as.matrix(expr_filtered))
  select <- order(rv, decreasing=TRUE)[seq_len(10000)]
  expr_most_var <- expr_filtered[select,]
  
  # return output
  return(expr_most_var)
}

# get umap output
get_umap_output <- function(expr_most_var){
  set.seed(100)
  umap_out <- uwot::umap(X = t(log2(expr_most_var+1)), n_neighbors = 21, n_components = 2, metric = "correlation", ret_nn = TRUE, n_sgd_threads = 123L)
  
  # add colnames/rownames to embeddings
  colnames(umap_out$embedding) <- c("UMAP1", "UMAP2")
  rownames(umap_out$embedding) <- colnames(expr_most_var)
  
  # add rownames to nearest neighbor
  rownames(umap_out$nn$correlation$idx) <- colnames(expr_most_var)
  rownames(umap_out$nn$correlation$dist) <- colnames(expr_most_var)
  
  # return output
  return(umap_out)
}

# extract nearest neighbor info
extract_umap_nearest_neighbor_table <- function(umap_out, expr_most_var, sampleInfo){
  corr <- as.data.frame(umap_out$nn$correlation$idx) # nn
  dist <- as.data.frame(umap_out$nn$correlation$dist) # distances
  corr <- t(apply(corr, MARGIN = 1, FUN = function(x) colnames(expr_most_var)[x]))
  nn_table <- data.frame(nearest_neighbor = as.character(corr[grep(sampleInfo$subjectID, rownames(corr)),]), 
                         distance = as.numeric(dist[grep(sampleInfo$subjectID, rownames(dist)),]))
  nn_table$distance <- round(nn_table$distance, digits = 3)
  
  # nearest neighbor table
  return(nn_table)
}

