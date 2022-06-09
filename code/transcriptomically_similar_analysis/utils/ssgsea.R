# ssGSEA of top 20 genomically similar patients using TPM
suppressPackageStartupMessages({
  library(msigdbr) ## Contains the msigDB gene sets
  library(GSVA) ## Performs GSEA analysis
  library(dplyr)
})    

ssgsea <- function(nn_tpm_input, patient_of_interest) {
  expression_data <- nn_tpm_input
  human_geneset <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") ## human REACTOME genes from `migsdbr` package. The loaded data is a tibble.
  
  # Prepare expression data: log2 transform re-cast as matrix
  ### Rownames are genes and column names are samples
  expression_data_log2_matrix <- as.matrix(log2(expression_data + 1))
  
  # Prepare REACTOME genes: Create a list of REACTOME gene sets, each of which is a list of genes
  human_geneset_twocols <- human_geneset %>% dplyr::select(gs_name, human_gene_symbol)
  human_geneset_list    <- base::split(human_geneset_twocols$human_gene_symbol, list(human_geneset_twocols$gs_name))
  
  # We then calculate the Gaussian-distributed scores
  gsea_scores <- GSVA::gsva(expression_data_log2_matrix,
                            human_geneset_list,
                            method = "ssgsea",
                            min.sz = 1, max.sz = 1500,
                            mx.diff = TRUE, ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
                            verbose = FALSE)        
  
  ### Clean scoring into tidy format
  gsea_scores_df <- as.data.frame(gsea_scores) %>%
    rownames_to_column(var = "geneset_name")
  
  #first/last_bs needed for use in gather (we are not on tidyr1.0)
  first_bs <- head(colnames(gsea_scores), n=1)
  last_bs  <- tail(colnames(gsea_scores), n=1)
  
  gsea_scores_df_tidy <- gsea_scores_df %>%
    tidyr::gather(Kids_First_Biospecimen_ID, gsea_score, !!first_bs : !!last_bs) %>%
    dplyr::select(Kids_First_Biospecimen_ID, geneset_name, gsea_score)
  
  # add sample id
  gsea_scores_df_tidy[,"IsSample"] <- ifelse(grepl(patient_of_interest, gsea_scores_df_tidy$Kids_First_Biospecimen_ID), T, F)
  
  # calculate median score
  gsea_scores_df_tidy <- gsea_scores_df_tidy %>%
    group_by(geneset_name) %>%
    arrange(geneset_name, gsea_score) %>%
    dplyr::mutate(gsea_score_median = median(gsea_score))
  
  # top 50 pathways
  top50 <- gsea_scores_df_tidy %>%
    filter(IsSample) %>%
    dplyr::summarise(abs.score = abs(gsea_score - gsea_score_median)) %>%
    arrange(desc(abs.score)) %>%
    slice_head(n = 50)
  gsea_scores_df_tidy <- gsea_scores_df_tidy %>%
    filter(geneset_name %in% top50$geneset_name)
  
  # save output
  return(gsea_scores_df_tidy)
} 

