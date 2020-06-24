# ssGSEA of top 20 genomically similar patients
# takes TPM data
library(msigdbr) ## Contains the hallmark data sets
library(GSVA)    ## Performs GSEA analysis

ssGSEA <- function(topCor, fname) {
  
  if(!file.exists(fname)) {
    expression_data <- topCor
    human_hallmark  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") ## human hallmark genes from `migsdbr` package. The loaded data is a tibble.
    
    # Prepare expression data: log2 transform re-cast as matrix
    ### Rownames are genes and column names are samples
    expression_data_log2_matrix <- as.matrix( log2(expression_data + 1) )
    
    # Prepare hallmark genes: Create a list of hallmarks, each of which is a list of genes
    human_hallmark_twocols <- human_hallmark %>% dplyr::select(gs_name, human_gene_symbol)
    human_hallmark_list    <- base::split(human_hallmark_twocols$human_gene_symbol, list(human_hallmark_twocols$gs_name))
    
    # We then calculate the Gaussian-distributed scores
    gsea_scores <- GSVA::gsva(expression_data_log2_matrix,
                              human_hallmark_list,
                              method = "gsva",
                              min.sz = 1, max.sz = 1500,
                              mx.diff = TRUE, ## Setting this argument to TRUE computes Gaussian-distributed scores (bimodal score distribution if FALSE)
                              verbose = FALSE)        
    
    ### Clean scoring into tidy format
    gsea_scores_df <- as.data.frame(gsea_scores) %>%
      rownames_to_column(var = "hallmark_name")
    
    #first/last_bs needed for use in gather (we are not on tidyr1.0)
    first_bs <- head(colnames(gsea_scores), n=1)
    last_bs  <- tail(colnames(gsea_scores), n=1)
    
    gsea_scores_df_tidy <- gsea_scores_df %>%
      tidyr::gather(Kids_First_Biospecimen_ID, gsea_score, !!first_bs : !!last_bs) %>%
      dplyr::select(Kids_First_Biospecimen_ID, hallmark_name, gsea_score)
    
    # sort by median
    gsea_scores_df_tidy <- gsea_scores_df_tidy %>%
      group_by(hallmark_name) %>%
      arrange(hallmark_name, gsea_score) %>%
      mutate(gsea_score_median = median(gsea_score))
    
    # factorize by median
    tmp <- gsea_scores_df_tidy %>%
      dplyr::select(hallmark_name, gsea_score_median) %>%
      arrange(desc(gsea_score_median)) %>%
      unique() %>%
      .$hallmark_name
    gsea_scores_df_tidy$hallmark_name <- factor(gsea_scores_df_tidy$hallmark_name, levels = tmp)
    gsea_scores_df_tidy[,"IsSample"] <- ifelse(grepl(sampleInfo$subjectID, gsea_scores_df_tidy$Kids_First_Biospecimen_ID), T, F)
    write.table(gsea_scores_df_tidy, file = fname, sep = "\t", row.names = F)
  } else {
    gsea_scores_df_tidy <- read.delim(fname, check.names = F)
  }
  # plot as boxplot
  p <- ggplot(gsea_scores_df_tidy, aes(hallmark_name, gsea_score)) + 
    geom_boxplot(outlier.shape = NA) +  
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 75, hjust = 1))
  raw.scoresSample <- gsea_scores_df_tidy[gsea_scores_df_tidy$IsSample == T,]
  p <- p + 
    geom_point(data = raw.scoresSample, aes(hallmark_name, gsea_score), colour = "red", size = 3, shape = "triangle") +
    theme(axis.text = element_text(size = 8, face = "bold"), 
          axis.title = element_blank())
  return(p)  
}
