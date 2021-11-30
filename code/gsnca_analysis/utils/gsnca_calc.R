
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("GSAR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("EGSEA"))
suppressPackageStartupMessages(library("DGCA"))


######################### prepare for GSNCA tests-filter out lowly expressed genes 
filter_low_expr_df <- function(expr_df){
  # filter lowly expressed genes by DGSA
  expr_df_filtered <- filterGenes(expr_df, 
                                  filterTypes = c("central", "dispersion"),
                                  filterDispersionType = "cv", 
                                  filterDispersionPercentile = 0.2,
                                  sequential= TRUE)
  
  # also additionally filter on sd to avoid error message on GSNCA step
  expr_df_filtered_sd <- apply(expr_df_filtered, 1, sd, na.rm = TRUE)
  expr_df_filtered_sd_filter <- which(expr_df_filtered_sd <= 0.015) %>% as.data.frame() %>% rownames()
  expr_df_filtered <- expr_df_filtered[!(row.names(expr_df_filtered) %in% expr_df_filtered_sd_filter),]
  
  return(expr_df_filtered)
}

####### get the function of plotMST2.pathway and modify for return of hub genes
return_hub_gene <- function(object, group, name=NULL, cor.method="pearson", min.sd=1e-3){
    nv <- ncol(object)
    object <- object[,c(which(group == 1), which(group == 2))]
    nv1 <- sum(group == 1)
    if(length(rownames(object)) < nrow(object)) 
      gnames <- as.character(c(1:nrow(object))) else gnames <- rownames(object)
    
    objt <- aperm(object, c(2,1))
    group1 <- objt[1:nv1,]
    group2 <- objt[(nv1+1):nv,]
    
    cormat1 <- abs(cor(group1, method=cor.method))
    cormat2 <- abs(cor(group2, method=cor.method))
    e1 <- eigen(cormat1)
    e2 <- eigen(cormat2)
    p1 <- matrix(abs(e1$vectors[,1]))
    p2 <- matrix(abs(e2$vectors[,1]))
    p1 <- p1 * norm(p1)
    p2 <- p2 * norm(p2)
    colnames(p1) <- "class1"
    colnames(p2) <- "class2"
    rownames(p1) <- rownames(p2) <- gnames
    major1.val <- max(p1)
    major2.val <- max(p2)
    major1.ind <- which.max(p1)
    major2.ind <- which.max(p2)
    MST2.group1 <- findMST2(object[,c(1:nv1)], cor.method, min.sd, TRUE)
    MST2.group2 <- findMST2(object[,c((nv1+1):nv)], cor.method, min.sd, TRUE)
    
    return(gnames[major1.ind])
  }


#### Get the pathway information and select genes in the pathways---------------
build_pathways <- function(gene_list) {
  # get entrezID for gene in expression pathway
  genes_df <- as.data.frame(gene_list) 
  colnames(genes_df) <- "symbol"
  
  genes_df <- genes_df %>%
    dplyr::mutate(eg = mapIds(org.Hs.eg.db, symbol, "ENTREZID", "SYMBOL")) 
  
  # get GSCollectionSet object
  pathway_build <- EGSEA::buildMSigDBIdx(entrezIDs = genes_df$eg, 
                                         geneSets="c2",
                                         species = "Homo sapiens")
  
  # get annotation to filter out disease related
  pathway_build_anno <- pathway_build[["c2"]]@anno %>% as.data.frame() %>% 
    dplyr::select(c("ID", "GeneSet", "Description")) 
  
  # get gene set IDs
  pathway_build_genesets <- pathway_build_anno %>% pull(GeneSet)
  
  # build df with pathways and entrezID of genes
  pathway_build_genesets_list <- lapply(pathway_build_genesets, function(x){
    pathway_df <- pathway_build[["c2"]]@idx[[x]] %>% as.data.frame() %>%
      mutate(GeneSet = x)
  })
  pathway_build_genesets_df <- do.call(rbind, pathway_build_genesets_list) %>%
    dplyr::left_join(pathway_build_anno) 
  colnames(pathway_build_genesets_df) <- c("eg", "description", "pathway", "type")
  pathway_build_genesets_df$eg <- as.character(pathway_build_genesets_df$eg)
  
  # annotate gene symbol and ensemble IDs 
  pathway_build_genesets_df <-pathway_build_genesets_df %>%
    dplyr::left_join(genes_df) %>%
    # filter the pathway file to contain only symbols available in expression
    dplyr::filter(!is.na(symbol))
  
  return(pathway_build_genesets_df)
}

#### function to create barplots for top pathways -----------------------------
pathway_barplots <- function(dat, title){
  # calculate log score
  dat <- dat %>% 
    mutate(log_score = (-1)*log10(pvalue)) %>%
    arrange(log_score, descending = TRUE)
  
  dat[,"pathway_description"] <- factor(dat[,"pathway_description"], levels = unique(dat[,"pathway_description"]))
  
  p <- ggplot(dat, aes(x = pathway_description, 
                       y = log_score,
                       fill = log_score)) + 
    geom_bar(stat="identity") + coord_flip() + theme_bw() +
    xlab("") + 
    ylab("-log10 P-Value") + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
    ggtitle(title) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
    guides(fill = "none")
  
  return(p)
}


######################### Run GSNCA and output plots and text files --------------
gsnca_analysis_plot <- function(similar_subjects_expr_df, ref_expr_df, ref_name, top_bar=20, top_net=5){
  ##### get annotation df for comparison groups
  # 20+1 POI group
  similar_subjects_anno <- as.data.frame(colnames(similar_subjects_expr_df)) 
  colnames(similar_subjects_anno) <- "subject_id"
  similar_subjects_anno <- similar_subjects_anno %>%
    dplyr::mutate(group="1")

  # reference group
  ref_subjects_anno <- as.data.frame(colnames(ref_expr_df)) 
  colnames(ref_subjects_anno) <- "subject_id"
  ref_subjects_anno<- ref_subjects_anno %>%
    dplyr::mutate(group="2")
  
  # combine both 
  combined_anno <- bind_rows(similar_subjects_anno, ref_subjects_anno)
  
  # generate combined matrix 
  similar_subjects_expr_df <- similar_subjects_expr_df %>% tibble::rownames_to_column("geneID")
  ref_expr_df <- ref_expr_df %>% tibble::rownames_to_column("geneID")
  
  # combine to contain only genes that pass filter for both expression matrix
  genes_in_common <- intersect(similar_subjects_expr_df$geneID, ref_expr_df$geneID) %>% unique()
  combined_matrix <- left_join(similar_subjects_expr_df[(similar_subjects_expr_df$geneID %in% genes_in_common), ], 
                               ref_expr_df[(ref_expr_df$geneID %in% genes_in_common),]) %>%
    tibble::column_to_rownames("geneID") 
  
  # generate pathway df 
  path_build_df <- build_pathways(rownames(combined_matrix))
  
  # only look at pathways that has 10-500 members
  pathway_list <- path_build_df %>% 
    group_by(pathway) %>% 
    mutate(n=n()) %>% 
    filter(n>=10 & n <=500) %>% 
    pull(pathway) %>% 
    unique()
  
  # define matrix to store results
  gsnca_results <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(gsnca_results) <- c("pathway_id", "pathway_description", "pvalue", "hub_gene")
  
  for(i in 1:length(pathway_list)){
    # iterate through path list 
    pathway_of_interest <- pathway_list[i]
    
    # get the description 
    description <- path_build_df %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(description) %>% unique()
    
    # find genes in pathway of interest
    genes_in_pathway <- path_build_df %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(symbol) %>% unique()
    
    # filter to genes in target pathway
    combined_matrix_per_pathway <- combined_matrix[row.names(combined_matrix) %in% genes_in_pathway,]
    
    # run GSNCA test on filtered 
    result_pval<-GSNCAtest(object=as.matrix(combined_matrix_per_pathway), 
                           # since the matrix is selected by order of row, the group will match
                           group=combined_anno$group, 
                           nperm=1000, 
                           cor.method="spearman", 
                           check.sd=TRUE, 
                           min.sd=1e-3, 
                           max.skip=10
    )
    
    # store the pathway name and results in the results table
    gsnca_results[i,1] <- pathway_of_interest
    gsnca_results[i,2] <- description
    gsnca_results[i,3] <- result_pval
    
    # output the hubgene for each pathway analysis
    gsnca_results[i,4] <- return_hub_gene(object=as.matrix(combined_matrix_per_pathway),
                                           # since the matrix is selected by order of row, the group will match
                                           group=combined_anno$group,
                                           cor.method="spearman")
  }
  
  # write out results
  gsnca_results <- gsnca_results %>% 
    dplyr::filter(!is.na(pathway_id)) %>%
    arrange(pvalue, descending = FALSE) %>%
    # filter to pval < 0.05
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::mutate(comparison = ref_name)
  
  gsnca_results %>% 
    readr::write_tsv(file.path(output_dir, paste0("top20_similar_vs_", ref_name, "_GSNCA_analysis.tsv")))

  # we plot out top n for networks
  gsnca_top_net <- gsnca_results %>%
    slice_head(n = top_net) 
  
  pdf(file = file.path(output_dir, paste0("top20_similar_vs_", ref_name, "_GSNCA_plots.pdf")))
  
  #### For pathways with top 5 pval, we plot out the network plots for them
  for(j in 1:nrow(gsnca_top_net)){
    # gather genes in the pathway of interest 
    pathway_of_interest <- gsnca_top_net[j,1]
    description <- gsnca_top_net[j,2]
    # find genes in pathway of interest
    genes_in_pathway <- path_build_df %>% 
      filter(pathway == pathway_of_interest) %>% 
      pull(symbol) %>% unique()
    
    # filter to genes in target pathway
    combined_matrix_per_pathway <- combined_matrix[row.names(combined_matrix) %in% genes_in_pathway,]
    
    plotMST2.pathway(object=as.matrix(combined_matrix_per_pathway),
                     # since the matrix is selected by order of row, the group will match
                     group=combined_anno$group,
                     cor.method="spearman",
                     group1.name="Patients of Interest",
                     group2.name=ref_name,
                     legend.size=0.9,
                     label.size=1.2,
                     name=paste0(pathway_of_interest, " ", description)) 
    
  }
  
  #### For pathways with top n pval, we plot out as bar plot
  gsnca_top_bar <- gsnca_results %>%
    slice_head(n = top_bar) 
  barplot <- pathway_barplots(gsnca_top_bar, title = paste0("GSNCA: TOP ", top_bar, " Similar Subjects vs. ", ref_name))
  print(barplot)
  
  dev.off()

}



