# capture upregulated pathways for all genomically similar patients

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# cancer Genes
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))

# genesets 
gene_set <- getGmt(file.path(data_dir, 'msigdb', 'c2.cp.reactome.v6.0.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)

# function to get upregulated pathways from top genomically similar patients  
pathway_analysis <- function(nn_table, normal_tissue, adult_cancer, pediatric_cancer,
                             filtered_cnv,
                             output_dir, prefix, comparison, patient_of_interest) {
  
  # source function for RNA-seq diffexpr & pathway analysis
  source(file.path(root_dir, "code", "rnaseq_analysis", "utils", "rnaseq_analysis_edgeR.R"))
  
  # Dataset1: GTex normals
  # counts
  normal_tissue_dir <- file.path(root_dir, "data", "normal_data")
  normal_counts <- list.files(path = normal_tissue_dir, pattern = "*counts.rds", full.names = T)
  normal_counts <- readRDS(normal_counts)
  
  # Dataset2: Adult cancer
  # counts
  adult_cancer_dir <- file.path(root_dir, "data", "adult_data")
  adult_counts <- list.files(path = adult_cancer_dir, pattern = "*counts.rds", full.names = T)
  adult_counts <- readRDS(adult_counts)
  
  # Dataset3: Pediatric cancer
  # counts
  pediatric_cancer_dir <- file.path(root_dir, "data", "pediatric_data")
  pediatric_counts <- list.files(path = pediatric_cancer_dir, pattern = "*counts.rds", full.names = T)
  pediatric_counts <- readRDS(pediatric_counts)
  
  # now convert counts to long format to compute single sample pathway analysis
  if(comparison == "pediatric"){
    expr_counts <- pediatric_counts %>% 
      rownames_to_column("gene_symbol") %>%
      gather(sample, expected_count, -c('gene_symbol')) %>%
      filter(sample %in% nn_table$nearest_neighbor)
    
    # tpm
    expr_tpm <- list.files(path = pediatric_cancer_dir, pattern = "*tpm.rds", full.names = T)
  } else if(comparison == "adult") {
    expr_counts <- adult_counts %>% 
      rownames_to_column("gene_symbol") %>%
      gather(sample, expected_count, -c('gene_symbol')) %>%
      filter(sample %in% nn_table$nearest_neighbor)
    
    # tpm
    expr_tpm <- list.files(path = adult_cancer_dir, pattern = "*tpm.rds", full.names = T)
  }
  
  # convert tpm expression to long format to compute single sample pathway analysis
  expr_tpm <- readRDS(expr_tpm)
  expr_tpm <- expr_tpm %>% 
    rownames_to_column("gene_symbol") %>%
    gather(sample, tpm, -c('gene_symbol')) %>%
    filter(sample %in% nn_table$nearest_neighbor)
  
  # function to read existing gsea output, process and add current subject to it if not already present
  gsea_enrichment <- function(gsea_output, exp.data.counts, exp.data.tpm, refData.counts, gene_set, comparison, cancer_genes) {
    
    # sample of interest
    patient_of_interest <- unique(as.character(exp.data.counts$sample))
    
    if(file.exists(gsea_output)){
      existing_gsea_output <- readRDS(gsea_output)
    } else {
      existing_gsea_output <- list()
    }
    
    # run only if sample does not exist
    if(length(intersect(names(existing_gsea_output), patient_of_interest)) == 0){
      res <- run_rnaseq_analysis_edger(exp.data.counts = exp.data.counts, 
                                       exp.data.tpm = exp.data.tpm, 
                                       refData.counts = refData.counts, 
                                       gene_set = gene_set, 
                                       comparison = comparison,
                                       cancer_genes = cancer_genes,
                                       single_sample = FALSE,
                                       sample_name = patient_of_interest)
      existing_gsea_output[[patient_of_interest]] <- res
      saveRDS(existing_gsea_output, file = gsea_output)
    } 
  }
  
  # compute pathway enrichment for 20 nearest neighbors vs normal tissue
  gsea_vs_normal <- file.path(output_dir, paste0(prefix, '_nn_vs_normals_gsea.rds'))
  if(!file.exists(gsea_vs_normal)){
    plyr::d_ply(.data = expr_counts, .variables = "sample", .fun = function(x)
      gsea_enrichment(gsea_output = gsea_vs_normal, 
                      exp.data.counts = x,
                      exp.data.tpm = expr_tpm,
                      refData.counts = normal_counts,
                      gene_set = gene_set,
                      comparison = paste0("Normal: ", normal_tissue),
                      cancer_genes = cancer_genes))
  }
  gsea_vs_normal <- readRDS(gsea_vs_normal)
  gsea_vs_normal <- plyr::ldply(gsea_vs_normal, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  
  # compute pathway enrichment for 20 nearest neighbors vs pediatric tumors
  gsea_vs_pediatric <- file.path(output_dir, paste0(prefix, '_nn_vs_pediatric_gsea.rds'))
  if(!file.exists(gsea_vs_pediatric)){
    plyr::d_ply(.data = expr_counts, .variables = "sample", .fun = function(x)
      gsea_enrichment(gsea_output = gsea_vs_pediatric, 
                      exp.data.counts = x,
                      exp.data.tpm = expr_tpm,
                      refData.counts = pediatric_counts,
                      gene_set = gene_set,
                      comparison = paste0("Pediatric: ", pediatric_cancer),
                      cancer_genes = cancer_genes))
  }
  gsea_vs_pediatric <- readRDS(gsea_vs_pediatric)
  gsea_vs_pediatric <- plyr::ldply(gsea_vs_pediatric, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  
  # compute pathway enrichment for 20 nearest neighbors vs adult tumors
  gsea_vs_adult <- file.path(output_dir, paste0(prefix, '_nn_vs_adult_gsea.rds'))
  if(!file.exists(gsea_vs_adult)){
    plyr::d_ply(.data = expr_counts, .variables = "sample", .fun = function(x)
      gsea_enrichment(gsea_output = gsea_vs_adult, 
                      exp.data.counts = x,
                      exp.data.tpm = expr_tpm,
                      refData.counts = adult_counts,
                      gene_set = gene_set,
                      comparison = paste0("Adult: ", adult_cancer),
                      cancer_genes = cancer_genes))
  }
  gsea_vs_adult <- readRDS(gsea_vs_adult)
  gsea_vs_adult <- plyr::ldply(gsea_vs_adult, .fun = function(x) return(x[[1]]), .id = 'sample_name')
  
  # combine all nearest neighbor gsea comparisons
  nn_comparisons <- rbind(gsea_vs_normal, gsea_vs_pediatric, gsea_vs_adult)
  
  # patient gsea comparisons
  patient_vs_normal <- readRDS(file.path(patient_dir, "output", "rnaseq_analysis", "patient_vs_normals_gsea.rds"))
  patient_vs_normal <- patient_vs_normal$pathways
  patient_vs_ped <- readRDS(file.path(patient_dir, "output", "rnaseq_analysis", "patient_vs_pediatric_gsea.rds"))
  patient_vs_ped <- patient_vs_ped$pathways
  patient_vs_adult <- readRDS(file.path(patient_dir, "output", "rnaseq_analysis", "patient_vs_adult_gsea.rds"))
  patient_vs_adult <- patient_vs_adult$pathways
  patient_comparisons <- rbind(patient_vs_normal, patient_vs_ped, patient_vs_adult)
  patient_comparisons$sample_name  <- patient_of_interest
  
  # now combine patient comparisons and other comparisons
  shared_pathways <- rbind(patient_comparisons, nn_comparisons)
  
  # different cutoff for pediatric and adult tumors
  n_perc = 12 # 60%
  
  # highly significant up/down pathways only (adj. pvalue < 0.01)
  # pathways which are seen misregulated in n.perc genomically similar samples
  shared_pathways <- shared_pathways %>%
    filter(padj < 0.01) %>% 
    dplyr::select(-c(ES, NES, size)) %>%
    mutate(pval = scientific(pval, digits = 3),
           padj = scientific(padj, digits = 3)) %>%
    group_by(pathway, comparison, direction) %>%
    mutate(sample_count_per_pathway = n()) %>%
    filter(sample_count_per_pathway >= n_perc) %>%
    arrange(desc(sample_name), desc(sample_count_per_pathway)) %>%
    as.data.frame()
  
  # now create table2 in which we will have genes, pathway and copy number info
  # this is only for patient of interest
  if(!is.null(filtered_cnv)){
    cnv_mapping <- shared_pathways %>%
      filter(sample_name == patient_of_interest) %>%
      ungroup() %>%
      mutate(pathway = paste0(pathway,' (', direction,')')) %>%
      dplyr::select(sample_name, pathway, genes, direction, comparison) %>%
      separate_rows(genes) %>%
      filter(genes != 1) %>%
      group_by(sample_name, genes, comparison) %>%
      summarise(pathway = toString(pathway)) %>%
      inner_join(filtered_cnv, by = c("genes"="hgnc_symbol"))
  } else {
    cnv_mapping <- data.frame()
  }
  
  return(list(shared_pathways = shared_pathways, cnv_mapping = cnv_mapping))
}
