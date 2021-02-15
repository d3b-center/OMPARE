# Author: Komal S. Rathi
# Date: 04/25/2020
# Function: CNV + Expression heatmap for PBTA + PNOC008

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# source functions
source(file.path(patient_level_analyses_utils, 'pubTheme.R'))
source(file.path(patient_level_analyses_utils, 'get_zscore.R'))
source(file.path(patient_level_analyses_utils, 'quiet.R'))
source(file.path(patient_level_analyses_utils, 'create_copy_number.R'))
source(file.path(patient_level_analyses_utils, 'batch_correct.R'))

# merge cnv pnoc008
merge_cnv_pnoc <- function(cnvData, gene.list){
  sample_name <- gsub(".*PNOC", "PNOC", cnvData)
  sample_name <- gsub('/.*', '', sample_name)
  cnvData <- data.table::fread(cnvData, header = T, check.names = T)
  
  cnvData <- cnvData %>%
    mutate(chr = as.character(chr))

  # map coordinates to gene symbol
  cnvOut <- create_copy_number(cnvData = cnvData, ploidy = NULL) 
  cnvOut <- cnvOut %>%
    filter(hgnc_symbol %in% gene.list) %>% # filter to gene list
    mutate(sample_name = sample_name) %>% # add PNOC008 patient id
    mutate(copy.number = ifelse(copy.number >=5, 5, copy.number)) # anything > 4 is a deep amplification
}

# merge cnv pbta
merge_cnv_pbta <- function(cnvData, gene.list){
  sample_name <- unique(cnvData$sample_id)
  cnvData$Kids_First_Biospecimen_ID <- NULL
  ploidy <- unique(cnvData$tumor_ploidy)
  
  # map coordinates to gene symbol
  cnvOut <- create_copy_number(cnvData = cnvData, ploidy = ploidy) 
  cnvOut <- cnvOut %>%
    filter(hgnc_symbol %in% gene.list) %>% # filter to gene list
    mutate(sample_name = sample_name) %>% # add PNOC008 patient id
    mutate(copy.number = ifelse(copy.number >=5, 5, copy.number)) # anything > 4 is a deep amplification
}

create_heatmap <- function(fname, genelist, plot.layout = "h"){
  
  genelist <- unique(genelist)
  
  ## Expression
  # PNOC008 clinical
  pnoc.clin <- pnoc008_clinical
  pnoc.clin <- pnoc.clin %>%
    mutate(disease = "HGG",
           disease_subtype = tumorType,
           sample_id = subjectID) %>%
    dplyr::select(subjectID, sample_id, disease, disease_subtype, sex, ethnicity, library_name)
  
  # PNOC008 mRNA
  pnoc.expr <- pnoc008_tpm
  
  # PBTA clinical
  pbta.clin <- read.delim(file.path(ref_dir, 'PBTA', 'pbta-histologies.tsv'))
  pbta.clin <- pbta.clin %>%
    filter(short_histology == "HGAT",
           experimental_strategy %in% c("WGS", "RNA-Seq")) %>%
    mutate(disease = "HGG", 
           disease_subtype = pathology_diagnosis,
           subjectID = Kids_First_Biospecimen_ID,
           sex = germline_sex_estimate,
           library_name = RNA_library) %>%
    dplyr::select(subjectID, sample_id, disease, disease_subtype, sex, ethnicity, experimental_strategy, library_name)
  
  # sample ids with unique WGS + RNA-seq mapping (n = 48)
  sids <- pbta.clin %>%
    group_by(sample_id, experimental_strategy) %>%
    summarise(count = n()) %>%
    group_by(sample_id) %>%
    mutate(sum = sum(count)) %>%
    filter(sum == 2  & count  == 1)
  pbta.clin <- pbta.clin %>%
    filter(sample_id %in% sids$sample_id)
  
  # separate RNA and WGS clinical file
  pbta.rna.clin <- pbta.clin %>%
    filter(experimental_strategy == "RNA-Seq") %>%
    column_to_rownames("subjectID")
  pbta.cnv.clin <-  pbta.clin %>%
    filter(experimental_strategy == "WGS") %>%
    column_to_rownames("subjectID")
  
  # PBTA HGG mRNA expression (n = 186)
  pbta.expr <- pbta_full_tpm
  rna.sids <- intersect(colnames(pbta.expr), rownames(pbta.rna.clin))
  pbta.rna.clin  <- pbta.rna.clin[rna.sids,]
  pbta.expr <- pbta.expr[,rna.sids]
  if(identical(colnames(pbta.expr), rownames(pbta.rna.clin))){
    colnames(pbta.expr) <- pbta.rna.clin$sample_id
  }
  
  # combine PBTA + PNOC008 expression matrix (protein coding only)
  expr <- pnoc.expr %>%
    rownames_to_column("gene_symbol") %>%
    inner_join(pbta.expr %>%
                rownames_to_column("gene_symbol"), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol") 
  
  # now combine clinical files for heatmap
  pbta.clin <- pbta.rna.clin %>%
    as.tibble() %>%
    mutate(subjectID = sample_id,
           study_id = "PBTA") %>%
    column_to_rownames("subjectID") %>%
    dplyr::select(-c(experimental_strategy))
  pnoc.clin <- pnoc.clin %>%
    as.tibble() %>% 
    mutate(study_id = "PNOC008") %>%
    column_to_rownames("subjectID")
  clin <- rbind(pbta.clin, pnoc.clin)
  
  # batch correct
  clin$batch <- paste0(clin$study_id, '_', clin$library_name)
  expr <- expr[,rownames(clin)]
  expr <- quiet(batch.correct(mat = expr, clin = clin))
  
  # subset to genelist of interest
  genelist.expr <- expr %>%
    as.data.frame() %>%
    rownames_to_column("gene_symbol") %>%
    filter(gene_symbol %in% genelist) %>%
    column_to_rownames("gene_symbol")
  
  # convert TPM matrix to z-score matrix
  genelist.expr <- as.data.frame(t(apply(genelist.expr, FUN = get_zscore, MARGIN = 1)))
  
  ## Copy number
  # PBTA
  pbta_cnv <- pbta_cnv %>%
    inner_join(pbta.cnv.clin %>% 
                 dplyr::select(sample_id) %>% 
                 rownames_to_column("subjectID"), by = c("Kids_First_Biospecimen_ID" = "subjectID"))
  
  # subset to genelist
  pbta.cnv.genelist <- plyr::ddply(.data = pbta_cnv, 
                                   .variables = 'Kids_First_Biospecimen_ID', 
                                   .fun = function(x) merge_cnv_pbta(cnvData = x, gene.list = genelist))
  
  # PNOC008
  cnv.files <- list.files(path = getwd(), pattern = "*.CNVs.p.value.txt", recursive = TRUE, full.names = T)
  pnoc.cnv.genelist <- lapply(cnv.files, FUN = function(x) merge_cnv_pnoc(cnvData = x, gene.list = genelist))
  pnoc.cnv.genelist <- data.table::rbindlist(pnoc.cnv.genelist)
  
  # merge PBTA and PNOC
  genelist.cnv <- rbind(pbta.cnv.genelist %>%
                          dplyr::select(sample_name, copy.number, hgnc_symbol), 
                        pnoc.cnv.genelist %>%
                          dplyr::select(sample_name, copy.number, hgnc_symbol))
  
  # convert to matrix
  genelist.cnv <- genelist.cnv %>%
    spread(sample_name, copy.number) %>%
    column_to_rownames("hgnc_symbol")
  
  # only keep NANT sample for PNOC008-5
  genelist.cnv <- genelist.cnv[,grep('CHOP', colnames(genelist.cnv), invert = T)]
  colnames(genelist.cnv)  <- gsub("-NANT", "", colnames(genelist.cnv))
  
  # plot with ComplexHeatmap
  # make cnv matrix consistent with expression using common genes and samples
  expr.mat <- genelist.expr
  cnv.mat <- genelist.cnv
  common.samples <- intersect(colnames(expr.mat), colnames(cnv.mat))
  common.genes <- intersect(rownames(expr.mat), rownames(cnv.mat))
  cnv.mat <- cnv.mat[common.genes, common.samples]
  expr.mat <- expr.mat[common.genes, common.samples]
  clin <- clin %>%
    filter(sample_id %in% common.samples)
  
  # define color scheme
  # disease_subtype: assign ggplot2 default colors 
  n = length(unique(clin$disease_subtype))
  subtype_cols = gg_color_hue(n)
  names(subtype_cols) <- unique(clin$disease_subtype)
  # gender
  gender_cols <- c("Male" = "steelblue", "Female" = "palevioletred1")
  
  # set heatmap options
  ht_opt(heatmap_column_title_gp = gpar(fontsize = 18, fontface = "bold"),
         heatmap_column_names_gp = gpar(fontsize = 18, fontface = "bold"),
         heatmap_row_names_gp = gpar(fontsize = 18, fontface = "bold"),
         legend_border = "black",
         heatmap_border = TRUE,
         annotation_border = TRUE)
  
  # set heatmap legend params
  heatmap_legend_param_global <- list(color_bar = 'continuous', 
                                      direction = 'horizontal',
                                      grid_height = unit(0.7, "cm"), 
                                      grid_width = unit(3, "mm"),
                                      title_gp = gpar(fontsize = 18, fontface = "bold"),
                                      labels_gp = gpar(fontsize = 18, fontface = "bold"))
  
  # set annotation legend params
  # subtype
  subtype_annotation_legend_params <- list(nrow = 10, 
                                           grid_height = unit(0.7, "cm"), 
                                           grid_width = unit(3, "mm"),
                                           title_gp = gpar(fontsize = 18, fontface = "bold"),
                                           labels_gp = gpar(fontsize = 18, fontface = "bold"))
  # gender
  gender_annotation_legend_params <- list(nrow = 4, 
                                          grid_height = unit(0.7, "cm"), 
                                          grid_width = unit(3, "mm"),
                                          title_gp = gpar(fontsize = 18, fontface = "bold"),
                                          labels_gp = gpar(fontsize = 18, fontface = "bold"))
  
  
  png(filename = fname, height = 20, width = 40, units = "in", res = 300)
  # horizontal or vertical layout
  if(plot.layout == "v"){
    # create topannotation
    ha1 = HeatmapAnnotation(Subtype = clin$disease_subtype, 
                            Gender = clin$sex,
                            col = list(Gender = gender_cols,
                                       Subtype = subtype_cols),
                            annotation_legend_param = list(
                              Subtype = subtype_annotation_legend_params,
                              Gender = gender_annotation_legend_params))
    
    ht1 <- Heatmap(as.matrix(expr.mat), cluster_rows = FALSE, 
                   top_annotation = ha1,
                   name = "RNA",  row_title = "RNA",
                   heatmap_legend_param = heatmap_legend_param_global)
    ht2 <- Heatmap(as.matrix(cnv.mat), cluster_rows = FALSE, cluster_columns = FALSE,
                   top_annotation = ha1,
                   name = "CNV", row_title = "CNV (WXS)",
                   heatmap_legend_param = heatmap_legend_param_global)
    draw(ht1 %v% ht2, heatmap_legend_side = "right", annotation_legend_side = "right")
  } else if(plot.layout == "h"){
    # create row annotation
    ha1 = rowAnnotation(Subtype = clin$disease_subtype, 
                        Gender = clin$sex,
                        col = list(Gender = gender_cols,
                                   Subtype = subtype_cols),
                        annotation_legend_param = list(
                          Subtype = subtype_annotation_legend_params,
                          Gender = gender_annotation_legend_params))
    
    ht1 <- Heatmap(t(expr.mat), cluster_columns = FALSE, 
                   left_annotation = ha1,
                   name = "RNA", column_title = "RNA", 
                   heatmap_legend_param = heatmap_legend_param_global)
    
    ht2 <- Heatmap(t(cnv.mat), cluster_rows = FALSE, cluster_columns = FALSE,
                   name = "CNV", column_title = "CNV (WXS)",
                   heatmap_legend_param = heatmap_legend_param_global)
    draw(ht1 + ht2, heatmap_legend_side = "right", padding = unit(c(2, 10, 2, 2), "mm"))
  }
  dev.off()
}