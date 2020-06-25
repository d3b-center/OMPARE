# Author: Komal S. Rathi
# Date: 04/25/2020
# Function: CNV + Expression heatmap for PBTA + PNOC008

# function to calculate z-score
getZ <- function(x) {
  x <- log2(x+1)
  out <- (x-mean(x))/sd(x)
  return(out)
}

# function to read cnv, filter by genes and merge
merge.cnv <- function(cnvData, gene.list){
  if(length(grep("Kids_First_Biospecimen_ID", colnames(cnvData))) == 1){
    # PBTA
    sample_name <- unique(cnvData$sample_id)
    cnvData$Kids_First_Biospecimen_ID <- NULL
    ploidy <- unique(cnvData$tumor_ploidy)
  } else {
    # PNOC
    sample_name <- gsub(".*PNOC", "PNOC", cnvData)
    sample_name <- gsub('/.*', '', sample_name)
    sample_name <- gsub('-[0]+', '-', sample_name)
    cnvData <- data.table::fread(cnvData, header = T, check.names = T)
    ploidy <- NULL
  }
  cnvData <- cnvData %>% 
    dplyr::select(chr, start, end, copy.number, 
                  status, WilcoxonRankSumTestPvalue) %>%
    as.data.frame()
  cnvOut <- createCopyNumber(cnvData = cnvData, ploidy = ploidy) # map coordinates to gene symbol
  cnvOut <- cnvOut %>%
    filter(Gene %in% gene.list) %>% # filter to gene list
    mutate(sample_name = sample_name) %>% # add PNOC008 patient id
    mutate(CNA = ifelse(CNA >=5, 5, CNA)) # anything > 4 is a deep amplification
  return(cnvOut)
}

create.heatmap <- function(fname, genelist, plot.layout = "h"){
  
  genelist <- unique(genelist)
  
  ## Expression
  # PNOC008 clinical
  pnoc.clin <- pnoc008.clinData
  pnoc.clin <- pnoc.clin %>%
    mutate(disease = "HGG",
           disease_subtype = tumorType,
           sample_id = subjectID) %>%
    dplyr::select(subjectID, sample_id, disease, disease_subtype, sex, ethnicity)
  
  # PNOC008 mRNA
  pnoc.expr <- pnoc008.data
 
  # PBTA clinical
  pbta.clin <- read.delim('data/Reference/PBTA/pbta-histologies.tsv')
  pbta.clin <- pbta.clin %>%
    filter(short_histology == "HGAT",
           experimental_strategy %in% c("WGS", "RNA-Seq")) %>%
    mutate(disease = "HGG", 
           disease_subtype = pathology_diagnosis,
           subjectID = Kids_First_Biospecimen_ID,
           sex = germline_sex_estimate) %>%
    dplyr::select(subjectID, sample_id, disease, disease_subtype, sex, ethnicity, experimental_strategy)
 
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
  pbta.expr <- pbta.full
  rna.sids <- intersect(colnames(pbta.expr), rownames(pbta.rna.clin))
  pbta.rna.clin  <- pbta.rna.clin[rna.sids,]
  pbta.expr <- pbta.expr[,rna.sids]
  if(identical(colnames(pbta.expr), rownames(pbta.rna.clin))){
    colnames(pbta.expr) <- pbta.rna.clin$sample_id
  }
  
  # combine PBTA + PNOC008 expression matrix
  expr <- pnoc.expr %>%
    rownames_to_column("gene_symbol") %>%
    full_join(pbta.expr %>%
                rownames_to_column("gene_symbol"), by = "gene_symbol") %>%
    column_to_rownames("gene_symbol") 
  
  
  # subset to genelist of interest
  genelist.expr <- expr %>%
    rownames_to_column("gene_symbol") %>%
    filter(gene_symbol %in% genelist) %>%
    column_to_rownames("gene_symbol")
  
  # convert TPM matrix to z-score matrix
  genelist.expr <- as.data.frame(t(apply(genelist.expr, FUN = getZ, MARGIN = 1)))
  
  ## Copy number
  # PBTA
  pbta.cnv <- pbta.cnv %>%
    inner_join(pbta.cnv.clin %>% 
                 dplyr::select(sample_id) %>% 
                 rownames_to_column("subjectID"), by = c("Kids_First_Biospecimen_ID" = "subjectID"))
  
  # subset to genelist
  pbta.cnv.genelist <- plyr::ddply(.data = pbta.cnv, 
                              .variables = 'Kids_First_Biospecimen_ID', 
                              .fun = function(x) merge.cnv(cnvData = x, gene.list = genelist))
  
  # PNOC008
  cnv.files <- list.files(path = getwd(), pattern = "*.CNVs.p.value.txt", recursive = TRUE, full.names = T)
  pnoc.cnv.genelist <- lapply(cnv.files, FUN = function(x) merge.cnv(cnvData = x, gene.list = genelist))
  pnoc.cnv.genelist <- data.table::rbindlist(pnoc.cnv.genelist)
  
  # merge PBTA and PNOC
  genelist.cnv <- rbind(pbta.cnv.genelist %>%
                     dplyr::select(-c(Kids_First_Biospecimen_ID)), pnoc.cnv.genelist)
  
  # convert to matrix
  genelist.cnv <- genelist.cnv %>%
    dplyr::select(-c(Pvalue, Status)) %>%
    spread(sample_name, CNA) %>%
    column_to_rownames("Gene")
  
  # only keep CHOP sample for PNOC008-5
  genelist.cnv <- genelist.cnv[,grep('NANT', colnames(genelist.cnv), invert = T)]
  colnames(genelist.cnv)  <- gsub("-CHOP", "", colnames(genelist.cnv))
  
  # now combine clinical files for heatmap
  pbta.clin <- pbta.rna.clin %>%
    as.tibble() %>%
    mutate(subjectID = sample_id) %>%
    column_to_rownames("subjectID") %>%
    dplyr::select(-c(experimental_strategy))
  pnoc.clin <- pnoc.clin %>% 
    as.tibble() %>% 
    column_to_rownames("subjectID")
  clin <- rbind(pbta.clin, pnoc.clin)
  
  # plot with ComplexHeatmap
  # make cnv matrix consistent with expression
  expr.mat <- genelist.expr
  cnv.mat <- genelist.cnv
  cnv.mat <- cnv.mat[rownames(expr.mat),colnames(expr.mat)]
  
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
  
  
  png(filename = fname, height = 15, width = 35, units = "in", res = 300)
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