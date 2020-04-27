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
  } else {
    # PNOC
    sample_name <- gsub(".*PNOC", "PNOC", cnvData)
    sample_name <- gsub('/.*', '', sample_name)
    sample_name <- gsub('-[0]+', '-', sample_name)
    cnvData <- data.table::fread(cnvData)
  }
  cnvOut <- createCopyNumber(cnvData = cnvData) # map coordinates to gene symbol
  cnvOut <- cnvOut %>%
    filter(Gene %in% gene.list) %>% # filter to gene list
    mutate(sample_name = sample_name) %>% # add PNOC008 patient id
    mutate(CNA = ifelse(CNA >=5, 5, CNA)) # anything > 4 is a deep amplification
  return(cnvOut)
}

create.heatmap <- function(fname, genelist, plot.layout = "h"){
  
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
    filter(integrated_diagnosis == "High-grade glioma",
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
  
  # PBTA HGG mRNA expression (n = 112)
  pbta.expr <- pbta.polya %>%
    rownames_to_column("gene_symbol") %>%
    full_join(pbta.stranded %>%
                rownames_to_column("gene_symbol"), by = 'gene_symbol') %>%
    column_to_rownames("gene_symbol")
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
  
  # convert FPKM matrix to z-score matrix
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
  
  # assign ggplot2 default colors to disease_subtype
  n = length(unique(clin$disease_subtype))
  cols = gg_color_hue(n) 
  names(cols) <- unique(clin$disease_subtype)
  
  # set heatmap options
  ht_opt(heatmap_column_names_gp = gpar(fontsize = 10),
         heatmap_row_names_gp = gpar(fontsize = 10),
         heatmap_column_title_gp = gpar(fontsize = 12),
         legend_border = "black",
         heatmap_border = TRUE,
         annotation_border = TRUE
  )
  
  png(fname, height = 17, width = 20, units = "in", res = 300)
  # horizontal or vertical layout
  if(plot.layout == "v"){
    # create topannotation
    ha1 = HeatmapAnnotation(disease_subtype = clin$disease_subtype, 
                            gender = clin$sex,
                            col = list(gender = c("Male" = "steelblue", "Female" = "palevioletred1"),
                                       disease_subtype = cols),
                            annotation_legend_param = list(
                              disease_subtype = list(direction = "horizontal"),
                              gender = list(direction = "horizontal")))
    ht1 <- Heatmap(as.matrix(expr.mat), cluster_rows = FALSE, top_annotation = ha1,
                   name = "RNA",  row_title = "RNA",
                   heatmap_legend_param= list(color_bar= 'continuous', direction = "horizontal"))
    ht2 <- Heatmap(as.matrix(cnv.mat), cluster_rows = FALSE, cluster_columns = FALSE,
                   top_annotation = ha1,
                   name = "CNV", row_title = "CNV (WXS)",
                   heatmap_legend_param = list(color_bar= 'continuous', direction = "horizontal"))
    draw(ht1 %v% ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  } else if(plot.layout == "h"){
    # create row annotation
    ha1 = rowAnnotation(disease_subtype = clin$disease_subtype, 
                        gender = clin$sex,
                        col = list(gender = c("Male" = "steelblue", "Female" = "palevioletred1"),
                                   disease_subtype = cols),
                        annotation_legend_param = list(
                          disease_subtype = list(direction = "horizontal"),
                          gender = list(direction = "horizontal")))
    ht1 <- Heatmap(t(expr.mat), cluster_columns = FALSE, left_annotation = ha1,
                   name = "RNA",  column_title = "RNA",
                   heatmap_legend_param = list(color_bar= 'continuous', direction = "horizontal"))
    ht2 <- Heatmap(t(cnv.mat), cluster_rows = FALSE, cluster_columns = FALSE,
                   name = "CNV", column_title = "CNV (WXS)",
                   heatmap_legend_param = list(color_bar= 'continuous', direction = "horizontal"))
    draw(ht1 + ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  }
  dev.off()
}