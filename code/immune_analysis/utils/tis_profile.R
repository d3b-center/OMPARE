# tumor inflammation signature
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# reference directories
pbta_dir <- file.path(data_dir, "pbta")
tcga_dir <- file.path(data_dir, "tcga")

tis_profile <- function(patient_clinical, sampleInfo){
  
  # TCGA counts
  tcga <- readRDS(file.path(tcga_dir, "tcga_matrix_counts.rds"))
  
  # PBTA counts  (polyA + stranded count data collapsed to gene symbols)
  pbta.full <- readRDS(file.path(pbta_dir, 'pbta-gene-expression-rsem-counts-collapsed.polya.stranded.rds'))
  
  # PNOC008 expression
  pnoc008 <- expData.counts[,sampleInfo$subjectID, drop=FALSE]
  
  # read TIS signature
  tis <- read.delim(file.path(data_dir, 'TIS_geneset.txt'), stringsAsFactors = F)
  
  # merge on common genes from TIS signature
  common.genes <- intersect(intersect(rownames(tcga), rownames(pbta.full)), rownames(pnoc008))
  tcga <- tcga[common.genes,]
  pbta.full <- pbta.full[common.genes,]
  pnoc008 <- pnoc008[common.genes, , drop = FALSE]
  total <- cbind(tcga, pbta.full, pnoc008)
  
  # now read meta data
  tcga.meta <- readRDS(file.path(tcga_dir, 'tcga_clinical.rds'))
  tcga.meta <- tcga.meta %>%
    rownames_to_column("sample_id") %>%
    mutate(Type = "Adult",
           study_id = "TCGA") %>%
    dplyr::select(sample_id, disease, Type, study_id, library_name)
  pbta.meta <- read.delim(file.path(pbta_dir, 'pbta-histologies.tsv'))
  pbta.meta <- pbta.meta %>%
    filter(experimental_strategy  == "RNA-Seq",
           Kids_First_Biospecimen_ID %in% colnames(pbta.full)) %>%
    mutate(sample_id = Kids_First_Biospecimen_ID, 
           disease = short_histology,
           Type = "Pediatric",
           study_id = "PBTA",
           library_name = RNA_library) %>%
    dplyr::select(sample_id, disease, Type, study_id, library_name)
  pnoc.meta <- patient_clinical %>%
    filter(subjectID == sampleInfo$subjectID) %>%
    mutate(sample_id = sampleInfo$subjectID, 
           disease = "HGAT",
           Type = "Pediatric",
           study_id = study_id, 
           library_name = sampleInfo$library_name) %>%
    dplyr::select(all_of(colnames(pbta.meta)))
  total.meta <- rbind(tcga.meta, pbta.meta, pnoc.meta)
  total.meta$batch <- paste0(total.meta$study_id, '_',  total.meta$library_name)
  
  # count number of samples per hist and only keep >= 20
  total.meta <- total.meta  %>%
    group_by(disease) %>%
    mutate(n = n()) %>%
    filter(n >= 20) %>%
    dplyr::select(-c(n)) %>%
    as.data.frame()
  total.meta <- total.meta[order(total.meta$sample_id),] # order rows
  total <- total[,colnames(total) %in% total.meta$sample_id]
  total <- total[,order(colnames(total))] # order columns
  
  # quantile normalize 
  normalize.mat <- function(mat, meta, method, genelist){
    if(method == "voom"){
      # create design to correct for batch
      var <- factor(meta$batch)
      design <- model.matrix(~0+var)
      colnames(design) <- levels(var)
      rownames(design) <- meta$sample_id
      
      # voom normalize
      v <- voom(counts = mat, design = design, plot = FALSE, normalize.method = "quantile")
      total.norm <- v$E
    }  else {
      mat <-  data.matrix(log2(mat + 1))
      total.norm <- normalize.quantiles(mat, copy = FALSE)
      total.norm <- as.data.frame(total.norm)
    }
    
    # filter by genelist and format
    total.norm <- total.norm[rownames(total.norm) %in% genelist,]
    total.sums <- colSums(total.norm)
    total.avg <- colMeans(total.norm)
    total.norm <- data.frame(subject_id = names(total.sums), score_sum = total.sums, score_avg = total.avg)
    total.norm <- total.norm[order(total.norm$subject_id),]
    
    # merge with meta file
    total.norm <- meta %>% 
      dplyr::select(disease, Type) %>% 
      cbind(total.norm) %>%
      as.data.frame()
    
    return(total.norm)
  }
  # normalize.quantiles is faster than voom
  total <- normalize.mat(mat = total, meta = total.meta, genelist = tis$Genes, method = "quantile") 
  return(total)
}

  
