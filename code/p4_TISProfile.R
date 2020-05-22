# Author: Komal S. Rathi
# Date:  05/11/2020
# Function:  TCGA + OpenPBTA + PNOC Quantile Norm -> Avg.

tisProfile <- function(fname, score){
  
  if(!file.exists(fname)){
    # TCGA
    tcga <- readRDS("data/Reference/TCGA/TCGA_matrix_counts.RDS")
    
    # PBTA (stranded)
    pbta.stranded <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-counts-collapsed.stranded.rds')
    
    # PNOC008 expression
    pnoc008 <- expData.counts[,sampleInfo$subjectID, drop=FALSE]

    # read TIS signature
    tis <- read.delim('data/Reference/TIS_geneset.txt', stringsAsFactors = F)
    
    # merge on common genes from TIS signature
    common.genes <- intersect(intersect(rownames(tcga), rownames(pbta.stranded)), rownames(pnoc008))
    tcga <- tcga[common.genes,]
    pbta.stranded <- pbta.stranded[common.genes,]
    pnoc008 <- pnoc008[common.genes, , drop = FALSE]
    total <- cbind(tcga, pbta.stranded, pnoc008)
    
    # now read meta data
    tcga.meta <- readRDS('data/Reference/TCGA/TCGA_meta.RDS')
    tcga.meta <- tcga.meta %>%
      rownames_to_column("sample_id") %>%
      mutate(Type = "Adult") %>%
      dplyr::select(sample_id, disease, Type)
    pbta.meta <- read.delim('data/Reference/PBTA/pbta-histologies.tsv')
    pbta.meta <- pbta.meta %>%
      filter(experimental_strategy  == "RNA-Seq", RNA_library == "stranded") %>%
      mutate(sample_id = Kids_First_Biospecimen_ID, 
             disease = short_histology,
             Type = "Pediatric") %>%
      dplyr::select(sample_id, disease, Type)
    pnoc.meta <- data.frame(sample_id =  sampleInfo$subjectID, disease = "HGAT", Type = "Pediatric")
    total.meta <- rbind(tcga.meta, pbta.meta, pnoc.meta)
    
    # count number of samples per hist and only keep >= 20
    total.meta <- total.meta  %>%
      group_by(disease) %>%
      mutate(n = n()) %>%
      filter(n >= 20) %>%
      dplyr::select(-c(n))
    total.meta <- total.meta[order(total.meta$sample_id),] # order rows
    total <- total[,colnames(total) %in% total.meta$sample_id]
    total <- total[,order(colnames(total))] # order columns
    
    # quantile normalize 
    normalize.mat <- function(mat, meta, method, genelist){
      if(method == "voom"){
        # create design
        var <- factor(meta[,'Type'])
        design <- model.matrix(~0+var)
        colnames(design) <- levels(var)
        rownames(design) <- meta$sample_id
        
        # voom normalize
        v <- voom(counts = mat, design = design, plot = FALSE, normalize.method = "quantile")
        total.norm <- v$E
      }  else {
        mat <-  data.matrix(log2(mat + 1))
        total.norm <- normalize.quantiles(mat)
        rownames(total.norm) <- rownames(mat)
        colnames(total.norm) <- colnames(mat)
        total.norm <- as.data.frame(total.norm)
      }
      
      # filter by genelist and format
      total.norm <- total.norm[rownames(total.norm) %in% genelist,]
      total.sums <- colSums(total.norm)
      total.avg <- colMeans(total.norm)
      total.norm <- data.frame(sample_barcode = names(total.sums), scoreSum = total.sums, scoreAvg = total.avg)
      total.norm <- total.norm[order(total.norm$sample_barcode),]
      
      # merge with meta file
      total.norm <- meta %>% 
        dplyr::select(disease, Type) %>% 
        cbind(total.norm) %>%
        as.data.frame()
      
      return(total.norm)
    }
    total <- normalize.mat(mat = total, meta = total.meta, genelist = tis$Genes, method = "quantile")
    write.table(total, file = fname, sep = "\t", quote = F, row.names = F)
  } else {
    total  <- read.delim(fname, stringsAsFactors = F)
  }
  # use sum or average
  pnoc008.scoreSum <- total[grep('PNOC008', total$sample_barcode),'scoreSum']
  pnoc008.scoreAvg <- total[grep('PNOC008', total$sample_barcode),'scoreAvg']
  total <- total[grep('PNOC008', total$sample_barcode, invert = T),]
  if(score == 'sum'){
    total <- total %>% 
      mutate(score = scoreSum)
    yint <- pnoc008.scoreSum
    ylab <- "TIS Signature Score (Sum)"
  } else {
    total <- total %>% 
      mutate(score = scoreAvg)
    yint <- pnoc008.scoreAvg
    ylab <- "TIS Signature Score (Avg.)"
  }
  
  # order
  disease.order <- total %>%
    group_by(Type, disease) %>%
    summarise(median = median(score), count = n()) %>%
    arrange(desc(median)) %>%
    .$disease
  total$disease <- factor(total$disease, levels = disease.order)
  
  # plot
  p <- ggplot(total, aes(disease, score, fill = Type)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values = c("blue", "red")) + 
    xlab("Disease") + ylab(ylab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = yint, linetype = 2, color = 'gray30') +
    annotate("text", x = 40, y = max(total$score) - 1, 
             label = "- - - Patient SigScore", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}
