# Author: Komal S. Rathi
# Date:  05/11/2020
# Function:  TCGA + OpenPBTA + PNOC GSVA

tisProfile <- function(fname, score){
  
  if(!file.exists(fname)){
    # TCGA
    tcga <- readRDS("data/Reference/TCGA/TCGA_matrix_FPKM.RDS")
    
    # PBTA (stranded)
    pbta.stranded <- readRDS('data/Reference/PBTA/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds')
    
    # PNOC008 expression
    pnoc008 <- expData[,sampleInfo$subjectID, drop=FALSE]
    
    # read TIS signature
    tis <- read.delim('data/Reference/TIS_geneset.txt', stringsAsFactors = F)
    
    # merge on common genes from TIS signature
    common.genes <- intersect(intersect(rownames(tcga), rownames(pbta.stranded)), rownames(pnoc008))
    common.genes <- intersect(tis$Genes, common.genes)
    tcga <- tcga[common.genes,]
    pbta.stranded <- pbta.stranded[common.genes,]
    pnoc008 <- pnoc008[common.genes, , drop = FALSE]
    total <- cbind(tcga, pbta.stranded, pnoc008)
    
    # z-score data
    total <- t(apply(total, 1, getZ))
    total.sums <- colSums(total)
    total.avg <- colMeans(total)
    total <- data.frame(sample_barcode = names(total.sums), scoreSum = total.sums, scoreAvg = total.avg)
    total <- total[order(total$sample_barcode),]
    
    # now read meta data
    tcga.meta <- readRDS('data/Reference/TCGA/TCGA_meta.RDS')
    tcga.meta <- tcga.meta %>%
      mutate(Type = "Adult") %>%
      dplyr::select(sample_barcode, disease, Type)
    pbta.meta <- read.delim('data/Reference/PBTA/pbta-histologies.tsv')
    pbta.meta <- pbta.meta %>%
      filter(experimental_strategy  == "RNA-Seq", RNA_library == "stranded") %>%
      mutate(sample_barcode = Kids_First_Biospecimen_ID, 
             disease = short_histology,
             Type = "Pediatric") %>%
      dplyr::select(sample_barcode, disease, Type)
    pnoc.meta <- data.frame(sample_barcode =  sampleInfo$subjectID, disease = "", Type = "")
    total.meta <- rbind(tcga.meta, pbta.meta, pnoc.meta)
    total.meta <- total.meta[order(total.meta$sample_barcode),]
    
    # merge both
    total <- total.meta %>% 
      dplyr::select(disease, Type) %>% 
      cbind(total)
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
    ylab <- "TIS Signature Score (Sum of Z-scores)"
  } else {
    total <- total %>% 
      mutate(score = scoreAvg)
    yint <- pnoc008.scoreAvg
    ylab <- "TIS Signature Score (Avg of Z-scores)"
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
    annotate("text", x = 50, y = max(total$score) - 1, 
             label = "- - - Patient SigScore", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}
