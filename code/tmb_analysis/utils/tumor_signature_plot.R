suppressPackageStartupMessages({
  library(data.table)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(decompTumor2Sig)
  library(ggplot2)
})

tumor_signature_plot <- function(patient_dir, signatures) {
  
  # have to write out to read in
  tmb_analysis_output_dir <- file.path(patient_dir, "output", "tmb_analysis")
  fname <- file.path(tmb_analysis_output_dir, "consensus_mpf_output.txt")
  
  # read "all" for mutational signatures
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = patient_dir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  mutFiles <- grep('consensus', mutFiles, invert = TRUE, value = TRUE)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData.all <- data.table::rbindlist(mutFiles, fill = TRUE)
    mutData.all <- as.data.frame(mutData.all)
    mutData.all <- unique(mutData.all)
  }
  
  mpfData <- mutData.all %>% 
    dplyr::select(Chromosome, Start_Position, Match_Norm_Seq_Allele1, Tumor_Seq_Allele2)
  mpfData <- data.frame("Patient", mpfData)
  write.table(mpfData, fname, sep = "\t", row.names = F, quote = F, col.names = F)
  
  # load the reference genome and the transcript annotation database
  refGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  transcriptAnno <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Read in the genome
  genomes <- readGenomesFromMPF(fname, numBases = 3, type = "Alexandrov", trDir = F, refGenome = refGenome, verbose=F)
  exposure <- decomposeTumorGenomes(genomes, signatures)$Patient
  exposure <- data.frame('Signature' = names(exposure), 'Value' = exposure)
  exposure$Signature <- gsub("[.]", "-", exposure$Signature)
  exposure$Signature <- factor(exposure$Signature, exposure[order(exposure$Value),'Signature'])
  
  p <- ggplot(exposure, aes(Signature, Value, fill = Value)) + 
    guides(fill = F) +
    geom_bar(stat="identity") + theme_bw() + coord_flip() + 
    xlab("Mutational Signature") + ylab("Exposures (percent contribution)")
  return(p)
}
