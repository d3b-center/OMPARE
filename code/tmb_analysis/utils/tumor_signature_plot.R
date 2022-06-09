suppressPackageStartupMessages({
  library(data.table)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(decompTumor2Sig)
  library(ggplot2)
})

tumor_signature_plot <- function(maf_data, signatures, output_dir) {
  
  # have to write out to read in
  fname <- file.path(output_dir, "consensus_mpf_output.txt")
  
  mpfData <- maf_data %>% 
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
    guides(fill = "none") +
    geom_bar(stat="identity") + theme_bw() + coord_flip() + 
    xlab("Mutational Signature") + ylab("Exposures (percent contribution)")
  return(p)
}
