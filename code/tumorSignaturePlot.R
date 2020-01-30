tumorSignaturePlot <- function() {
  
  # Have to write out to read in?
  fname <- paste0(topDir, '/MutationsMAF/mpfDataFormat.txt')
  if(!file.exists(fname)){
    
    # read "all" for mutational signatures
    somatic.mut.pattern <- '*.maf'
    mutFiles <- list.files(path = topDir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
    mutFiles <- grep('consensus', mutFiles, invert = TRUE, value = TRUE)
    if(length(mutFiles) >= 1){
      mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
      mutData.all <- data.table::rbindlist(mutFiles)
      mutData.all <- as.data.frame(mutData.all)
      mutData.all <- unique(mutData.all)
      # assign("mutData.all", mutData.all, envir = globalenv())
    }
    
    mpfData <- mutData.all[, c("Chromosome", "Start_Position", "Match_Norm_Seq_Allele1", "Tumor_Seq_Allele2")]
    mpfData <- data.frame("Patient", mpfData)
    write.table(mpfData, fname, sep="\t", row.names=F, quote=F, col.names=F)
  }
  
  # load the reference genome and the transcript annotation database
  refGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  transcriptAnno <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Read in the genome
  genomes <- readGenomesFromMPF(fname, numBases=3, type="Alexandrov", trDir=F, refGenome=refGenome, verbose=F)
  
  exposure <- decomposeTumorGenomes(genomes, signatures)[[1]]
  
  exposure <- data.frame(names(exposure), exposure)
  colnames(exposure) <- c("Signature", "Value")
  exposure[,1] <- gsub("\\.", "-", exposure[,1])
  exposure[,1] <- factor(exposure[,1], levels=exposure[,1])
  p <- ggplot(exposure, aes(Signature, Value))+geom_bar(stat="identity")+theme_bw()+coord_flip()
  p <- p+xlab("Mutational Signature")+ylab("Exposures (percent contribution)")
  return(p)
  
}