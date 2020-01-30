############################
# TMB Profile
############################

tmb.calculate <- function(myTMB = TMBFileBED) {
  
  # read mutect2 for TMB profile
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = topDir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  mutFiles <- grep("mutect2", mutFiles, value = TRUE)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData.mutect2 <- data.table::rbindlist(mutFiles)
    mutData.mutect2 <- as.data.frame(mutData.mutect2)
    mutData.mutect2 <- unique(mutData.mutect2)
    # assign("mutData.mutect2", mutData.mutect2, envir = globalenv())
  }
  
  # filter to nonsense and missense
  myMutData <- mutData.mutect2 %>%
    filter(Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation")) %>%
    dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
  
  # intersect with bed file
  subject <- with(myTMB, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(myMutData, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res <- data.frame(myMutData[queryHits(res),], myTMB[subjectHits(res),])
  
  # return the number of missense+nonsense overlapping with the bed file
  return(nrow(res))
}


tmbProfile <- function(pedTMBScores = pedTMB, adultTMBScores = adultTMB, TMB) {
  
  TMB <- tmb.calculate()/TMB
  pedTMBScores$Type <- "Pediatric"
  adultTMBScores$Type <- "Adult"
  tmbScores <- rbind(pedTMBScores, adultTMBScores)
  
  # count median per histology
  disease.order <- tmbScores %>%
    group_by(Type, Diseasetype) %>%
    summarise(median = median(TMBscore), count = n()) %>%
    filter(count > 2) %>%
    arrange(desc(median)) %>%
    .$Diseasetype
  
  # subset to disease type with > 2 samples
  tmbScores <- tmbScores %>%
    filter(Diseasetype %in% disease.order)
  tmbScores$Diseasetype <- factor(tmbScores$Diseasetype, levels = disease.order)
  
  # Plot it
  p <- ggplot(tmbScores, aes(Diseasetype, TMBscore, fill = Type)) + 
    geom_boxplot() + theme_bw() + scale_y_log10(breaks = c(.25, 1, 10, 100, 500)) + 
    scale_fill_manual(values = c("blue", "red")) + 
    xlab("Disease") + ylab("Mutations per MB") + 
    theme(axis.text.x = element_text(angle = -90, hjust = (0))) +
    geom_hline(yintercept = TMB, linetype = 2, color = 'gray30') +
    annotate("text", x = 35, y = max(tmbScores$TMBscore) - 50, 
             label = "- - - Patient TMB", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}
