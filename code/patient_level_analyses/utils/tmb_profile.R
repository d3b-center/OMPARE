tmb.calculate <- function(myTMB = tmb_bed_file, var_class = c('Missense_Mutation', 
                                                              'Nonsense_Mutation', 
                                                              'Frame_Shift_Del', 
                                                              'Frame_Shift_Ins', 
                                                              'In_Frame_Del', 
                                                              'In_Frame_Ins'), 
                          vaf_cutoff = 0.05, var_count = 3, tumor_depth = 25) {
  
  # read mutect2 for TMB profile
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = topDir, pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  mutFiles <- grep("mutect2", mutFiles, value = TRUE)
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData.mutect2 <- data.table::rbindlist(mutFiles)
    mutData.mutect2 <- as.data.frame(mutData.mutect2)
    mutData.mutect2 <- unique(mutData.mutect2)
  }
  
  # filter to nonsense and missense
  myMutData <- mutData.mutect2 %>%
    mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
    filter(Variant_Classification %in% var_class,
           t_depth >= tumor_depth,
           vaf <= vaf_cutoff) %>%
  dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
  
  # intersect with bed file
  subject <- with(myTMB, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(myMutData, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res <- data.frame(myMutData[queryHits(res),], myTMB[subjectHits(res),])
  
  # return the number of missense+nonsense overlapping with the bed file
  return(nrow(res))
}


tmb_profile <- function(pedTMBScores = ped_tmb, adultTMBScores = adult_tmb, TMB, tmb_bed_file) {
  
  TMB <- tmb.calculate(myTMB = tmb_bed_file)/TMB
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
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = TMB, linetype = 2, color = 'gray30') +
    annotate("text", x = 35, y = max(tmbScores$TMBscore) - 50, 
             label = "- - - Patient TMB", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}
