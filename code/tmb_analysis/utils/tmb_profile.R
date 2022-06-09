tmb_profile <- function(patient_dir, pbta_tmb, tcga_in_pbta_tmb, tcga_not_in_pbta_tmb, tmb_bed_file) {
  
  # Friends of Cancer Research Project cutoffs
  var_class = c('Missense_Mutation', 'Nonsense_Mutation',
                'Frame_Shift_Del', 'Frame_Shift_Ins', 
                'In_Frame_Del', 'In_Frame_Ins')
  vaf_cutoff = 0.05
  var_count = 3
  tumor_depth = 25
  
  # calculate the length of the bed file
  bed_length <- 0
  for (i in 1:nrow(tmb_bed_file)) {
    distance <- as.numeric(tmb_bed_file[i,3]) - as.numeric(tmb_bed_file[i,2])
    bed_length <- bed_length + distance
  }
  
  # read mutect2 for TMB profile
  mutect2_maf <- list.files(path = patient_dir, pattern = "mutect2", recursive = TRUE, full.names = T)
  mutect2_maf <- data.table::fread(mutect2_maf)
  
  # apply filters using Friends of Cancer Research Project Standards 
  mutect2_maf <- mutect2_maf %>%
    mutate(vaf = t_alt_count/(t_alt_count+t_ref_count)) %>%
    filter(Variant_Classification %in% var_class,
           t_depth >= tumor_depth,
           vaf >= vaf_cutoff,
           t_alt_count >= var_count) %>%
    dplyr::select(Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
  
  # intersect with bed file
  subject <- with(tmb_bed_file, GRanges(chr, IRanges(start = start, end = end)))
  query <- with(mutect2_maf, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res <- data.frame(mutect2_maf[queryHits(res),], tmb_bed_file[subjectHits(res),])
  
  # number of missense+nonsense overlapping with the bed file
  TMB <- (nrow(res))*1000000/bed_length
  
  pbta_tmb$Type <- "PBTA"
  tcga_in_pbta_tmb$Type <- "TCGA in PBTA"
  tcga_not_in_pbta_tmb$Type <- "TCGA not in PBTA"
  tmbScores <- rbind(pbta_tmb, tcga_in_pbta_tmb, tcga_not_in_pbta_tmb)
  
  # count median per histology
  disease.order <- tmbScores %>%
    group_by(Type, Diseasetype) %>%
    dplyr::summarise(median = median(TMBscore), count = n()) %>%
    filter(count > 2) %>%
    arrange(desc(median)) %>%
    .$Diseasetype
  
  # subset to disease type with > 2 samples
  tmbScores <- tmbScores %>%
    filter(Diseasetype %in% disease.order)
  tmbScores$Diseasetype <- factor(tmbScores$Diseasetype, levels = disease.order)
  
  # plot it
  p <- ggplot(tmbScores, aes(Diseasetype, TMBscore, fill = Type)) + 
    geom_boxplot() + theme_bw() + scale_y_log10(breaks = c(.25, 1, 10, 100, 500)) + 
    scale_fill_manual(values = c("PBTA" = "blue", "TCGA in PBTA" = "orange", "TCGA not in PBTA" = "red")) + 
    xlab("") + ylab("Mutations per MB") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = TMB, linetype = 2, color = 'gray30') +
    annotate("text", x = 35, y = max(tmbScores$TMBscore) - 50, 
             label = "- - - Patient TMB", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}

