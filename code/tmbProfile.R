############################
# TMB Profile
############################

tmbProfile <- function(pedTMBScores = pedTMB, adultTMBScores = adultTMB, TMB) {
  
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
    geom_hline(yintercept = TMB, linetype = 2)
  return(p)
}