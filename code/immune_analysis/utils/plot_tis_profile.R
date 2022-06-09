suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
})

# function to take tis input data and plot
plot_tis_profile <- function(tis_profile_output, score_type){
  
  # use sum or average
  patient_score_sum <- tis_profile_output %>%
    filter(type == "Patient") %>%
    .$score_sum
  patient_score_avg <- tis_profile_output %>%
    filter(type == "Patient") %>%
    .$score_avg
  tis_profile_output <- tis_profile_output %>%
    filter(type != "Patient")
  if(score_type == 'sum'){
    tis_profile_output <- tis_profile_output %>% 
      dplyr::mutate(score = score_sum)
    yint <- patient_score_sum
    ylab <- "TIS Signature Score (Sum)"
  } else {
    tis_profile_output <- tis_profile_output %>% 
      dplyr::mutate(score = score_avg)
    yint <- patient_score_avg
    ylab <- "TIS Signature Score (Avg.)"
  }
  
  # order
  tis_profile_output <- tis_profile_output %>%
    mutate(label = paste0(cohort, ": ", short_histology)) 
  disease.order <- tis_profile_output %>%
    group_by(type, label) %>%
    dplyr::summarise(median = median(score), count = n()) %>%
    arrange(desc(median)) %>%
    .$label
  tis_profile_output$label <- factor(tis_profile_output$label, levels = disease.order)
  
  # remove subject identifier as saving the plot takes lots of space
  tis_profile_output <- tis_profile_output %>%
    dplyr::select(-c(score_sum, score_avg, subject_id))
  
  # plot
  p <- ggplot(tis_profile_output, aes(label, score, fill = type)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values = c("Adult" = "blue", "Pediatric" = "red")) + 
    xlab("Disease") + ylab(ylab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = yint, linetype = 2, color = 'gray30') +
    annotate("text", x = length(unique(tis_profile_output$label))/2, y = max(tis_profile_output$score), 
             label = "- - - Patient SigScore", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}
