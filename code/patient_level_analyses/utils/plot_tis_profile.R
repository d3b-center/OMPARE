library(tidyverse)
library(dplyr)
library(ggplot2)

# function to take tis input data and plot
plot_tis_profile <- function(tis_input, score, sampleInfo){
  
  # use sum or average
  patient_score_sum <- tis_input %>%
    filter(subject_id == sampleInfo$subjectID) %>%
    .$score_sum
  patient_score_avg <- tis_input %>%
    filter(subject_id == sampleInfo$subjectID) %>%
    .$score_avg
  tis_input <- tis_input %>%
    filter(subject_id != sampleInfo$subjectID)
  if(score == 'sum'){
    tis_input <- tis_input %>% 
      mutate(score = score_sum)
    yint <- patient_score_sum
    ylab <- "TIS Signature Score (Sum)"
  } else {
    tis_input <- tis_input %>% 
      mutate(score = score_avg)
    yint <- patient_score_avg
    ylab <- "TIS Signature Score (Avg.)"
  }
  
  # order
  disease.order <- tis_input %>%
    group_by(Type, disease) %>%
    summarise(median = median(score), count = n()) %>%
    arrange(desc(median)) %>%
    .$disease
  tis_input$disease <- factor(tis_input$disease, levels = disease.order)
  
  # remove subject identifier as saving the plot takes lots of space
  tis_input <- tis_input %>%
    dplyr::select(-c(score_sum, score_avg, subject_id))
  
  # plot
  p <- ggplot(tis_input, aes(disease, score, fill = Type)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values = c("blue", "red")) + 
    xlab("Disease") + ylab(ylab) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    geom_hline(yintercept = yint, linetype = 2, color = 'gray30') +
    annotate("text", x = 40, y = max(tis_input$score) - 1, 
             label = "- - - Patient SigScore", size = 4, 
             fontface = 'italic', color = "gray30")
  return(p)
}