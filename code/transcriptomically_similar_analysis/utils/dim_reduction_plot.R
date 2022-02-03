#############################
# Function for UMAP plotting
#############################

# Top 10000 most variable genes
dim_reduction_plot <- function(dat, clindata, study, patient, title) {
  set.seed(100)
  
  # add clinical data
  colnames(dat)[1:2] <- c("Dim1", "Dim2")
  dat <- dat %>%
    rownames_to_column('tmp')  %>%
    inner_join(clindata %>%
                 rownames_to_column('tmp'), by = 'tmp') %>%
    column_to_rownames('tmp')
  
  # shape
  dat$label <- ifelse(dat$subject_id == patient, 'patient_of_interest', dat$study_id)
  dat$label <- as.factor(dat$label)
  dat$label <- relevel(dat$label, ref = 'patient_of_interest')
  
  # plot
  p <- ggplot(dat, aes(Dim1, Dim2,
                       color = short_histology,
                       shape = label,
                       text = paste0("Sample:", subject_id,
                                     "\nShort_histology:", short_histology,
                                     "\nBroad_histology:", broad_histology,
                                     "\nPathology_diagnosis:", pathology_diagnosis,
                                     "\nIntegrated_diagnosis:", integrated_diagnosis))) +
    geom_jitter(width = 0.5, height = 0.5, size = 2, alpha = 0.8) +
    theme_bw() + ggtitle(title) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)) +
    guides(shape = "none", color = "none", size = "none", alpha = "none")
  return(p)
}
