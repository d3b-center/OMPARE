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
  
  # reverse factors to set PNOC to lowest
  # this is for shape
  dat$label <- dat$study_id
  dat$label <- ifelse(dat$sample_barcode == patient, patient, dat$label)
  dat$label <- factor(dat$label, levels = c(study, "PNOC008", patient))
  
  # this is for size
  dat$study_id <- as.factor(dat$study_id)
  dat$study_id <- relevel(dat$study_id, ref = "PNOC008")
  dat$study_id <- fct_rev(dat$study_id)
  # if(study != "PBTA"){
  #   dat$pathology_diagnosis <- dat$short_histology
  #   dat$integrated_diagnosis <- dat$broad_histology
  # }
  
  # plot
  p <- ggplot(dat, aes(Dim1, Dim2,
                       color = pathology_diagnosis,
                       size = study_id,
                       shape = label,
                       text = paste0("Sample:", sample_barcode,
                                     "\nShort_histology:", short_histology,
                                     "\nBroad_histology:", broad_histology,
                                     "\nPathology_diagnosis:", pathology_diagnosis,
                                     "\nIntegrated_diagnosis:", integrated_diagnosis))) +
    geom_jitter(width = 0.5, height = 0.5) +
    theme_bw() + ggtitle(title) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title=element_text(size = 8),
          legend.text=element_text(size = 8)) +
    guides(size = FALSE, shape = F, color = F)
  ggplotly(p, tooltip = "text")
}
