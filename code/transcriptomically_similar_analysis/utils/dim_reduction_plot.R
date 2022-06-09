#############################
# Function for UMAP plotting
#############################

# Top 10000 most variable genes
dim_reduction_plot <- function(umap_embedding, clindata, patient_of_interest) {
  set.seed(100)
  
  # add clinical data
  dat <- umap_embedding %>%
    rownames_to_column('Kids_First_Biospecimen_ID')  %>%
    inner_join(clindata %>%
                 rownames_to_column('Kids_First_Biospecimen_ID'), by = 'Kids_First_Biospecimen_ID') 
  
  # shape
  dat$label <- ifelse(dat$Kids_First_Biospecimen_ID == patient_of_interest, "patient_of_interest", dat$cohort)
  dat$label <- as.factor(dat$label)
  dat$label <- relevel(dat$label, ref = "patient_of_interest")
  
  # plot
  p <- ggplot(dat, aes(UMAP1, UMAP2,
                       color = short_histology,
                       shape = label,
                       text = paste0("Sample:", Kids_First_Biospecimen_ID,
                                     "\nShort_histology:", short_histology,
                                     "\nBroad_histology:", broad_histology,
                                     "\nPathology_diagnosis:", pathology_diagnosis,
                                     "\nIntegrated_diagnosis:", integrated_diagnosis,
                                     "\nMolecular Subtype:", molecular_subtype))) +
    geom_jitter(width = 0.5, height = 0.5, aes(size = 5, alpha = 0.8)) +
    theme_bw() + ggtitle("UMAP clustering") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title=element_text(size = 8),
          legend.text=element_text(size = 8)) +
    guides(shape = "none", color = "none", size = "none", alpha = "none")
  return(p)
}
