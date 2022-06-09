suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

kaplan_meier <- function(nn_table, surv_data) {
  
  # survival with 20 nearest neighbors
  surv_data <- surv_data %>%
    rownames_to_column("Kids_First_Biospecimen_ID")
  surv_data$group <- surv_data$Kids_First_Biospecimen_ID  %in% nn_table$nearest_neighbor
  surv_data$group <- ifelse(surv_data$group, "Cluster With Patient", "Cluster away from Patient")
  surv_data$OS_status <- ifelse(surv_data$OS_status == "DECEASED", 0, 1)
  fit <- survival::survfit(Surv(as.numeric(OS_days), OS_status) ~ group, data = surv_data)
  ggsurvplot(fit,  data = surv_data, 
             risk.table = TRUE,
             pval = TRUE,   
             conf.int = F, 
             ggtheme = theme_minimal(),
             risk.table.y.text.col = T, 
             risk.table.y.text = FALSE)  
  
}
