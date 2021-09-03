kaplan_meier <- function(all_cor, surv_data) {
  
  # survival with 20 nearest neighbors
  surv_data$group <- surv_data$subject_id %in% all_cor$nearest_neighbor
  surv_data$group <- ifelse(surv_data$group, "Cluster With Patient", "Cluster away from Patient")
  fit <- survival::survfit(Surv(OS_days, OS_status) ~ group, data = surv_data)
  ggsurvplot(fit,  data = surv_data, 
             risk.table = TRUE,
             pval = TRUE,   
             conf.int = F, 
             xlim = c(0,1000),  
             break.time.by = 100, 
             ggtheme = theme_minimal(),
             risk.table.y.text.col = T, 
             risk.table.y.text = FALSE)  
  
}
