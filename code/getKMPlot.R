getKMPlot <- function(allCor, survData) {
  
  # survival with 20 nearest neighbors
  survData$group <- survData$sample_barcode %in% allCor$nearest_neighbor
  survData$group <- ifelse(survData$group, "Cluster With Patient", "Cluster away from Patient")
  fit <- survfit(Surv(OS_days, OS_status) ~ group, data = survData)
  ggsurvplot(fit,  data = survData, 
             risk.table = TRUE,
             pval = TRUE,   
             conf.int = F, 
             xlim = c(0,1000),  
             break.time.by = 100, 
             ggtheme = theme_minimal(),
             risk.table.y.text.col = T, 
             risk.table.y.text = FALSE)  
  
}
