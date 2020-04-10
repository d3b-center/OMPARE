getKMPlot <- function(allCor, survData, numNeighbors = 20) {
  
  # first get all data with survival
  allCor <- allCor[intersect(rownames(allCor),survData$sample_barcode),]
  
  # subset to top 20 most correlated samples that have survival data
  colnames(allCor)[1] <- 'Correlation'
  allCor <- allCor[order(allCor$Correlation, decreasing = T),]
  mySamps <- allCor[1:numNeighbors,"sample_barcode"]
  survData$group <- survData$sample_barcode %in% mySamps
  survData$group <- ifelse(survData$group, "Cluster With Patient", "Cluster away from Patient")
  survData <- survData[survData$sample_barcode %in% rownames(allCor),]
  fit <- survfit(Surv(OS_days, OS_status) ~ group, data = survData)
  ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = survData,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimaes of survival curves.
    xlim = c(0,1000),        # present narrower X axis, but not affect
    # survival estimates.
    break.time.by = 100,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )  
  
}
