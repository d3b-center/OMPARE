##################################
# Function to get similar patients
##################################

getSimilarPatients <- function(allCor, clinData, numNeighbors = 20) {
  colnames(allCor)[1] <- 'Correlation'
  
  # return top 20 most correlated if numNeighbors is not NULL
  if(!is.null(numNeighbors)){
    allCor <- allCor[1:numNeighbors,]
    allCor <- merge(clinData, allCor, by =  'sample_barcode')
    allCor <- allCor[order(allCor$Correlation, decreasing = T),]
  } else {
    # else return anything with pearson correlation > 0.5
    allCor <- allCor %>%
      filter(Correlation > 0.5) %>%
      inner_join(clinData, by = c("sample_barcode")) %>%
      arrange(desc(Correlation))
  }
  return(allCor)
}
