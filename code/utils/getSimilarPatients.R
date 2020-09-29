##################################
# Function to get similar patients
##################################

getSimilarPatients <- function(allCor, clinData) {
  allCor <- clinData %>%
    inner_join(allCor, by = c("sample_barcode" = "nearest_neighbor")) %>%
    arrange(distance)
  return(allCor)
}
