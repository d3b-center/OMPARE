##################################
# Function to get similar patients
##################################

getSimilarPatients <- function(numNeighbors=15) {
  if(exists('expData')){
    mySamps <- allCor[1:numNeighbors,"samps"]
    clinDataTmp <- clinData.full[clinData.full$Kids_First_Biospecimen_ID %in% mySamps,]
  } else {
    clinDataTmp <- data.frame()
  }
  return(clinDataTmp)
}