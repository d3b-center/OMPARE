##################################
# Code to pull sample/patient info
##################################

# assign patient clinical info to global env
lapply(colnames(sampleInfo), FUN = function(x) assign(x, sampleInfo[,x], envir = globalenv()))

##################
#Header Information
##################

#Get Patient Name
getPatientName <- function() {
  return(patientName)
}

# Get Report Date
getReportDate <- function() {
  return(reportDate)
}

# Get Report Version
getReportVersion <- function() {
  return(reportVersion)
}

##################
#Patient/Sample Information
##################

# Get Subject ID
getSubjectID<- function() {
  return(subjectID)
}

# Get Sex
getSex <- function() {
  return(sex)
}

# Get DOB
getDOB <- function() {
  return(dob)
}

# Get ethnicity
getEthnicity <- function() {
  return(ethnicity)
}

##################

# Get Medical Facility
getMedicalFacility <- function() {
  return(medicalFacility)
}

# Get Primary Physician
getPrimPhysician <- function() {
  return(primPhysician)
}

# Get Pathologist
getPathologist <- function() {
  return(pathologist)
}

# Get Collection Date
getLabDirector <- function() {
  return(labDirector)
}

##################

# Get Collection Date
getCollectionDate <- function() {
  return(collectionDate)
}

# Get Tumor Location
getTumorLocation <- function() {
  return(tumorLocation)
}

# Get Tumor Type
getTumorType <- function() {
  return(tumorType)
}

# Get P/R
getPrimRelapse <- function() {
  return(primRelapse)
}





