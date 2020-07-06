##################################
# Code to pull sample/patient info
##################################

# assign patient clinical info to global env
lapply(colnames(sampleInfo), FUN = function(x) assign(x, sampleInfo[,x], envir = globalenv()))

##################
#Header Information
##################

#Get Patient Name
getPatientID <- function() {
  return(KF_ParticipantID)
}

# Get Report Date
getReportDate <- function() {
  return(reportDate)
}

# Get Report Version
# getReportVersion <- function() {
#   return(reportVersion)
# }

##################
#Patient/Sample Information
##################

# Get Subject ID
getSubjectID<- function() {
  return(subjectID)
}

# sex
getSex <- function() {
  return(sex)
}

# age at diagnosis
getAge <- function() {
  return(age_diagnosis_days)
}

# ethnicity
getEthnicity <- function() {
  return(ethnicity)
}

##################

# Get Medical Facility
# getMedicalFacility <- function() {
#   return(medicalFacility)
# }

# Get Primary Physician
# getPrimPhysician <- function() {
#   return(primPhysician)
# }
# 
# # Get Pathologist
# getPathologist <- function() {
#   return(pathologist)
# }
# 
# # Get Collection Date
# getLabDirector <- function() {
#   return(labDirector)
# }

##################

# age at collection
getCollectionDate <- function() {
  return(age_collection_days)
}

# tumor location
getTumorLocation <- function() {
  return(tumorLocation)
}

# tumor type
getTumorType <- function() {
  return(tumorType)
}

# # Get P/R
# getPrimRelapse <- function() {
#   return(primRelapse)
# }





