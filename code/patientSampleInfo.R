###############################
#Purpose: Code to pull sample / patient information
#Date: 2/11/2019
#Author: Pichai Raman
###############################

#All hardcoded, in future will read from somewhere and set

#Set variables
patientName <- " n/a ";
reportDate <- "5/07/2019"
reportVersion <- "V1.01"
primRelapse <- "Primary"
tumorType <- "High-grade glioma/astrocytoma (WHO grade III/IV)"
tumorLocation <- "Brain"
collectionDate <- "3/08/2019"
labDirector <- "n/a"
pathologist <- "n/a"
primPhysician <- "n/a"
medicalFacility <- "UCSF"
ethnicity <- "Unknown"
dob <- "n/a"
sex <- "Male"
subjectID <- "PNOC008-2"

##################
#Header Information
##################

#Get Patient Name
getPatientName <- function()
{
  return(patientName)
}

#Get Report Date
getReportDate <- function()
{
  return(reportDate)
}

#Get Report Version
getReportVersion <- function()
{
  return(reportVersion)
}

##################
#Patient/Sample Information
##################

#Get Subject ID
getSubjectID<- function()
{
  return(subjectID)
}

#Get Sex
getSex <- function()
{
  return(sex)
}

#Get DOB
getDOB <- function()
{
  return(dob)
}

#Get ethnicity
getEthnicity <- function()
{
  return(ethnicity)
}

##################

#Get Medical Facility
getMedicalFacility <- function()
{
  return(medicalFacility)
}

#Get Primary Physician
getPrimPhysician <- function()
{
  return(primPhysician)
}

#Get Pathologist
getPathologist <- function()
{
  return(pathologist)
}

#Get Collection Date
getLabDirector <- function()
{
  return(labDirector)
}

##################

#Get Collection Date
getCollectionDate <- function()
{
  return(collectionDate)
}

#Get Tumor Location
getTumorLocation <- function()
{
  return(tumorLocation)
}

#Get Tumor Type
getTumorType <- function()
{
  return(tumorType)
}

#Get P/R
getPrimRelapse <- function()
{
  return(primRelapse)
}





