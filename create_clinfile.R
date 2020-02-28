# Author: Komal S. Rathi
# Date: 02/28/2020
# Function: Script to read from google sheets and create clinical file 

library(googlesheets4)
library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("-s", "--sheet"), type = "character",
              help = "Link to PNOC008 Google sheet"),
  make_option(c("-d", "--dir"), type = "character",
              help = "Path to PNOC008 patient folder (top directory)"),
  make_option(c("-p", "--patient"), type = "character",
              help = "Patient identifier for PNOC008. e.g. PNOC008-1, PNOC008-10 etc")
)

opt <- parse_args(OptionParser(option_list = option_list))
sheet <- opt$sheet
dir <- opt$dir
patient <- opt$patient

# read from google sheets (would need authentication the first time)
dat <- read_sheet(sheet)
dat <- dat %>%
  dplyr::filter(subjectID == patient)
colnames(dat)[3] <- "KF_ParticipantID"

# create clinical file
df <- data.frame(matrix(ncol = 15, nrow = 1))
colnames(df) <- c("subjectID","reportDate","reportVersion","primRelapse","tumorType","tumorLocation","collectionDate","labDirector","pathologist","primPhysician","medicalFacility","ethnicity","age_years","sex","KF_ParticipantID")
df$subjectID <- dat$subjectID
df$reportDate <- Sys.Date()
df$reportVersion <- "V1.01"
df$primRelapse <- "Primary"
df$tumorType <- dat$tumorType
df$tumorLocation <- dat$tumorLocation
df$collectionDate <- "n/a"
df$labDirector  <- "n/a"
df$pathologist <- "n/a"
df$primPhysician <- "n/a"
df$medicalFacility <- "UCSF"
df$ethnicity <- dat$ethnicity
df$age_years <- dat$AgeAtCollection[[1]]
df$sex <- dat$sex
df$KF_ParticipantID <- dat$KF_ParticipantID

# write out
fname <- file.path(dir, "Clinical", "patient_report.txt")
write.table(df, file = fname, sep = "\t", quote = F, row.names = F)
