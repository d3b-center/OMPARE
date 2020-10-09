# Author: Komal S. Rathi
# Date: 02/28/2020
# Function: Script to read from google sheets and create clinical file 
# This will be called from within run_OMPARE.R

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

option_list <- list(
  make_option(c("-s", "--sheet"), type = "character",
              help = "PNOC008 Manifest file (.xlsx)"),
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
dat <- readxl::read_xlsx(sheet, sheet = 1)
colnames(dat) <- gsub('[() ]', '.', colnames(dat))
dat <- dat %>%
  filter_all(any_vars(!is.na(.))) %>%
  mutate(PNOC.Subject.ID = gsub('P-','PNOC008-', PNOC.Subject.ID)) %>%
  dplyr::filter(PNOC.Subject.ID == patient)

# create clinical file
df <- dat %>%
  mutate(KF_ParticipantID = PNOC.Subject.ID,
         subjectID = PNOC.Subject.ID,
         reportDate = Sys.Date(),
         tumorType = Diagnosis.a,
         tumorLocation = Primary.Site.a,
         ethnicity = Ethnicity,
         age_diagnosis_days = Age.at.Diagnosis..in.days.,
         age_collection_days = Age.at.Collection..in.days.,
         sex = Gender) %>%
  dplyr::select(subjectID, reportDate, tumorType,	tumorLocation, ethnicity, sex, age_diagnosis_days, age_collection_days, KF_ParticipantID, library_name)

# write out
fname <- file.path(dir, "clinical", "patient_report.txt")
write.table(df, file = fname, sep = "\t", quote = F, row.names = F)
