# Author: Komal S. Rathi
# Function: script to pull rsem files from all existing patients and create a matrix of expression
# script to pull clinical files from all existing patients and create metadata
# this is to be run everytime a new patient comes in - before generating the report

setwd('~/Projects/OMPARE/')
library(dplyr)
library(stringr)
library(tidyverse)
library(googlesheets4)

merge.res <- function(nm){
  sample_name <- gsub(".*PNOC","PNOC",nm)
  sample_name <- gsub('/.*','',sample_name)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# list of all PNOC patients
pat.expDat <- list.files(path = 'data/', pattern = "*.genes.results*", recursive = TRUE, full.names = T)
pat.expDat <- pat.expDat[grep('PNOC008-',  pat.expDat)]
pat.expr.mat <- lapply(pat.expDat, FUN = function(x) merge.res(x))
pat.expr.mat <- data.table::rbindlist(pat.expr.mat)

# separate gene_id and gene_symbol
pat.expr.mat <- pat.expr.mat %>% 
  mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
  separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
  unique()

# uniquify gene_symbol
pat.expr.mat <- pat.expr.mat %>% 
  group_by(sample_name) %>%
  arrange(desc(FPKM)) %>% 
  distinct(gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(gene_symbol, sample_name, FPKM) %>%
  spread(sample_name, FPKM) %>%
  column_to_rownames('gene_symbol')
colnames(pat.expr.mat) <- gsub('-[0]+','-',colnames(pat.expr.mat))

# now merge all clinical data for all patients
# read from google sheets (would need authentication the first time)
pat.clinData <- read_sheet(ss = "https://docs.google.com/spreadsheets/d/1cZgdMhIi53eNZWZMCYqRS_ThqNLb9w2TsRh2Rdg7aHY/edit#gid=0")
pat.clinData <- pat.clinData[-1,]
colnames(pat.clinData)[3] <- "KF_ParticipantID"
pat.clinData <- pat.clinData %>%
  mutate(study_id = "PNOC008") %>%
  dplyr::select(subjectID, KF_ParticipantID, tumorType, tumorLocation, ethnicity, sex, AgeAtCollection, study_id) %>%
  as.data.frame()
rownames(pat.clinData) <- pat.clinData$subjectID

common.pat <- intersect(rownames(pat.clinData), colnames(pat.expr.mat))
pat.clinData <- pat.clinData[common.pat,]
pat.expr.mat <- pat.expr.mat[,common.pat]

system('mkdir -p data/Reference/PNOC008/')
saveRDS(pat.expr.mat, file = 'data/Reference/PNOC008/PNOC008_matrix.RDS')
saveRDS(pat.clinData, file = "data/Reference/PNOC008/PNOC008_clinData.RDS")
