# TMB profile
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
module_dir <- file.path(root_dir, "code", "tmb_analysis")
patient_dir <- file.path(root_dir, "results", patient)
output_dir <- file.path(patient_dir, "output", "tmb_analysis")
dir.create(output_dir, recursive = T, showWarnings = F)

# source functions
source(file.path(module_dir, "utils", 'tmb_profile.R')) 

# output file
fname <- file.path(output_dir, "tmb_profile_output.rds")

# reference data: TMB from PBTA and TCGA
pbta_tmb <- data.table::fread(file.path(data_dir, 'tmb', 'PBTA-TMBscores_withdiseastype.txt')) %>% 
  dplyr::select(Diseasetype, Samplename, TMBscore)
tcga_not_in_pbta_tmb <- data.table::fread(file.path(data_dir, 'tmb', 'TCGA_not_in_pbta_diseasetypes_and_samples_TMBscores.txt')) %>%
  dplyr::select(Diseasetype, Samplename, TMBscore)
tcga_in_pbta_tmb <- data.table::fread(file.path(data_dir, 'tmb', 'TCGA_in_pbta_diseasetypes_and_samples_TMBscores.txt')) %>%
  dplyr::select(Diseasetype, Samplename, TMBscore)
tmb_bed_file <- data.table::fread(file.path(data_dir, 'tmb', 'ashion_confidential_exome_v2_2nt_pad.Gh38.bed'))
colnames(tmb_bed_file)  <- c('chr', 'start', 'end')
tmb_bed_file$chr <- paste0("chr", tmb_bed_file$chr)

# tmb profile
tmb_profile_output <- tmb_profile(patient_dir = patient_dir,
                                  pbta_tmb = pbta_tmb, 
                                  tcga_in_pbta_tmb = tcga_in_pbta_tmb, 
                                  tcga_not_in_pbta_tmb = tcga_not_in_pbta_tmb,
                                  tmb_bed_file = tmb_bed_file)

# save output
saveRDS(tmb_profile_output, file = fname)


