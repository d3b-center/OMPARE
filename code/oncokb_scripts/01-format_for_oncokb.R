# filter and format cnv and fusion files for oncokb annotator
library(tidyverse)
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

patient <- 'PNOC008-28'
output_dir <- file.path(results_dir, patient, 'output')

# CNV - cnaAnnotator.py - expects log2 ratio data as input (leverage cnvkit's *gainloss.txt, which already is a table of gene names and log2 ratio). Include log2 ratios of 0.2<=x<0.7 for single copy gain, x>=0.7 for amplification, -1.1<x<=-0.25 for single copy loss, and x<-1.1 for homozygous loss. Also requires reformatting cnv data per the oncokb specification
cnv_dir <- file.path(results_dir, patient, 'copy-number-variations')
cnv_data <- list.files(path = cnv_dir, pattern = "*gainloss*", full.names = T)
cnv_data <- read.delim(cnv_data)
cnv_data <- cnv_data %>%
  dplyr::select(-c(chromosome, start, end, depth, weight, n_bins))

# write out matrix of GISTIC like output
cnv_data$value <- 0
cnv_data$value[cnv_data$log2 >= 0.7] <- 2
cnv_data$value[cnv_data$log2 >= 0.2 & cnv_data$log2 < 0.7] <- 1
cnv_data$value[cnv_data$log2 > -1.1 & cnv_data$log2 <= -0.25] <- -1
cnv_data$value[cnv_data$log2 < -1.1] <- -2
cnv_data <- cnv_data %>%
  mutate(`Gene Symbol` = gene,
         `Locus ID` = 0,
         Cytoband = 0)  %>%
  dplyr::select(c(`Gene Symbol`, `Locus ID`, Cytoband, value)) %>%
  rename(!!patient := value)
write.table(cnv_data, file = file.path(output_dir, 'oncokb_cnv.txt'), quote = F, sep = "\t", row.names = F)

# Fusion - take arriba output, filter for high confidence fusions only, and reformat fusion data to format specified by oncokb (e.g. to id fields + hyphenated fusion name). 
# Take STAR fusion output, reformat fusion names (e.g. replace double hyphen with single hyphen), take those only with FFPM>0.1 AND LargeAnchorSupport (YES_LDAS), reformat to oncokb specifications linked above, and merge with arriba data. 
# Run FusionAnnotator.py on non-redundant arriba + star-fusion calls.
fusion_dir <- file.path(results_dir, patient, 'gene-fusions')

# arriba
arriba_data <- list.files(path = fusion_dir, pattern = "*arriba*", full.names = T)
arriba_data <- read.delim(arriba_data)
arriba_data <- arriba_data %>%
  filter(confidence == "high") %>%
  mutate(Fusion = paste0(X.gene1, '-', gene2),
         Tumor_Sample_Barcode = patient) %>%
  dplyr::select(Tumor_Sample_Barcode, Fusion)

# star fusion
star_data <- list.files(path = fusion_dir, pattern = "*STAR*", full.names = T)
star_data <- read.delim(star_data)
star_data <- star_data %>%
  mutate(Fusion = gsub('--','-',X.FusionName),
         Tumor_Sample_Barcode = patient) %>%
  filter(FFPM > 0.1,
         LargeAnchorSupport == "YES_LDAS") %>%
  dplyr::select(Tumor_Sample_Barcode, Fusion)

# combine both and write out
fusion_data <- unique(rbind(star_data, arriba_data))
write.table(fusion_data, file = file.path(output_dir, 'oncokb_fusion.txt'), quote = F, sep = "\t", row.names = F)
  
  