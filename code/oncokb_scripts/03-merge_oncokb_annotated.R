# read all annotated output files and merge in one
library(tidyverse)
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

patient <- 'PNOC008-28'
output_dir <- file.path(results_dir, patient, 'output')
annotated_files <- list.files(path = output_dir, pattern = "oncokb_.*annotated.txt", full.names = T)
maf_files <- grep('fusion|cnv', annotated_files, invert = T, value = T)
cnv_file <- grep('cnv', annotated_files, value = T)
fusion_file <- grep('fusion', annotated_files, value = T)

# mutation
for(i in 1:length(maf_files)){
  print(maf_files[i])
  dat <- read.delim(maf_files[i])
  dat <- dat %>%
    mutate(GENE = Hugo_Symbol,
           ALTERATION = gsub('^p.', '', HGVSp_Short)) %>%
    dplyr::select(GENE, ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL)
  if(i == 1){
    mutation <- dat
  }  else {
    mutation <- rbind(mutation, dat)
  }
}

# cnv
cnv <- read.delim(cnv_file)
cnv <- cnv %>%
  mutate(GENE = HUGO_SYMBOL) %>%
  dplyr::select(GENE,	ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL)

# fusion
fusion <- read.delim(fusion_file)
fusion <- fusion %>%
  mutate(GENE = strsplit(as.character(Fusion), "-")) %>% 
  unnest(GENE) %>%
  mutate(ALTERATION = 'Fusion') %>%
  dplyr::select(GENE,	ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL)

# rbind all and write out
total <- rbind(mutation, cnv, fusion)
total <- unique(total)
total <- total %>%
  filter(GENE_IN_ONCOKB == 'True')
write.table(total, file = file.path(output_dir, 'oncokb_merged_annotated.txt'), quote = F, sep = "\t", row.names = F)

# remove all oncokb files from output directory
annotated_files <- c(annotated_files, cnv_file, fusion_file)
lapply(annotated_files, FUN = function(x) system(paste0('rm ', x)))
