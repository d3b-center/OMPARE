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

# read cancer genes
cancer_genes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))
cancer_genes <- unique(cancer_genes$Gene_Symbol)

# mutation
for(i in 1:length(maf_files)){
  print(maf_files[i])
  dat <- read.delim(maf_files[i])
  dat <- dat %>%
    mutate(GENE = Hugo_Symbol,
           ALTERATION = gsub('^p.', '', HGVSp_Short), # add HGVSp short
           ALTERATION = ifelse(ALTERATION == "", NA, ALTERATION), # replace "" with NAs
           ALTERATION_v2 = ifelse(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', # add Oncogenic mutations
                                                             'Frame_Shift_Del', 'Frame_Shift_Ins', 
                                                             'In_Frame_Del', 'In_Frame_Ins', 'CDS_position') & 
                                 IMPACT %in% c('MODERATE', 'HIGH'), 
                               yes = 'Oncogenic Mutations', 
                               no = NA)) %>%
    unite(., col = "ALTERATION",  ALTERATION, ALTERATION_v2, na.rm=TRUE, sep = ", ") %>% # create comma separated values in case of HGVSp and "Oncogenic Mutations"
    mutate(ALTERATION = strsplit(ALTERATION, ", ")) %>% 
    unnest(ALTERATION) %>% # split into new rows
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
  mutate(ALTERATION_v2 = ifelse(GENE %in% cancer_genes, 'Oncogenic Mutations', NA)) %>% # use cancer genes to define Oncogenic Mutations
  unite(., col = "ALTERATION",  ALTERATION, ALTERATION_v2, na.rm=TRUE, sep = ", ") %>% # create comma separated values of "Amplification/Deletion" and "Oncogenic Mutations"
  mutate(ALTERATION = strsplit(ALTERATION, ", ")) %>% 
  unnest(ALTERATION) %>% # split into new rows
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
total[is.na(total)] <- ''
write.table(total, file = file.path(output_dir, 'oncokb_merged_annotated.txt'), quote = F, sep = "\t", row.names = F)

# remove all oncokb files from output directory
# annotated_files <- c(annotated_files, cnv_file, fusion_file)
# lapply(annotated_files, FUN = function(x) system(paste0('rm ', x)))
