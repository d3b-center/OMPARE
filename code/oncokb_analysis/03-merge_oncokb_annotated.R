# read all annotated output files and merge in one
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# arguments
option_list <- list(
  make_option(c("-s", "--snv_caller"), type = "character",
              help = "SNV caller pattern: lancet, vardict, consensus, strelka2, mutect2 and all"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
snv_caller <- opt$snv_caller
output_dir <- opt$output_dir

# set directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# input files
maf_file <- file.path(output_dir, paste0('oncokb_', snv_caller, '_annotated.txt'))
cnv_file <- file.path(output_dir, 'oncokb_cnv_annotated.txt')
fusion_file <- file.path(output_dir, 'oncokb_fusion_annotated.txt')

# read cancer genes and split into oncogenes and tsgs
cancer_genes <- readRDS(file.path(data_dir, 'cancer_gene_list.rds'))
oncogenes <- cancer_genes %>%
  filter(type %in% c("Is.Oncogene", "Oncogene")) %>%
  .$Gene_Symbol %>% unique
tsgs <- cancer_genes %>%
  filter(type %in% c("Is.Tumor.Suppressor.Gene", "TumorSuppressorGene")) %>%
  .$Gene_Symbol %>% unique

# mutation
mutation <- data.table::fread(maf_file)
mutation <- mutation %>%
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
  unite(., col = "CITATIONS",  DX_CITATIONS, PX_CITATIONS, TX_CITATIONS, na.rm=TRUE, sep = ", ") %>% # create comma separated values in case of HGVSp and "Oncogenic Mutations"
  mutate(CITATIONS = strsplit(CITATIONS, ", ")) %>% 
  unnest(CITATIONS) %>% # split into new rows
  dplyr::select(GENE, ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL) %>%
  unique()

# cnv
cnv <- read.delim(cnv_file)
cnv <- cnv %>%
  mutate(GENE = HUGO_SYMBOL) %>%
  mutate(ALTERATION_v2 = ifelse((GENE %in% oncogenes & ALTERATION == "Amplification") | 
                                  (GENE %in% tsgs & ALTERATION == "Deletion"), yes = 'Oncogenic Mutations', no = NA)) %>% # use cancer genes to define Oncogenic Mutations
  unite(., col = "ALTERATION",  ALTERATION, ALTERATION_v2, na.rm=TRUE, sep = ", ") %>% # create comma separated values of "Amplification/Deletion" and "Oncogenic Mutations"
  mutate(ALTERATION = strsplit(ALTERATION, ", ")) %>% 
  unnest(ALTERATION) %>% # split into new rows
  unite(., col = "CITATIONS",  DX_CITATIONS, PX_CITATIONS, TX_CITATIONS, na.rm=TRUE, sep = ", ") %>% # create comma separated values in case of HGVSp and "Oncogenic Mutations"
  mutate(CITATIONS = strsplit(CITATIONS, ", ")) %>% 
  unnest(CITATIONS) %>% # split into new rows
  dplyr::select(GENE,	ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL) %>%
  unique()

# fusion
fusion <- read.delim(fusion_file)
fusion <- fusion %>%
  mutate(GENE = strsplit(as.character(Fusion), "-")) %>% 
  unnest(GENE) %>%
  mutate(ALTERATION = 'Fusion') %>%
  unite(., col = "CITATIONS",  DX_CITATIONS, PX_CITATIONS, TX_CITATIONS, na.rm=TRUE, sep = ", ") %>% # create comma separated values in case of HGVSp and "Oncogenic Mutations"
  mutate(CITATIONS = strsplit(CITATIONS, ", ")) %>% 
  unnest(CITATIONS) %>% # split into new rows
  dplyr::select(GENE,	ALTERATION, GENE_IN_ONCOKB, VARIANT_IN_ONCOKB, ONCOGENIC, MUTATION_EFFECT, HIGHEST_LEVEL, CITATIONS, HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL) %>%
  unique()

# rbind all and write out
output_file <- file.path(output_dir, paste0('oncokb_merged_', snv_caller, '_annotated.txt'))
total <- rbind(mutation, cnv, fusion)
total <- unique(total)
total <- total %>%
  mutate(GENE_IN_ONCOKB = toupper(GENE_IN_ONCOKB),
         VARIANT_IN_ONCOKB = toupper(VARIANT_IN_ONCOKB)) %>%
  filter(GENE_IN_ONCOKB == 'TRUE') %>%
  mutate_all(as.character)
total[is.na(total)] <- ""
write.table(total, file = output_file, quote = F, sep = "\t", row.names = F)

