# Prepare PBTA data for Mutational analysis
# This needs to be run with every updated version of OpenPBTA (current v18)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
pbta.dir <- file.path(ref_dir, 'pbta')

# PBTA adapt histology file
# remove unwanted samples
pbta_hist <- read.delim(file.path(pbta.dir, 'pbta-histologies-base-adapt.tsv'))
pbta_hist <- pbta_hist %>%
  filter(kf_visibility == "True")

# only use ids where both RNA-seq and WGS are present
sids <- pbta_hist %>%
  group_by(sample_id, experimental_strategy) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(sum = sum(count)) %>%
  filter(sum == 2  & count  == 1)
pbta_hist_cnv_mut <- pbta_hist %>%
  filter(sample_id %in% sids$sample_id)
pbta_hist_cnv_mut <- pbta_hist_cnv_mut %>%
  mutate(SampleID = sample_id) %>%
  dplyr::select(SampleID, Kids_First_Biospecimen_ID)

# subset to cancer genes 
cancerGenes <- readRDS(file.path(ref_dir, "cancer_gene_list.rds"))
gene.list <- unique(cancerGenes$Gene_Symbol)

# function to map coordinates to gene symbol
source(file.path(patient_level_analyses_utils, 'create_copy_number.R'))

# chr coordinates to gene symbol map
chrMap <- read.delim(file.path(ref_dir, "mart_export_genechr_mapping.txt"), stringsAsFactors =F)

# read cnvs
pbta.cnv <- data.table::fread(file.path(pbta.dir, 'pbta-cnv-controlfreec.tsv.gz'))
pbta.cnv <- pbta.cnv %>%
   filter(Kids_First_Biospecimen_ID %in% pbta_hist_cnv_mut$Kids_First_Biospecimen_ID)

# CNV analysis (chr coordinates to gene symbol map)
chr_map <- read.delim(file.path(ref_dir, 'mart_export_genechr_mapping.txt'), stringsAsFactors = F, check.names = F)
colnames(chr_map) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")

# function to read pbta cnv, map to gene symbol and merge
merge.cnv <- function(cnvData, genelist){
  sample_name <- unique(cnvData$Kids_First_Biospecimen_ID)
  cnvData$Kids_First_Biospecimen_ID <- NULL
  ploidy <- unique(cnvData$tumor_ploidy)
  cnvData <- cnvData %>% 
    dplyr::select(chr, start, end, copy.number, status, WilcoxonRankSumTestPvalue) %>%
    filter(WilcoxonRankSumTestPvalue < 0.05) %>%
    as.data.frame()
  cnvOut <- create_copy_number(cnvData = cnvData, ploidy = ploidy) # map coordinates to gene symbol
  print(head(cnvOut))
  cnvOut <- cnvOut %>%
    filter(hgnc_symbol %in% genelist,
           status != "neutral") %>%
    mutate(sample_name = sample_name) 
  if(nrow(cnvOut) > 0){
    return(cnvOut)
  }
}

pbta.cnv <- plyr::ddply(.data = pbta.cnv, 
                        .variables = 'Kids_First_Biospecimen_ID', 
                        .fun = function(x) merge.cnv(cnvData = x, genelist = gene.list))
pbta.cnv <- pbta.cnv %>%
  mutate(Alteration_Datatype = "CNV",
         Alteration_Type = stringr::str_to_title(status),
         Alteration = paste0('Copy Number Value:', copy.number),
         Study = "PBTA",
         Gene = hgnc_symbol) %>%
  inner_join(pbta_hist_cnv_mut, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pbta.cnv, file = file.path(pbta.dir, 'pbta-cnv-controlfreec-filtered.rds'))

# read fusions (oncogenic arriba fusions)
# no need to filter fusions by pbta.clin as all are RNA-seq
pbta.fusions <- read.delim(file.path(pbta.dir, 'pbta-fusion-putative-oncogenic.tsv'))
pbta.fusions <- pbta.fusions[grep('ARRIBA', pbta.fusions$CalledBy),]
pbta.fusions <- pbta.fusions  %>%
  unite(col = "Gene", Gene1A, Gene1B, Gene2A, Gene2B, sep = ", ", na.rm = T)  %>%
  mutate(Alteration_Datatype = "Fusion",
         Alteration_Type = stringr::str_to_title(Fusion_Type),
         Alteration = FusionName,
         Kids_First_Biospecimen_ID = Sample,
         Study = "PBTA") %>%
  inner_join(pbta.clin, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  separate_rows(Gene, convert = TRUE) %>%
  filter(Gene %in% gene.list) %>%
  unique()
saveRDS(pbta.fusions, file = file.path(pbta.dir, 'pbta-fusion-putative-oncogenic-filtered.rds'))

# read mutations (consensus)
pbta.mutations <- fread(file.path(pbta.dir, 'pbta-snv-consensus-mutation.maf.tsv.gz'))
keepVC <- c("Nonsense_Mutation", "Missense_Mutation", 
            "Splice_Region", "Splice_Site",
            "3'UTR", "5'UTR", 
            "In_Frame_Ins", "In_Frame_Del",
            "Frame_Shift_Ins", "Frame_Shift_Del")
keepVI <- c("MODIFIER", "MODERATE", "HIGH")
pbta.mutations <- pbta.mutations %>%
  filter(BIOTYPE == "protein_coding",
         Variant_Classification %in% keepVC,
         IMPACT %in% keepVI,
         Hugo_Symbol %in% gene.list) %>%
  mutate(Gene = Hugo_Symbol,
         Alteration_Datatype = "Mutation",
         Alteration_Type = Variant_Classification,
         Alteration = HGVSp_Short,
         Kids_First_Biospecimen_ID = Tumor_Sample_Barcode,
         Study = "PBTA") %>%
  inner_join(pbta.clin, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pbta.mutations, file = file.path(pbta.dir, 'pbta-snv-consensus-mutation-filtered.rds'))

