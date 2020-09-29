# Prepare PBTA data for Mutational analysis
library(tidyverse)
library(data.table)

# PBTA histology file
# only use ids where both RNA-seq and WGS are present
pbta.clin <- read.delim('data/Reference/PBTA/pbta-histologies.tsv')
sids <- pbta.clin %>%
  group_by(sample_id, experimental_strategy) %>%
  summarise(count = n()) %>%
  group_by(sample_id) %>%
  mutate(sum = sum(count)) %>%
  filter(sum == 2  & count  == 1)
pbta.clin <- pbta.clin %>%
  filter(sample_id %in% sids$sample_id)
pbta.clin <- pbta.clin %>%
  mutate(SampleID = sample_id) %>%
  dplyr::select(SampleID, Kids_First_Biospecimen_ID)

# subset to cancer genes 
cancerGenes <- read.delim("data/Reference/CancerGeneList.tsv", stringsAsFactors = F)
cancerGenes <- cancerGenes %>%
  filter(Gene_Symbol != "") %>%
  dplyr::select(-Count) %>%
  gather(key = "file", value = "type", -Gene_Symbol) %>%
  mutate(type = file)
geneListRef <- read.delim("data/Reference/genelistreference.txt", stringsAsFactors = F)
geneListRef <- subset(geneListRef, type == "TumorSuppressorGene" | type == "CosmicCensus" | type == "Oncogene")
cancerGenes <- rbind(cancerGenes, geneListRef)
gene.list <- unique(cancerGenes$Gene_Symbol)

# function to map coordinates to gene symbol
source('code/utils/createCopyNumber.R')

# chr coordinates to gene symbol map
chrMap <- read.delim("data/Reference/mart_export_genechr_mapping.txt", stringsAsFactors =F)

# read cnvs
pbta.cnv <- data.table::fread('data/Reference/PBTA/pbta-cnv-controlfreec.tsv.gz')

# function to read pbta cnv, map to gene symbol and merge
merge.cnv <- function(cnvData, genelist){
  sample_name <- unique(cnvData$Kids_First_Biospecimen_ID)
  cnvData$Kids_First_Biospecimen_ID <- NULL
  ploidy <- unique(cnvData$tumor_ploidy)
  cnvData <- cnvData %>% 
    dplyr::select(chr, start, end, copy.number, status, WilcoxonRankSumTestPvalue) %>%
    filter(WilcoxonRankSumTestPvalue < 0.05) %>%
    as.data.frame()
  cnvOut <- createCopyNumber(cnvData = cnvData, ploidy = ploidy) # map coordinates to gene symbol
  cnvOut <- cnvOut %>%
    filter(Gene %in% genelist,
           Status != "neutral") %>%
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
         Alteration_Type = stringr::str_to_title(Status),
         Alteration = paste0('Copy Number Value:', CNA),
         Study = "PBTA") %>%
  inner_join(pbta.clin, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pbta.cnv, file = 'data/Reference/PBTA/pbta-cnv-controlfreec-filtered.rds')

# read fusions (oncogenic arriba fusions)
# no need to filter fusions by pbta.clin as all are RNA-seq
pbta.fusions <- read.delim('data/Reference/PBTA/pbta-fusion-putative-oncogenic.tsv')
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
saveRDS(pbta.fusions, file = 'data/Reference/PBTA/pbta-fusion-putative-oncogenic-filtered.rds')

# read mutations (consensus)
pbta.mutations <- fread('data/Reference/PBTA/pbta-snv-consensus-mutation.maf.tsv.gz')
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
saveRDS(pbta.mutations, file = 'data/Reference/PBTA/pbta-snv-consensus-mutation-filtered.rds')

