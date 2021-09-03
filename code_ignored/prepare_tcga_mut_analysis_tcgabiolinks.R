library(TCGAbiolinks)
library(tidyverse)
library(reshape2)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))
tcga.dir <- file.path(ref_dir, 'tcga')

# clinical file
tcga.clin <- readRDS(file.path(tcga.dir, 'tcga_gbm_clinical.rds'))

# subset to cancer genes 
cancerGenes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))
gene_list <- unique(cancerGenes$Gene_Symbol)

# mutations
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Simple Nucleotide Variation", 
                  access="open", 
                  legacy = F, 
                  data.type = "Masked Somatic Mutation", 
                  workflow.type = "MuTect2 Variant Aggregation and Masking")
GDCdownload(query)
mut_data <- GDCprepare(query, add.gistic2.mut = T)
mut_data$Tumor_Sample_Barcode <- gsub('D.*|W.*', '',mut_data$Tumor_Sample_Barcode)
mut_data <- mut_data %>%
  filter(Tumor_Sample_Barcode %in% tcga.clin$sample_barcode)

# filters
keepVC <- c("Nonsense_Mutation", "Missense_Mutation", 
            "Splice_Region", "Splice_Site",
            "3'UTR", "5'UTR", 
            "In_Frame_Ins", "In_Frame_Del",
            "Frame_Shift_Ins", "Frame_Shift_Del")
keepVI <- c("MODIFIER", "MODERATE", "HIGH")
tcga_mutations <- mut_data %>%
  filter(Variant_Classification %in% keepVC,
         Hugo_Symbol %in% gene_list) %>%
  mutate(Gene = Hugo_Symbol,
         Alteration_Datatype = "Mutation",
         Alteration_Type = Variant_Classification,
         Alteration = HGVSp_Short,
         Study = "TCGA GBM",
         SampleID = Tumor_Sample_Barcode,
         Kids_First_Biospecimen_ID = SampleID) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(tcga_mutations, file = file.path(tcga.dir, 'tcga_gbm_mutation_filtered.rds'))

# query expression data to get gene id and symbol mapping
query1 <- GDCquery(project = "TCGA-GBM",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts",
                   experimental.strategy = "RNA-Seq", legacy = F)
GDCdownload(query1)
data1 <- GDCprepare(query1)
annot <- SummarizedExperiment::rowRanges(data1)
annot <- as.data.frame(annot)
annot <- annot %>%
  dplyr::select(original_ensembl_gene_id, external_gene_name) %>%
  unique()

# copy number
query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number Scores",              
                  access="open", legacy = F)
GDCdownload(query)
cnv_data <- GDCprepare(query)

# convert matrix to long format
cnv_data <- cnv_data %>%
  dplyr::select(-c(`Gene ID`, `Cytoband`))
cnv_data <- melt(cnv_data)

# format identifiers so we can map to RNA
cnv_data$SampleID <- gsub('D.*|W.*', '', cnv_data$variable)

# only keep gain/loss
cnv_data <- cnv_data %>%
  filter(SampleID %in% tcga.clin$sample_barcode,
         value != 0) %>%
  mutate(Alteration_Type = ifelse(value == 1, 'Gain', 'Loss'),
         Alteration_Datatype = 'CNV',
         Alteration = '',
         Kids_First_Biospecimen_ID = SampleID,
         Study = 'TCGA GBM') %>%
  unique()

# filter to cancer genes
tcga_cnv_filtered <- cnv_data %>%
  inner_join(annot, by = c('Gene Symbol' = 'original_ensembl_gene_id')) %>%
  mutate(Gene = external_gene_name) %>%
  filter(Gene %in% gene_list) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(tcga_cnv_filtered, file = file.path(tcga.dir, 'tcga_gbm_cnv_filtered.rds'))

