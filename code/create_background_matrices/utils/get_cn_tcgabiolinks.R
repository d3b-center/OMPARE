suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(tidyverse)
  library(biomaRt)
})

# get copy number from GDC
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Copy Number Variation",
                  data.type = "Gene Level Copy Number")
GDCdownload(query)

# get ploidy info from clinical data
cnv_data <- GDCprepare(query, summarizedExperiment = T)
clin_data <- colData(cnv_data) %>%
  as.data.frame() %>%
  mutate(sample_id = paste0(sample_submitter_id, "-", sample_type_id),
         barcode = gsub(",.*", "", barcode)) %>%
  dplyr::select(barcode, sample_id, paper_ABSOLUTE.ploidy) %>%
  filter(!is.na(paper_ABSOLUTE.ploidy))

# get absolute copy number scores (ASCAT)
cnv_data <- GDCprepare(query, summarizedExperiment = F)
# there are values with min_copy_number and max_copy_number but I am not using them
cnv_data <- cnv_data %>%
  dplyr::select(-c(grep('_min_copy_number|_max_copy_number', colnames(cnv_data))))
# format column names to only keep barcode
colnames(cnv_data) <- gsub(",.*", "", colnames(cnv_data))

# convert to long format
cnv_data <- cnv_data %>%
  dplyr::select(-c(chromosome, start, end)) %>%
  dplyr::rename("hgnc_symbol" = "gene_name",
                "ensembl" = "gene_id") %>%
  gather('Kids_First_Biospecimen_ID', 'copy_number', -c("hgnc_symbol", "ensembl")) %>%
  filter(!is.na(copy_number))
cnv_data <- cnv_data %>%
  inner_join(clin_data, by = c("Kids_First_Biospecimen_ID" = "barcode"))

# map status using absolute copy number and ploidy
# taken from https://cancer.sanger.ac.uk/cosmic/help/cnv/overview
cnv_data <- cnv_data %>%
  dplyr::rename("ploidy" = "paper_ABSOLUTE.ploidy") %>%
  mutate(status = ifelse(test = (ploidy <= 2.7 & copy_number >= 5) | (ploidy > 2.7 & copy_number >= 9), 
                         yes = "Gain", 
                         no = ifelse(test = (ploidy <= 2.7 & copy_number == 0) | (ploidy > 2.7 & copy_number < (ploidy-2.7)), 
                                     yes = "Loss", 
                                     no = "Neutral")))

# filter cnv to oncogenes/tsgs with gain and loss respectively
# doing this before getting cytoband info from biomart reduces search space
cancer_genes <- readRDS('data/cancer_gene_list.rds')
source('code/create_background_matrices/utils/filter_cnv.R')
cnv_data <- filter_cnv(myCNVData = cnv_data, myCancerGenes = cancer_genes)

# map cytoband information to gene symbols
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
my_regions <- getBM(c("hgnc_symbol", "band"),
                    filters = c("hgnc_symbol"),
                    values = list(hgnc_symbol = unique(cnv_data$hgnc_symbol)),
                    mart = ensembl)
my_regions <- my_regions %>%
  filter(band != "") %>%
  dplyr::rename("cytoband" = "band") 
cnv_data <- cnv_data %>%
  left_join(my_regions, by = "hgnc_symbol")

# map histology from open targets
hist_file <- read_tsv('data/OpenPedCan-analysis/data/histologies.tsv') %>%
  filter(cancer_group == "Glioblastoma Multiforme") %>%
  mutate(sample_id = gsub('[R]-[0-9A-Z]{4}-[0-9]{2}', '', Kids_First_Biospecimen_ID)) %>%
  dplyr::select(sample_id, cohort_participant_id, cohort)

# final file
cnv_data <- cnv_data %>%
  inner_join(hist_file, by = "sample_id") %>%
  dplyr::select(Kids_First_Biospecimen_ID, status, copy_number, ploidy, hgnc_symbol, cytoband, cohort_participant_id, cohort, sample_id) %>%
  unique()

  
