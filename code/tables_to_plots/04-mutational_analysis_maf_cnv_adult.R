suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(TCGAbiolinks)
  library(optparse)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier, i.e. PNOC008-XX"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
pnoc008_sample_of_interest <- opt$patient
output_dir <- opt$output_dir

# Define directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
adult_dir <- file.path(data_dir, "tcga")
pnoc008_dir <- file.path(data_dir, "pnoc008")

## P6 Adult Tumor Analysis - Mutational Analysis
### Query for TCGA MAF file and filter using histology file
# Query with TCGABiolinks
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Simple Nucleotide Variation", 
                  access = "open", 
                  legacy = F, 
                  data.type = "Masked Somatic Mutation", 
                  workflow.type = "MuTect2 Variant Aggregation and Masking")

# currently the query is not working so this is a workaround
if(exists('query')){
  GDCdownload(query)
  tcga_maf <- GDCprepare(query, add.gistic2.mut = T)
} else {
  # read pre-existing
  tcga_maf <- data.table::fread("GDCdata/TCGA-GBM/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/da904cd3-79d7-4ae3-b6c0-e7127998b3e6/TCGA.GBM.mutect.da904cd3-79d7-4ae3-b6c0-e7127998b3e6.DR-10.0.somatic.maf.gz")
}
tcga_maf$Tumor_Sample_Barcode <- gsub('D.*|W.*', '',tcga_maf$Tumor_Sample_Barcode)

# Filter to the samples that we have RNA information 
tcga_clin <- readRDS(file.path(adult_dir, "tcga_gbm_clinical.rds")) %>%
  dplyr::rename(Tumor_Sample_Barcode = sample_barcode)

tcga_maf <- tcga_maf %>%
  filter(Tumor_Sample_Barcode %in% tcga_clin$Tumor_Sample_Barcode)

### Read in MAF for PNOC008
consensus_maf_pnoc008 <- readRDS(file.path(pnoc008_dir, "pnoc008_consensus_mutation.rds"))

### Read in mutational analysis
mutational_analysis <- readRDS(file.path(output_dir, "mutational_analysis_adult.rds"))

### Genearte figures for recurrent mutations 
#### Prepare combined MAF for recurrent mutations
recurrent_subjects <- mutational_analysis$recurrent_alterations %>% 
  filter(Alteration_Datatype == "Mutation") %>% 
  pull(subject_id) %>% 
  append(pnoc008_sample_of_interest) %>%
  unique()

recurrent_alterations_pnoc008 <- recurrent_subjects[recurrent_subjects %in% consensus_maf_pnoc008$sample_name] 
recurrent_alterations_tcga <- recurrent_subjects[!recurrent_subjects %in% consensus_maf_pnoc008$sample_name]

stopifnot(recurrent_alterations_tcga %in% tcga_maf$Tumor_Sample_Barcode)

# Use the kids_first_biospecimen_id to filter TCGA combined consensus SNV
matched_maf_tcga_recurrent <- tcga_maf %>% 
  filter(Tumor_Sample_Barcode %in% recurrent_alterations_tcga)

# filter based on sample_name
matched_maf_pnoc008_recurrent <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% recurrent_alterations_pnoc008)

# combine the maf to prepare for MAF object
common_columns <- intersect(colnames(tcga_maf), colnames(consensus_maf_pnoc008))

matched_maf_tcga_recurrent <- matched_maf_tcga_recurrent %>% 
  select(common_columns)
matched_maf_pnoc008_recurrent <- matched_maf_pnoc008_recurrent %>% 
  select(common_columns)

combined_maf_recurrent <- matched_maf_tcga_recurrent %>%
  rbind(matched_maf_pnoc008_recurrent) %>% 
  dplyr::rename(AAChange = HGVSp_Short) 

#### Prepare combined histology for recurrent mutations
match_pnoc008 <- consensus_maf_pnoc008 %>% 
  select(sample_name, Tumor_Sample_Barcode) %>% 
  distinct() 

histology_pnoc008 <- readRDS(file.path(pnoc008_dir, "pnoc008_clinical.rds")) %>%
  left_join(match_pnoc008, by = c("subjectID" =  "sample_name")) %>%
  dplyr::rename(
    gender = sex,
    age_at_diagnosis_in_days = age_diagnosis_days
  ) 

histology_pnoc008_recurrent <- histology_pnoc008 %>% 
  filter(subjectID %in% recurrent_alterations_pnoc008) 

tcga_clin_recurrent <- tcga_clin %>% 
  filter(Tumor_Sample_Barcode %in% recurrent_alterations_tcga)

common_columns_histology <- intersect(colnames(histology_pnoc008), colnames(tcga_clin))

histology_pnoc008_recurrent <- histology_pnoc008_recurrent %>%
  select(common_columns_histology)
tcga_clin_recurrent <- tcga_clin_recurrent%>% 
  select(common_columns_histology)

combined_histology_recurrent <- rbind(histology_pnoc008_recurrent, tcga_clin_recurrent)

#### Clean up the MAF file and generate a MAF object

# only certain variant classifications are allowed
# https://github.com/PoisonAlien/maftools/blob/master/R/lollipopPlot.R#L130-L131
combined_maf_recurrent <- combined_maf_recurrent %>%
  filter(Variant_Classification %in% c("Nonstop_Mutation", "Missense_Mutation",
                                       "Nonsense_Mutation", "Splice_Site",
                                       "Frame_Shift_Ins", "Frame_Shift_Del",
                                       "In_Frame_Del", "In_Frame_Ins"))

maf_object_recurrent <- read.maf(maf = combined_maf_recurrent , clinicalData = combined_histology_recurrent)

#### Plot MAF object
pdf(file.path(output_dir, "mutational_recurrent_adult.pdf"))
plotmafSummary(maf = maf_object_recurrent, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#### Generate lollipop plot 
# First find the ten most mutated gene
summary_genes <- getGeneSummary(maf_object_recurrent) %>% 
  arrange(desc(AlteredSamples)) 

# find the most mutated 10 genes - PAML2-AKAP2 does not have structure available 
top10_gene_list <- summary_genes %>%
  head(10) %>% 
  pull(Hugo_Symbol) %>% 
  unique()

# check genes in db
domain_data_src <- system.file('extdata', 'prot_len.txt.gz', package = 'maftools')
domain_data_src <- data.table::fread(domain_data_src)
top10_gene_list <- top10_gene_list[top10_gene_list %in% domain_data_src$Hugo_Symbol]

top10_gene_list_position <- lapply(top10_gene_list, function(x) {
  protein_changes <- combined_maf_recurrent %>% 
    filter(Hugo_Symbol == x) %>% 
    filter(AAChange != ".") %>% 
    pull(AAChange)
  positions <- gsub(".*?([0-9]+).*", "\\1", protein_changes)
})

names(top10_gene_list_position) <- top10_gene_list 

pdf(file.path(output_dir, "lollipop_recurrent_adult.pdf"), onefile = TRUE)
lapply(top10_gene_list, function(x){
  lollipopPlot(
  maf = maf_object_recurrent,
  gene = x,
  showMutationRate = TRUE, 
  labelPos = top10_gene_list_position[[x]]
)
})
dev.off()

### Genearte figures for shared genes
#### Prepare combined MAF for shared mutations
shared_subjects <- mutational_analysis$shared_genes %>% 
  filter(Alteration_Datatype == "Mutation") %>% 
  pull(subject_id) %>% 
  append(pnoc008_sample_of_interest) %>%
  unique()

shared_subjects_pnoc008 <- shared_subjects[shared_subjects %in% consensus_maf_pnoc008$sample_name]
shared_subjects_tcga <- shared_subjects[!shared_subjects %in% consensus_maf_pnoc008$sample_name]

# stopifnot(shared_subjects_tcga%in% tcga_maf$Tumor_Sample_Barcode)

# Use the kids_first_biospecimen_id to filter TCGA combined consensus SNV
matched_maf_tcga_shared <- tcga_maf %>% 
  filter(Tumor_Sample_Barcode %in% shared_subjects_tcga)

# filter based on sample_name
matched_maf_pnoc008_shared <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% shared_subjects_pnoc008)

# combine the maf to prepare for MAF object
matched_maf_tcga_shared <- matched_maf_tcga_shared %>% 
  select(common_columns)
matched_maf_pnoc008_shared <- matched_maf_pnoc008_shared %>% 
  select(common_columns)

combined_maf_shared <- matched_maf_tcga_shared %>%
  rbind(matched_maf_pnoc008_shared) %>%
  dplyr::rename(AAChange = HGVSp_Short) 

#### Prepare combined histology for recurrent mutations
histology_pnoc008_shared <- histology_pnoc008 %>% 
  filter(subjectID %in% shared_subjects_pnoc008) 

tcga_clin_shared <- tcga_clin %>% 
  filter(Tumor_Sample_Barcode %in% shared_subjects_tcga)

histology_pnoc008_shared <- histology_pnoc008_shared %>% select(common_columns_histology)
tcga_clin_shared <- tcga_clin_shared %>% 
  select(common_columns_histology)

combined_histology_shared <- rbind(histology_pnoc008_shared, tcga_clin_shared)

#### Clean up the MAF file and generate a MAF object
# only certain variant classifications are allowed
# https://github.com/PoisonAlien/maftools/blob/master/R/lollipopPlot.R#L130-L131
combined_maf_shared <- combined_maf_shared %>%
  filter(Variant_Classification %in% c("Nonstop_Mutation", "Missense_Mutation",
                                       "Nonsense_Mutation", "Splice_Site",
                                       "Frame_Shift_Ins", "Frame_Shift_Del",
                                       "In_Frame_Del", "In_Frame_Ins"))

maf_object_shared <- read.maf(maf = combined_maf_shared, 
                              clinicalData = combined_histology_shared)

#### Plot MAF object
pdf(file.path(output_dir, "mutational_shared_adult.pdf"))
plotmafSummary(maf = maf_object_shared, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#### Generate lollipop plot 
# First find the ten most mutated gene
summary_genes <- getGeneSummary(maf_object_shared) %>% 
  arrange(desc(AlteredSamples)) 

# find the most mutated 10 genes 
top10_gene_list <- summary_genes %>%
  head(10) %>% 
  pull(Hugo_Symbol) %>% 
  unique()

# check genes in db
domain_data_src <- system.file('extdata', 'prot_len.txt.gz', package = 'maftools')
domain_data_src <- data.table::fread(domain_data_src)
top10_gene_list <- top10_gene_list[top10_gene_list %in% domain_data_src$Hugo_Symbol]

top10_gene_list_position <- lapply(top10_gene_list, function(x) {
  protein_changes <- combined_maf_shared %>% filter(Hugo_Symbol == x) %>% 
    filter(AAChange != ".") %>% 
    pull(AAChange)
  positions <- gsub(".*?([0-9]+).*", "\\1", protein_changes)
})

names(top10_gene_list_position) <- top10_gene_list 

pdf(file.path(output_dir, "lollipop_shared_adult.pdf"), onefile = TRUE)
lapply(top10_gene_list, function(x){
  lollipopPlot(
  maf = maf_object_shared,
  gene = x,
  showMutationRate = TRUE, 
  labelPos = top10_gene_list_position[[x]]
)
})
dev.off()

### Genearte Oncoplot for CNV

#### Read in CNV files 
pnoc008_cnv <- readRDS(file.path(pnoc008_dir, "pnoc008_cnv_filtered.rds"))
tcga_cnv <- readRDS(file.path(adult_dir, "tcga_gbm_cnv_filtered.rds"))

#### Prepare combined MAF for recurrent CNV
recurrent_subjects <- mutational_analysis$recurrent_alterations %>% 
  filter(Alteration_Datatype == "CNV") %>% 
  pull(subject_id) %>% 
  append(pnoc008_sample_of_interest) %>%
  unique()

recurrent_subjects_pnoc008 <- recurrent_subjects[recurrent_subjects %in% consensus_maf_pnoc008$subjectID]
recurrent_subjects_tcga <- recurrent_subjects[!recurrent_subjects %in% consensus_maf_pnoc008$subjectID]

# Use the kids_first_biospecimen_id to filter TCGA combined consensus SNV
matched_maf_tcga_recurrent <- tcga_maf %>% 
  filter(Tumor_Sample_Barcode %in% recurrent_subjects_tcga)

# filter based on sample_name
matched_maf_pnoc008_recurrent  <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% recurrent_subjects_pnoc008)

# combine the maf to prepare for MAF object
matched_maf_tcga_recurrent <- matched_maf_tcga_recurrent %>% 
  select(common_columns)
matched_maf_pnoc008_recurrent <- matched_maf_pnoc008_recurrent %>% 
  select(common_columns)

combined_maf_recurrent <- matched_maf_tcga_recurrent %>%
  rbind(matched_maf_pnoc008_recurrent) %>%
  dplyr::rename(AAChange = HGVSp_Short) 

#### Filter CNV files to contain only recurrent specimens
pnoc008_cnv_recurrent <- pnoc008_cnv %>% 
  filter(SampleID %in% recurrent_subjects_pnoc008)
  
tcga_cnv_recurrent <- tcga_cnv %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_subjects_tcga)

combined_cnv_recurrent <- rbind(pnoc008_cnv_recurrent, tcga_cnv_recurrent)

#### Clean up combined CNV recurrent 
combined_cnv_recurrent <- combined_cnv_recurrent %>% 
  dplyr::mutate(CN = case_when(
    Alteration_Type == "Gain" ~ "ShallowAmp", 
    Alteration_Type == "Amplification" ~ "Amp", 
    Alteration_Type == "Loss" ~ "Del"
  )) %>% 
  dplyr::select(-Alteration_Type) %>%
  dplyr::rename(Sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Gene, Sample_name, CN)

#### Generate MAF object with CNV data 
maf_object_recurrent_cnv = read.maf(maf = combined_maf_recurrent, 
                                    cnTable = combined_cnv_recurrent)

#### Plot Oncoplot 
pdf(file.path(output_dir, "mutational_cnv_recurrent_adult.pdf"))
oncoplot(maf = maf_object_recurrent_cnv, top=10)
dev.off()

#### Prepare combined MAF for shared CNV
shared_subjects <- mutational_analysis$shared_genes %>% 
  filter(Alteration_Datatype == "CNV") %>% 
  pull(subject_id) %>% 
  append(pnoc008_sample_of_interest) %>%
  unique() 

shared_subjects_pnoc008 <- shared_subjects[shared_subjects %in% consensus_maf_pnoc008$sample_name]
shared_subjects_tcga <- shared_subjects[!shared_subjects %in% consensus_maf_pnoc008$subjectID]

# Use the kids_first_biospecimen_id to filter TCGA combined consensus SNV
matched_maf_tcga_shared <- tcga_maf %>% filter(Tumor_Sample_Barcode %in% shared_subjects_tcga)

# filter based on sample_name
matched_maf_pnoc008_shared  <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% shared_subjects_pnoc008)

# combine the maf to prepare for MAF object
matched_maf_tcga_shared <- matched_maf_tcga_shared %>% 
  select(common_columns)
matched_maf_pnoc008_shared <- matched_maf_pnoc008_shared %>% 
  select(common_columns)

combined_maf_shared <- matched_maf_tcga_shared %>%
  rbind(matched_maf_pnoc008_shared) %>%
  dplyr::rename(AAChange = HGVSp_Short) 

#### Filter CNV files to contain only recurrent specimens
pnoc008_cnv_shared <- pnoc008_cnv %>% 
  filter(SampleID %in% shared_subjects_pnoc008 )
  
tcga_cnv_shared <- tcga_cnv %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_subjects_tcga)

combined_cnv_shared <- rbind(pnoc008_cnv_shared, tcga_cnv_shared)

#### Clean up combined CNV in shared gene samples
combined_cnv_shared <- combined_cnv_shared %>% 
  dplyr::mutate(CN = case_when(
    Alteration_Type == "Gain" ~ "ShallowAmp", 
    Alteration_Type == "Amplification" ~ "Amp", 
    Alteration_Type == "Loss" ~ "Del"
  )) %>% 
  dplyr::select(-Alteration_Type) %>%
  dplyr::rename(Sample_name = Kids_First_Biospecimen_ID) %>%
  dplyr::select(Gene, Sample_name, CN)

#### Generate MAF object with CNV data 
maf_object_shared_cnv = read.maf(maf = combined_maf_shared, cnTable = combined_cnv_shared)

#### Plot Oncoplot 
pdf(file.path(output_dir, "mutational_cnv_shared_adult.pdf"))
oncoplot(maf = maf_object_shared_cnv, top=10)
dev.off()
