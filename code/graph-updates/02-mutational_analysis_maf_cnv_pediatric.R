suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
  library(optparse)
})

option_list <- list(
  make_option(c("--patient"), type = "character",
              help = "Patient identifier, i.e. PNOC008-XX")
)

# parameters to pass
opt <- parse_args(OptionParser(option_list = option_list))
pnoc008_sample_of_interest <- opt$patient

# Define directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
ref_dir <- file.path(root_dir, "data", "reference")
pediatric_dir <- file.path(ref_dir, "pbta")
pnoc008_dir <- file.path(ref_dir, "pnoc008")
output_dir <- file.path(root_dir, "results", pnoc008_sample_of_interest, "output")

## P5 Pediatric Tumor Analysis - Mutational Analysis 
### Read in the data files
# read in PNOC008 combined maf 
consensus_maf_pnoc008 <- readRDS(file.path(pnoc008_dir, "pnoc008_consensus_mutation.rds"))

# read in PBTA maf 
consensus_maf_pbta <- readr::read_tsv(file.path(pediatric_dir, "pbta-snv-consensus-mutation.maf.tsv.gz"))

# Since the MAF object requires histology information, histology files for PNOC008 and pbta will be harmonized to a combined one - histology files downloaded from (s3://d3b-bix-dev-data-bucket/PNOC008/reference/pbta/pbta-histologies.tsv; s3://d3b-bix-dev-data-bucket/PNOC008/reference/pnoc008/pnoc008_clinical.rds)

# generate a combined histology for PNOC008 and PBTA
histology <- readr::read_tsv(file.path(pediatric_dir, "pbta-histologies.tsv"))
histology_pnoc008 <- readRDS(file.path(pnoc008_dir, "pnoc008_clinical.rds"))

histology_pnoc008 <- histology_pnoc008 %>% dplyr::rename(
  Participant_ID = subjectID,
  harmonized_diagnosis = tumorType,
  primary_site = tumorLocation,
  age_at_diagnosis_days = age_diagnosis_days, 
  reported_gender = sex, 
  cohort = study_id
) %>% select(-age_collection_days) %>% 
  select(-library_name)

histology <-  histology %>% 
  dplyr::rename(Participant_ID = Kids_First_Participant_ID) %>%
  select(colnames(histology_pnoc008))

### Generate plots for recurrent alterations
#### Filter combined consensus MAF to recurrent alternated samples 
mutational_analysis <- readRDS(file.path(output_dir, "mutational_analysis_pediatric.rds"))

recurrent_alterations_all <- mutational_analysis$recurrent_alterations %>% 
  filter(Alteration_Datatype == "Mutation") %>% 
  pull(kids_first_biospecimen_id) %>% 
  unique() 


recurrent_alterations_pbta <- recurrent_alterations_all[recurrent_alterations_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]
recurrent_alterations_pnoc008 <- recurrent_alterations_all[!recurrent_alterations_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]

# Use subject ID to filter pnoc008 combined consensus SNV
pnoc008_subject <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_alterations_pnoc008) %>%
  pull(Participant_ID) %>% 
  unique() %>%
  append(pnoc008_sample_of_interest)

# Check to see whether all the non-pbta subject are in consensus maf 
stopifnot(all(pnoc008_subject %in% consensus_maf_pnoc008$sample_name))

# Use the kids_first_biospecimen_id to filter pbta combined consensus SNV
matched_maf_pbta <- consensus_maf_pbta %>% 
  filter(Tumor_Sample_Barcode %in% recurrent_alterations_pbta)

# filter based on sample_name
matched_maf_pnoc008 <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% pnoc008_subject)

# merge the maf into one combined file
common_columns <- intersect(colnames(matched_maf_pbta), colnames(matched_maf_pnoc008))
matched_maf_pbta <- matched_maf_pbta %>% select(common_columns)
matched_maf_pnoc008 <- matched_maf_pnoc008 %>% select(common_columns)
combined_maf <- rbind(matched_maf_pbta, matched_maf_pnoc008)

#### Fix the histology for PNOC008 to use DNA samples as `Kids_First_Biospecimen_ID`
match_pnoc008 <- consensus_maf_pnoc008 %>% 
  select(sample_name, Tumor_Sample_Barcode) %>%
  distinct() %>%
  dplyr::rename(Participant_ID = sample_name)

histology_pnoc008_recurrent <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_alterations_pnoc008) %>%
  left_join(match_pnoc008) %>%
  select(-Kids_First_Biospecimen_ID)
  
histology_recurrent <- histology %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_alterations_pbta) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

dna_matched_histology <- rbind(histology_pnoc008_recurrent, histology_recurrent)

#### Clean up the MAF file and generate a MAF object
combined_maf  <- combined_maf %>% 
  dplyr::rename(AAChange = HGVSp_Short) 

# only certain variant classifications are allowed
# https://github.com/PoisonAlien/maftools/blob/master/R/lollipopPlot.R#L130-L131
combined_maf <- combined_maf %>%
  filter(Variant_Classification %in% c("Nonstop_Mutation", "Missense_Mutation",
                                       "Nonsense_Mutation", "Splice_Site",
                                       "Frame_Shift_Ins", "Frame_Shift_Del",
                                       "In_Frame_Del", "In_Frame_Ins"))

maf_object <- read.maf(maf = combined_maf , clinicalData = dna_matched_histology)

#### Plot MAF Summary 
pdf(file = file.path(output_dir, "mutational_recurrent_pediatric.pdf"))
plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#### Generate lollipop plot 
# First find the ten most mutated gene
summary_genes <- getGeneSummary(maf_object) %>% 
  arrange(desc(AlteredSamples)) 

# find the most mutated 10 genes
top10_gene_list <- summary_genes %>%
  head(10) %>% 
  pull(Hugo_Symbol) %>% 
  unique()

top10_gene_list_position <- lapply(top10_gene_list, function(x) {
  protein_changes <- combined_maf %>% 
    filter(Hugo_Symbol == x) %>% 
    filter(AAChange != ".") %>% 
    pull(AAChange)
  positions <- gsub(".*?([0-9]+).*", "\\1", protein_changes)
})

names(top10_gene_list_position) <- top10_gene_list 

pdf(file.path(output_dir, "lollipop_recurrent_pediatric.pdf"), onefile = TRUE)
lapply(top10_gene_list, function(x){
  lollipopPlot(
  maf = maf_object,
  gene = x,
  showMutationRate = TRUE, 
  labelPos = top10_gene_list_position[[x]]
)
})
dev.off()

### Generate plots for shared mutations 

#### Filter combined consensus MAF to shared genes samples 
shared_genes_all <- mutational_analysis$shared_genes %>% 
  filter(Alteration_Datatype == "Mutation") %>% 
  pull(kids_first_biospecimen_id) %>% unique() 

shared_genes_pbta <- shared_genes_all[shared_genes_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]
shared_genes_pnoc008 <- shared_genes_all[!shared_genes_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]

# Use the kids_first_biospecimen_id to filter pbta combined consensus SNV
matched_maf_pbta_shared <- consensus_maf_pbta %>% 
  filter(Tumor_Sample_Barcode %in% shared_genes_pbta)

# For PNOC008 samples, we need to match back to the subject ID to filter for the MAF file 
# Use subject ID to filter pnoc008 combined consensus SNV
pnoc008_subject <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_genes_pnoc008) %>%
  pull(Participant_ID) %>% unique() %>%
  append(pnoc008_sample_of_interest)

# Check to see whether all the non-pbta subject are in consensus maf 
stopifnot(all(pnoc008_subject %in% consensus_maf_pnoc008$sample_name))

# filter based on sample_name
matched_maf_pnoc008_shared <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% pnoc008_subject)

# merge the maf into one combined file
matched_maf_pbta_shared <- matched_maf_pbta_shared %>% select(common_columns)
matched_maf_pnoc008_shared <- matched_maf_pnoc008_shared %>% select(common_columns)
combined_maf_shared <- rbind(matched_maf_pbta_shared, matched_maf_pnoc008_shared)

#### Fix the histology for PNOC008 to use DNA samples as `Kids_First_Biospecimen_ID`
histology_pnoc008_shared <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_genes_pnoc008) %>%
  left_join(match_pnoc008) %>%
  select(-Kids_First_Biospecimen_ID)
  
histology_shared <- histology %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_genes_pbta) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

dna_matched_histology_shared <- rbind(histology_pnoc008_shared, histology_shared)


#### Clean up the MAF file and generate a MAF object
combined_maf_shared  <- combined_maf_shared %>% 
  dplyr::rename(AAChange = HGVSp_Short) 
maf_object <- read.maf(maf = combined_maf_shared , clinicalData = dna_matched_histology_shared)

#### Plot MAF Summary 
pdf(file.path(output_dir, "mutational_shared_pediatric.pdf"))
plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#### Generate lollipop plot 
# First find the ten most mutated gene
summary_genes <- getGeneSummary(maf_object) %>% 
  arrange(desc(AlteredSamples)) 

# find the most mutated 10 genes
top10_gene_list <- summary_genes %>%
  head(10) %>% 
  pull(Hugo_Symbol) %>% unique()

top10_gene_list_position <- lapply(top10_gene_list, function(x) {
  protein_changes <- combined_maf_shared %>% 
    filter(Hugo_Symbol == x) %>% 
    filter(AAChange != ".") %>% 
    pull(AAChange)
  positions <- gsub(".*?([0-9]+).*", "\\1", protein_changes)
})

names(top10_gene_list_position) <- top10_gene_list 

pdf(file.path(output_dir, "lollipop_shared_pediatric.pdf"))
lapply(top10_gene_list, function(x){
  lollipopPlot(
  maf = maf_object,
  gene = x,
  showMutationRate = TRUE, 
  labelPos = top10_gene_list_position[[x]]
)
})
dev.off()

### Generate oncoplots for CNV in recurrent alteration

#### First read in the CNV files
pnoc008_cnv <- readRDS(file.path(pnoc008_dir, "pnoc008_cnv_filtered.rds"))
pbta_cnv <- readRDS(file.path(pediatric_dir, "pbta-cnv-cnvkit-filtered.rds"))

#### Filter combined consensus MAF to recurrent CNV samples 
recurrent_cnv_all <- mutational_analysis$recurrent_alterations %>% 
  filter(Alteration_Datatype == "CNV") %>% 
  pull(kids_first_biospecimen_id) %>% 
  unique() 

recurrent_cnv_pbta <- recurrent_cnv_all[recurrent_cnv_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]
recurrent_cnv_pnoc008 <- recurrent_cnv_all[!recurrent_cnv_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]

# Use the kids_first_biospecimen_id to filter pbta combined consensus SNV
matched_maf_pbta_recurrent_cnv <- consensus_maf_pbta %>% 
  filter(Tumor_Sample_Barcode %in% recurrent_cnv_pbta)

# For PNOC008 samples, we need to match back to the subject ID to filter for the MAF file 
# Use subject ID to filter pnoc008 combined consensus SNV
pnoc008_subject <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_cnv_pnoc008) %>%
  pull(Participant_ID) %>% unique() %>%
  append(pnoc008_sample_of_interest)

# filter based on sample_name
matched_maf_pnoc008_recurrent_cnv <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% pnoc008_subject)

# merge the maf into one combined file
matched_maf_pbta_recurrent_cnv <- matched_maf_pbta_recurrent_cnv %>% 
  select(common_columns)
matched_maf_pnoc008_recurrent_cnv <- matched_maf_pnoc008_recurrent_cnv %>%
  select(common_columns)

combined_maf_recurrent_cnv <- rbind(matched_maf_pbta_recurrent_cnv, matched_maf_pnoc008_recurrent_cnv)

#### Filter CNV files to contain only recurrent specimens
pnoc008_cnv_recurrent <- pnoc008_cnv %>% 
  filter(SampleID %in% pnoc008_subject)
  
pbta_cnv_recurrent <- pbta_cnv %>% 
  filter(Kids_First_Biospecimen_ID %in% recurrent_cnv_pbta)

combined_cnv_recurrent <- rbind(pnoc008_cnv_recurrent, pbta_cnv_recurrent)

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
maf_object_recurrent_cnv = read.maf(maf = combined_maf_recurrent_cnv, cnTable = combined_cnv_recurrent)

#### Plot Oncoplot 
pdf(file.path(output_dir, "mutational_cnv_recurrent_pediatric.pdf"))
oncoplot(maf = maf_object_recurrent_cnv, top=10)
dev.off()

### Generate oncoplots for CNV in shared genes
#### Filter combined consensus MAF to shared CNV samples 
shared_cnv_all <- mutational_analysis$shared_genes %>% 
  filter(Alteration_Datatype == "CNV") %>% 
  pull(kids_first_biospecimen_id) %>% unique() 

shared_cnv_pbta <- shared_cnv_all[shared_cnv_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]
shared_cnv_pnoc008 <- shared_cnv_all[!shared_cnv_all %in% consensus_maf_pbta$Tumor_Sample_Barcode]

# Use the kids_first_biospecimen_id to filter pbta combined consensus SNV
matched_maf_pbta_shared_cnv <- consensus_maf_pbta %>% 
  filter(Tumor_Sample_Barcode %in% shared_cnv_pbta)

# For PNOC008 samples, we need to match back to the subject ID to filter for the MAF file 
# Use subject ID to filter pnoc008 combined consensus SNV
pnoc008_subject <- histology_pnoc008 %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_cnv_pnoc008 ) %>%
  pull(Participant_ID) %>% unique() %>%
  append(pnoc008_sample_of_interest)

# filter based on sample_name
matched_maf_pnoc008_shared_cnv <- consensus_maf_pnoc008 %>% 
  filter(sample_name %in% pnoc008_subject)

# merge the maf into one combined file
matched_maf_pbta_shared_cnv <- matched_maf_pbta_shared_cnv %>% 
  select(common_columns)
matched_maf_pnoc008_shared_cnv <- matched_maf_pnoc008_shared_cnv %>% 
  select(common_columns)
combined_maf_shared_cnv <- rbind(matched_maf_pbta_shared_cnv, matched_maf_pnoc008_shared_cnv)

#### Filter CNV files to contain only recurrent specimens
pnoc008_cnv_shared <- pnoc008_cnv %>% 
  filter(SampleID %in% pnoc008_subject)
  
pbta_cnv_shared <- pbta_cnv %>% 
  filter(Kids_First_Biospecimen_ID %in% shared_cnv_pbta)

combined_cnv_shared <- rbind(pnoc008_cnv_shared, pbta_cnv_shared)

#### Clean up combined CNV recurrent 
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
maf_object_shared_cnv = read.maf(maf = combined_maf_shared_cnv, cnTable = combined_cnv_shared)

#### Plot Oncoplot 
pdf(file.path(output_dir, "mutational_cnv_shared_pediatric.pdf"))
oncoplot(maf = maf_object_shared_cnv, top=10)
dev.off()
