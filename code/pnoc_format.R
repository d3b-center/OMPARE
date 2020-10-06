# Author: Komal S. Rathi
# Function: script to pull rsem files from all existing patients and do the following: 
# create a matrix of expression
# create metadata using clinical files
# create full summary data files of cnv, mutations and fusions
# This is to be run everytime a new patient comes in - before generating the report

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(GenomicRanges))

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
source(file.path(root_dir, "code", "utils", "define_directories.R"))

# create output directory
pnoc008.dir <- file.path(ref_dir, 'PNOC008')
dir.create(pnoc008.dir, showWarnings = F, recursive = T)

# source functions
source(file.path(utils_dir, 'createCopyNumber.R'))

# reference data
# cancer genes 
cancerGenes <- readRDS(file.path(ref_dir, 'cancer_gene_list.rds'))
gene.list <- unique(cancerGenes$Gene_Symbol)

# chr coordinates to gene symbol map
chrMap <- read.delim(file.path(ref_dir, "mart_export_genechr_mapping.txt"), stringsAsFactors =F)
colnames(chrMap) <- c("hgnc_symbol", "gene_start", "gene_end", "chromosome")

# function to merge expression
merge.res <- function(nm){
  sample_name <- gsub(".*PNOC","PNOC",nm)
  sample_name <- gsub('/.*','',sample_name)
  x <- data.table::fread(nm)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# function to merge degs
merge.excel <- function(nm){
  sample_name <- gsub(".*PNOC", "PNOC", nm)
  sample_name <- gsub('_.*','',sample_name)
  x <- readxl::read_xlsx(path = nm, sheet = 'DE_Genes_Up')
  y <- readxl::read_xlsx(path = nm, sheet = 'DE_Genes_Down')
  x <- x %>%
    filter(Comparison == "GTExBrain_1152") %>%
    dplyr::select(Gene_name, DE)
  y <- y %>%
    filter(Comparison == "GTExBrain_1152") %>%
    dplyr::select(Gene_name, DE)
  x <- rbind(x, y)
  if(nrow(x) > 1){
    x <- as.data.frame(x)
    x$sample_name <- sample_name
    return(x)
  } 
}

# function to read cnv, filter by genes and merge
merge.cnv <- function(cnvData, genelist){
  # PNOC
  sample_name <- gsub(".*PNOC", "PNOC", cnvData)
  sample_name <- gsub('/.*', '', sample_name)
  cnvData <- data.table::fread(cnvData, header = T, check.names = T)
  ploidy <- NULL
  cnvData <- cnvData %>% 
    dplyr::select(chr, start, end, copy.number, status, WilcoxonRankSumTestPvalue) %>%
    filter(WilcoxonRankSumTestPvalue < 0.05) %>%
    as.data.frame()
  cnvOut <- createCopyNumber(cnvData = cnvData, ploidy = ploidy) # map coordinates to gene symbol
  cnvOut <- cnvOut %>%
    filter(Gene %in% genelist,
           Status != "neutral") %>% # filter to gene list
    mutate(sample_name = sample_name) # add PNOC008 patient id
  return(cnvOut)
}

# list of all PNOC patients
pat.expDat <- list.files(path = data_dir, pattern = "*.genes.results*", recursive = TRUE, full.names = T)
pat.expDat <- pat.expDat[grep('PNOC008-',  pat.expDat)]
pat.expr.mat <- lapply(pat.expDat, FUN = function(x) merge.res(x))
pat.expr.mat <- data.table::rbindlist(pat.expr.mat)

# separate gene_id and gene_symbol
pat.expr.mat <- pat.expr.mat %>% 
  mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
  separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
  unique()

# uniquify gene_symbol
pat.expr.mat <- pat.expr.mat %>% 
  group_by(sample_name) %>%
  arrange(desc(TPM)) %>% 
  distinct(gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(gene_symbol, sample_name, TPM) %>%
  spread(sample_name, TPM) %>%
  column_to_rownames('gene_symbol')

# only keep NANT sample for PNOC008-5
pat.expr.mat <- pat.expr.mat[,grep('CHOP', colnames(pat.expr.mat), invert = T)]
colnames(pat.expr.mat)  <- gsub("-NANT", "", colnames(pat.expr.mat))

# now merge all clinical data for all patients
pat.clinData <- readxl::read_xlsx(file.path(ref_dir, 'Manifest' , 'PNOC008_Manifest.xlsx'), sheet = 1)
colnames(pat.clinData) <- gsub('[() ]', '.', colnames(pat.clinData))
pat.clinData <- pat.clinData %>%
  filter_all(any_vars(!is.na(.))) %>%
  mutate(PNOC.Subject.ID = gsub('P-','PNOC008-', PNOC.Subject.ID))
pat.clinData <- pat.clinData %>%
  mutate(study_id = "PNOC008",
         KF_ParticipantID = PNOC.Subject.ID,
         subjectID = PNOC.Subject.ID,
         reportDate = Sys.Date(),
         tumorType = Diagnosis.a,
         tumorLocation = Primary.Site.a,
         ethnicity = Ethnicity,
         age_diagnosis_days = Age.at.Diagnosis..in.days.,
         age_collection_days = Age.at.Collection..in.days.,
         sex = Gender) %>%
  dplyr::select(subjectID, KF_ParticipantID, tumorType, tumorLocation, ethnicity, sex, age_diagnosis_days, age_collection_days, study_id, library_name) %>%
  as.data.frame()
rownames(pat.clinData) <- pat.clinData$subjectID

common.pat <- intersect(rownames(pat.clinData), colnames(pat.expr.mat))
pat.clinData <- pat.clinData[common.pat,]
pat.expr.mat <- pat.expr.mat[,common.pat]

# save expression and clinical
saveRDS(pat.expr.mat, file = file.path(pnoc008.dir, 'PNOC008_TPM_matrix.RDS'))
saveRDS(pat.clinData, file = file.path(pnoc008.dir, "PNOC008_clinData.RDS"))

# copy number
cnv.files <- list.files(path = data_dir, pattern = "*.CNVs.p.value.txt", recursive = TRUE, full.names = T)
cnv.files <- cnv.files[grep('PNOC008-',  cnv.files)]
pnoc.cnv <- lapply(cnv.files, FUN = function(x) merge.cnv(cnvData = x, genelist = gene.list))
pnoc.cnv <- data.table::rbindlist(pnoc.cnv)
# only keep NANT sample for PNOC008-5
pnoc.cnv <- pnoc.cnv[grep('CHOP', pnoc.cnv$sample_name, invert = T),]
pnoc.cnv$sample_name  <- gsub("-NANT", "", pnoc.cnv$sample_name)

pnoc.cnv <- pnoc.cnv %>%
  mutate(Alteration_Datatype = "CNV",
         Alteration_Type = stringr::str_to_title(Status),
         Alteration = paste0('Copy Number Value:', CNA),
         Kids_First_Biospecimen_ID = sample_name,
         SampleID = sample_name,
         Study = "PNOC008") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pnoc.cnv, file = file.path(pnoc008.dir, "PNOC008_cnvData_filtered.rds"))

# fusions
fusion.files <- list.files(path = data_dir, pattern = "*.arriba.fusions.tsv", recursive = TRUE, full.names = T)
fusion.files <- fusion.files[grep('PNOC008-',  fusion.files)]
pnoc.fusions <- lapply(fusion.files, FUN = function(x) merge.res(x))
pnoc.fusions <- data.table::rbindlist(pnoc.fusions)
colnames(pnoc.fusions)[1] <- "gene1"
# only keep NANT sample for PNOC008-5
pnoc.fusions <- pnoc.fusions[grep('CHOP', pnoc.fusions$sample_name, invert = T),]
pnoc.fusions$sample_name  <- gsub("-NANT", "", pnoc.fusions$sample_name)

pnoc.fusions <- pnoc.fusions %>%
  separate_rows(gene1, gene2, sep = ",", convert = TRUE) %>%
  mutate(gene1 = gsub('[(].*', '', gene1),
         gene2 = gsub('[(].*',' ', gene2),
         reading_frame = ifelse(reading_frame == ".", "other", reading_frame)) %>%
  mutate(Alteration_Datatype = "Fusion",
         Alteration_Type = stringr::str_to_title(reading_frame),
         Alteration = paste0(gene1, '_',  gene2),
         Kids_First_Biospecimen_ID = sample_name,
         SampleID = sample_name,
         Study = "PNOC008") %>%
  unite(col = "Gene", gene1, gene2, sep = ", ", na.rm = T) %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  separate_rows(Gene, convert = TRUE) %>%
  filter(Gene %in% gene.list) %>%
  unique()
saveRDS(pnoc.fusions, file = file.path(pnoc008.dir, "PNOC008_fusData_filtered.rds"))

# mutations
mut.files <- list.files(path = data_dir, pattern = "*consensus_somatic.vep.maf", recursive = TRUE, full.names = T)
mut.files <- mut.files[grep('PNOC008-',  mut.files)]
pnoc.mutations <- lapply(mut.files, FUN = function(x) merge.res(x))
pnoc.mutations <- data.table::rbindlist(pnoc.mutations)
# only keep NANT sample for PNOC008-5
pnoc.mutations <- pnoc.mutations[grep('CHOP', pnoc.mutations$sample_name, invert = T),]
pnoc.mutations$sample_name  <- gsub("-NANT", "", pnoc.mutations$sample_name)

keepVC <- c("Nonsense_Mutation", "Missense_Mutation", 
            "Splice_Region", "Splice_Site",
            "3'UTR", "5'UTR", 
            "In_Frame_Ins", "In_Frame_Del",
            "Frame_Shift_Ins", "Frame_Shift_Del")
keepVI <- c("MODIFIER", "MODERATE", "HIGH")
pnoc.mutations <- pnoc.mutations %>%
  filter(BIOTYPE == "protein_coding",
         Variant_Classification %in% keepVC,
         IMPACT %in% keepVI,
         Hugo_Symbol %in% gene.list) %>%
  mutate(Gene = Hugo_Symbol,
         Alteration_Datatype = "Mutation",
         Alteration_Type = Variant_Classification,
         Alteration = HGVSp_Short,
         Kids_First_Biospecimen_ID = sample_name,
         SampleID = sample_name,
         Study = "PNOC008") %>%
  dplyr::select(Gene, Alteration_Datatype, Alteration_Type, Alteration, Kids_First_Biospecimen_ID, SampleID, Study) %>%
  unique()
saveRDS(pnoc.mutations, file = file.path(pnoc008.dir, "PNOC008_consensus_mutData_filtered.rds"))

# cohort level degs
deg.files <- list.files(path = data_dir, pattern = "*_summary.xlsx", recursive = TRUE, full.names = T)
deg.files <- deg.files[grep('PNOC008-',  deg.files)]
pnoc.deg <- lapply(deg.files, FUN = function (x) merge.excel(x))
pnoc.deg <- data.table::rbindlist(pnoc.deg)
# only keep NANT sample for PNOC008-5
pnoc.deg <- pnoc.deg[grep('CHOP', pnoc.deg$sample_name, invert = T),]
pnoc.deg$sample_name  <- gsub("-NANT", "", pnoc.deg$sample_name)
saveRDS(pnoc.deg, file = file.path(pnoc008.dir, "PNOC008_deg_GTExBrain.rds"))

# cohort level tmb scores
TMBFileBED <- data.table::fread(file.path(ref_dir, "xgen-exome-research-panel-targets_hg38_ucsc_liftover.100bp_padded.sort.merged.bed"))
colnames(TMBFileBED)  <- c("chr", "start", "end")

# read mutect2 for TMB profile
mutFiles <- list.files(path = data_dir, pattern = 'mutect2_somatic.vep.maf', recursive = TRUE, full.names = T)
mutFiles <- mutFiles[grep('PNOC008-',  mutFiles)]
pnoc.mutations <- lapply(mutFiles, FUN = function(x) merge.res(x))
pnoc.mutations <- data.table::rbindlist(pnoc.mutations, fill = T)
# only keep NANT sample for PNOC008-5
pnoc.mutations <- pnoc.mutations[grep('CHOP', pnoc.mutations$sample_name, invert = T),]
pnoc.mutations$sample_name  <- gsub("-NANT", "", pnoc.mutations$sample_name)

# mutect2 nonsense and missense mutations
pnoc.mutations <- pnoc.mutations %>%
  filter(Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation")) %>%
  dplyr::select(sample_name, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position) %>%
  unique()
  
# intersect with bed file
subject <- with(TMBFileBED, GRanges(chr, IRanges(start = start, end = end)))
query <- with(pnoc.mutations, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = Hugo_Symbol)))
tmb.res <- findOverlaps(query = query, subject = subject, type = "within")
tmb.res <- data.frame(pnoc.mutations[queryHits(tmb.res),], TMBFileBED[subjectHits(tmb.res),])
  
# return the number of missense + nonsense overlapping with the bed file/77.46
tmb.res <- tmb.res %>%
  group_by(sample_name) %>%
  mutate(num.mis.non = n()) %>%
  dplyr::select(sample_name, num.mis.non)  %>%
  unique() %>%
  mutate(tmb = num.mis.non/77.46) %>%
  dplyr::select(sample_name, tmb)
saveRDS(tmb.res, file = file.path(pnoc008.dir, "PNOC008_TMBscores.rds"))



