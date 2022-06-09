# tumor inflammation signature
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# protein coding genes and remove HISTONE genes
# read gencode from OpenPedCan
gencode_gtf <- rtracklayer::import(con = file.path(data_dir, 'OpenPedCan-analysis/data/gencode.v27.primary_assembly.annotation.gtf.gz'))
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_pc <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding",
         !grepl("^HIST", gene_name)) %>%
  unique()

# collapse rnaseq data
source(file.path(root_dir, "code", "utils", "collapse_rnaseq.R"))

tis_profile <- function(pediatric_dir, adult_dir, patient_of_interest, norm_method){
  
  # for this we need the full set of pediatric and adult data
  
  # pediatric tumors from OT (because master genomics does not have short histology information)
  ped_clinical <- read.delim(file.path(data_dir, "OpenPedCan-analysis", "data", "histologies.tsv"))
  ped_clinical <- ped_clinical %>%
    filter(experimental_strategy == "RNA-Seq",
           !cohort %in% c("GTEx", "TCGA"),
           sample_type == "Tumor")
  
  # pnoc tumors
  pnoc_clinical <- read.delim(file.path(data_dir, "pnoc008", "pnoc008_clinical.tsv"))
  ped_clinical <- ped_clinical %>%
    filter(!cohort_participant_id %in% pnoc_clinical$cohort_participant_id)
  common_cols <- intersect(colnames(ped_clinical), colnames(pnoc_clinical))
  ped_clinical <- rbind(pnoc_clinical[,common_cols], ped_clinical[,common_cols])
  
  # add short_histology information from open targets
  ped_clinical <- ped_clinical %>% 
    mutate(type = "Pediatric") %>%
    unique()
  
  # counts
  pediatric_counts <- list.files(path = pediatric_dir, pattern = "expected_count", full.names = T)
  pediatric_counts <- readRDS(pediatric_counts)
  pediatric_counts <- pediatric_counts[grep("^HIST", pediatric_counts$gene_id, invert = T),]
  pediatric_counts <- pediatric_counts %>%
    filter(gsub(".*_", "", gene_id) %in% gencode_pc$gene_name)
  pediatric_counts <- collapse_rnaseq(pediatric_counts)
  
  # use common ids
  common_cols <- intersect(colnames(pediatric_counts), ped_clinical$Kids_First_Biospecimen_ID)
  ped_clinical <- ped_clinical %>%
    filter(Kids_First_Biospecimen_ID %in% common_cols)
  pediatric_counts <- pediatric_counts[,ped_clinical$Kids_First_Biospecimen_ID]
  
  # adult tumors
  # histology
  adult_clinical <- list.files(path = adult_dir, pattern = "histologies", full.names = T)
  adult_clinical <- data.table::fread(adult_clinical)
  adult_clinical <- adult_clinical %>%
    filter(experimental_strategy == "RNA-Seq",
           cohort == "TCGA",
           sample_type == "Tumor") %>%
    mutate(type = "Adult")
  
  # counts
  adult_counts <- readRDS(file.path(adult_dir, "tcga-gene-counts-rsem-expected_count-collapsed.rds"))
  adult_counts <- adult_counts[grep("^HIST", rownames(adult_counts), invert = T),]
  adult_counts <- adult_counts[rownames(adult_counts) %in% gencode_pc$gene_name,]
  
  # use common ids
  common_cols <- intersect(colnames(adult_counts), adult_clinical$Kids_First_Biospecimen_ID)
  adult_clinical <- adult_clinical %>%
    filter(Kids_First_Biospecimen_ID %in% common_cols)
  adult_counts <- adult_counts[,adult_clinical$Kids_First_Biospecimen_ID]
  
  # combine counts on common genes 
  common_genes <- intersect(rownames(adult_counts), rownames(pediatric_counts))
  adult_counts <- adult_counts[common_genes,]
  pediatric_counts <- pediatric_counts[common_genes,]
  merged_counts <- cbind(adult_counts, pediatric_counts)
  
  # combine clinical
  common_cols <- intersect(colnames(adult_clinical), colnames(ped_clinical))
  merged_clinical <- adult_clinical %>%
    select(common_cols) %>%
    rbind(ped_clinical %>%
            select(common_cols))
  merged_clinical <- merged_clinical %>%
    mutate(batch = paste0(cohort,"_", RNA_library),
           type = ifelse(cohort_participant_id %in% patient_of_interest, "Patient", type))
  
  # count number of samples per cancer group and only keep >= 5
  merged_clinical <- merged_clinical  %>%
    group_by(short_histology) %>%
    dplyr::mutate(n = n()) %>%
    filter(n >= 5) %>%
    dplyr::select(-c(n)) %>%
    as.data.frame()
  merged_counts <- merged_counts %>%
    dplyr::select(merged_clinical$Kids_First_Biospecimen_ID)
  
  if(norm_method == "voom"){
    # create design to correct for batch
    var <- factor(merged_clinical$batch)
    design <- model.matrix(~0+var)
    colnames(design) <- levels(var)
    rownames(design) <- merged_clinical$sample_id
    
    # voom normalize
    v <- voom(counts = merged_counts, design = design, plot = FALSE, normalize.method = "quantile")
    norm_dat <- v$E
  } else if(norm_method == "quantile") {
    merged_counts <-  data.matrix(log2(merged_counts + 1))
    norm_dat <- preprocessCore::normalize.quantiles(merged_counts, copy = FALSE)
    norm_dat <- as.data.frame(norm_dat)
  }
  
  # read tumor inflammation signature
  tis <- read.delim(file.path(data_dir, 'tumor_inflammation_signatures.txt'), stringsAsFactors = F)
  
  # filter by genelist and format
  norm_dat <- norm_dat[rownames(norm_dat) %in% tis$Genes,]
  dat_sums <- colSums(norm_dat)
  dat_means <- colMeans(norm_dat)
  tis_output <- data.frame(Kids_First_Biospecimen_ID = names(dat_sums), score_sum = dat_sums, score_avg = dat_means)
  
  # merge with meta file
  tis_output <- merged_clinical %>% 
    dplyr::select(Kids_First_Biospecimen_ID, short_histology, cohort, type) %>% 
    inner_join(tis_output, by = "Kids_First_Biospecimen_ID") %>%
    dplyr::rename("subject_id" = "Kids_First_Biospecimen_ID") %>%
    as.data.frame()
  return(tis_output)
}
