suppressPackageStartupMessages({
  library(dplyr)
  library(maftools)
})

recurrent_alterations_plots <- function(ref_cancer_dir, patient_maf = filtered_maf, patient_cnv = filtered_cnv, mutational_analysis_output, prefix){
  
  # ref tumors maf
  ref_maf <- list.files(path = ref_cancer_dir, pattern = "mutation_filtered.rds", full.names = T)
  ref_maf <- readRDS(ref_maf)
  
  # common columns
  common_cols <- intersect(colnames(ref_maf), colnames(patient_maf))
  
  # create combined maf
  combined_maf <- patient_maf %>%
    dplyr::select(common_cols) %>%
    rbind(ref_maf %>%
            dplyr::select(common_cols)) %>%
    unique()
  
  # recurrent mutations
  recurrent_mutations_samples <- mutational_analysis_output$recurrent_alterations %>% 
    filter(Alteration_Datatype == "Mutation") %>% 
    pull(Kids_First_Biospecimen_ID) %>% 
    unique() 
  
  # filter combined maf to samples with recurrent mutations 
  combined_maf_recurrent <- combined_maf %>%
    filter(Tumor_Sample_Barcode %in% recurrent_mutations_samples) %>%
    dplyr::rename(AAChange = HGVSp_Short)
  
  # only certain variant classifications are allowed
  # https://github.com/PoisonAlien/maftools/blob/master/R/lollipopPlot.R#L130-L131
  combined_maf_recurrent <- combined_maf_recurrent %>%
    filter(Variant_Classification %in% c("Nonstop_Mutation", "Missense_Mutation",
                                         "Nonsense_Mutation", "Splice_Site",
                                         "Frame_Shift_Ins", "Frame_Shift_Del",
                                         "In_Frame_Del", "In_Frame_Ins"))
  
  maf_object <- read.maf(maf = combined_maf_recurrent)
  
  # 1) plot maf summary 
  pdf(file = file.path(output_dir, paste0("recurrent_mutations_", prefix, ".pdf")))
  plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  # 2) lollipop plot 
  # find the ten most mutated gene
  top10_gene_list <- getGeneSummary(maf_object) %>% 
    arrange(desc(AlteredSamples)) %>%
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
  
  pdf(file.path(output_dir, paste0("recurrent_mutations_lollipop_", prefix, ".pdf")))
  lapply(top10_gene_list, function(x){
    lollipopPlot(
      maf = maf_object,
      gene = x,
      showMutationRate = TRUE, 
      labelPos = top10_gene_list_position[[x]]
    )
  })
  dev.off()
  
  # 3) oncoplots for CNV in recurrent alteration
  # recurrent copy number changes
  recurrent_cnv_samples <- mutational_analysis_output$recurrent_alterations %>% 
    filter(Alteration_Datatype == "CNV") %>% 
    pull(Kids_First_Biospecimen_ID) %>% 
    unique() 
  
  # ref tumors cnv
  ref_cnv <- list.files(path = ref_cancer_dir, pattern = "cnv_filtered.rds", full.names = T)
  ref_cnv <- readRDS(ref_cnv)
  ref_cnv <- ref_cnv %>%
    mutate(CN = case_when(
      status == "Gain" ~ "ShallowAmp", 
      status == "Amplification" ~ "Amp", 
      status == "Loss" ~ "Del"),
      Gene = hgnc_symbol,
      Sample_name = Kids_First_Biospecimen_ID) %>%
    dplyr::select(Gene, Sample_name, CN) %>%
    unique()
  
  # patient cnv
  patient_cnv <- patient_cnv %>%
    mutate(CN = case_when(
      status == "Gain" ~ "ShallowAmp", 
      status == "Amplification" ~ "Amp", 
      status == "Loss" ~ "Del"),
      Gene = hgnc_symbol,
      Sample_name = Kids_First_Biospecimen_ID) %>%
    dplyr::select(Gene, Sample_name, CN) %>%
    unique()
  
  # create combined cnv
  combined_cnv <- patient_cnv %>% 
    rbind(ref_cnv) %>%
    unique()
  
  # filter combined cnv to samples with recurrent cnv 
  combined_cnv_recurrent <- combined_cnv %>%
    filter(Sample_name %in% recurrent_mutations_samples) 
  
  # filter combined maf to samples with recurrent cnv
  combined_maf_matched_recurrent_cnv <- combined_maf %>%
    filter(Tumor_Sample_Barcode %in% recurrent_mutations_samples) 
  maf_object_recurrent_cnv <- read.maf(maf = combined_maf_matched_recurrent_cnv, cnTable = combined_cnv_recurrent)
  
  # plot oncoplot 
  pdf(file.path(output_dir, paste0("recurrent_mutational_cnv_", prefix, ".pdf")))
  oncoplot(maf = maf_object_recurrent_cnv, top = 10)
  dev.off()
}
