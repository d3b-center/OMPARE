# Author: Komal Rathi
# read patient-specific data and save all data to global env

# libraries
suppressPackageStartupMessages({
  library(dplyr)
})

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
utils_dir <- file.path(root_dir, "code", "utils")
data_dir <- file.path(root_dir, "data")

# source functions
source(file.path(utils_dir, 'filter_fusions.R')) # filter fusions
source(file.path(utils_dir, 'filter_cnv.R')) # filter copy number variations
source(file.path(utils_dir, 'filter_mutations.R')) # filter mutations
source(file.path(utils_dir, 'annotate_mutations.R')) # annotate mutations

# gencode v27 reference
gencode_v27 <- read.delim(file.path(data_dir, 'pnoc008', 'gencode.v27.primary_assembly.annotation.txt'))
gencode_v27_pc <- gencode_v27 %>%
  filter(biotype == "protein_coding")

read_patient_data <- function(patient_dir, fusion_method = c("star","arriba"), snv_caller = "all"){
  
  # patient sample info (at most one with .txt extension)
  sampleInfo <- list.files(path = file.path(patient_dir, 'clinical'), pattern = "patient_report.txt", recursive = T, full.names = T)
  sampleInfo <- read.delim(sampleInfo, stringsAsFactors = F)
  assign("sampleInfo", sampleInfo, envir = globalenv())
  
  # get patient subject id
  patient_of_interest <- sampleInfo$subjectID
  
  # mutation data (can be multiple with .maf extension)
  # if snv_caller is all, then use all files except consensus
  somatic.mut.pattern <- '*.maf'
  mutFiles <- list.files(path = file.path(patient_dir, 'simple-variants'), pattern = somatic.mut.pattern, recursive = TRUE, full.names = T)
  if(snv_caller == "all"){ 
    mutFiles <- grep('consensus', mutFiles, invert = TRUE, value = TRUE)
  } else {
    mutFiles <- grep(snv_caller, mutFiles, value = TRUE)
  }
  if(length(mutFiles) >= 1){
    mutFiles <- lapply(mutFiles, data.table::fread, skip = 1, stringsAsFactors = F)
    mutData <- data.table::rbindlist(mutFiles, fill = TRUE)
    mutData <- as.data.frame(mutData)
    mutData <- unique(mutData)
    assign("mutData", mutData, envir = globalenv())
    
    # filter mutations
    mutDataFilt <- filter_mutations(myMutData = mutData, myCancerGenes = cancer_genes)
    assign("mutDataFilt", mutDataFilt, envir = globalenv())
    
    # annotate mutations
    mutDataAnnot <- annotate_mutations(myMutData = mutData, myCancerGenes = cancer_genes, cancer_hotspots = cancer_hotspots_v2)
    assign("mutDataAnnot", mutDataAnnot, envir = globalenv())
  } 
  
  # use cnv from CNVkit *gainloss.txt
  cnvData <- list.files(path = file.path(patient_dir, 'copy-number-variations'), pattern = "*.gainloss.txt", recursive = TRUE, full.names = T)
  if(length(cnvData) == 1){
    cnvData <- read.delim(cnvData)
    cnvData <- cnvData %>%
      dplyr::rename("hgnc_symbol" = "gene",
                    "chr" = "chromosome") %>%
      dplyr::select(-c(depth, weight, n_bins))
    
    # pull ploidy and purity from controlfreec
    controlfreec_info <- list.files(file.path(patient_dir, 'copy-number-variations'), pattern = "info.txt", full.names = T)
    controlfreec_info <- read.delim(controlfreec_info, header = F)
    ploidy <- controlfreec_info %>% filter(V1 == "Output_Ploidy") %>% .$V2 %>% as.numeric()
    purity <- controlfreec_info %>% filter(V1 == "Sample_Purity") %>% .$V2 %>% as.numeric()
    
    # calculate absolute copy number
    cnvData$copy.number <- 0
    if(ploidy == 2){
      # compute log2 ratio cutoffs
      cutoff <- log2((1 - purity) + purity * (0:3 + .5) / ploidy)
      cutoff <- min(cutoff)

      # compute absolute copy number
      cnvData$copy.number <- (((2^(cnvData$log2)-(1-purity)) * ploidy)/ purity) - 0.5
      cnvData <- cnvData %>%
        rowwise() %>%
        mutate(copy.number = ifelse(log2 < cutoff, round(copy.number), ceiling(copy.number)))
    } else {
      # compute log2 ratio cutoffs
      cutoff <- log2((1 - purity) + purity * (0:6 + .5) / ploidy)
      cutoff <- min(cutoff)

      # compute absolute copy number
      cnvData$copy.number <- (((2^(cnvData$log2)-(1-purity)) * ploidy)/ purity) - 0.5
      cnvData <- cnvData %>%
        rowwise() %>%
        mutate(copy.number = ifelse(log2 < cutoff, round(copy.number), ceiling(copy.number)))
    }
    
    # add copy number status
    cnvData <- cnvData %>%
      mutate(status = case_when(copy.number == 0 ~ "Complete Loss",
                                copy.number < ploidy & copy.number > 0 ~ "Loss",
                                copy.number == ploidy ~ "Neutral",
                                copy.number > ploidy & copy.number < ploidy + 3 ~ "Gain",
                                copy.number >= ploidy + 3 ~ "Amplification"))
    
    assign("cnvData", cnvData, envir = globalenv())
    
    # filter cnv
    cnvDataFilt <- filter_cnv(myCNVData = cnvData, myCancerGenes = cancer_genes)
    assign("cnvDataFilt", cnvDataFilt, envir = globalenv())
  }
  
  # fusion data (chose either star or arriba or both)
  if(fusion_method == "star"){
    fusPattern = "*star-fusion.fusion_candidates.final"
  } else if(fusion_method == "arriba") {
    fusPattern = "*.arriba.fusions.tsv"
  } else {
    fusPattern = "*star-fusion.fusion_candidates.final|*.arriba.fusions.tsv"
  }
  # read fusion files + filter and merge them
  fusFiles <- list.files(path = file.path(patient_dir, 'gene-fusions'), pattern = fusPattern, recursive = TRUE, full.names = T)
  if(length(fusFiles) >= 1){
    fusFiles <- lapply(fusFiles, filter_fusions, myCancerGenes = cancer_genes)
    fusData <- do.call('rbind', fusFiles)
    fusData <- unique(fusData)
    assign("fusData", fusData, envir = globalenv())
  }
  
  # expression data: TPM (can be only 1 per patient with .genes.results)
  expDat <- list.files(path = file.path(patient_dir, 'gene-expressions'), pattern = "*.genes.results*", recursive = TRUE, full.names = T)
  if(length(expDat) == 1){
    expData <- read.delim(expDat)
    expData.full <- expData %>% 
      mutate(gene_id = str_replace(gene_id, "_PAR_Y_", "_"))  %>%
      separate(gene_id, c("gene_id", "gene_symbol"), sep = "\\_", extra = "merge") %>%
      mutate(gene_id = gsub('[.].*', '', gene_id))  %>%
      unique()
    
    # filter HIST genes
    expData.full <- expData.full[grep('^HIST', expData.full$gene_symbol, invert = T),]
    
    # filter to protein coding genes
    expData.full <- expData.full %>%
      filter(gene_symbol %in% gencode_v27_pc$gene_symbol)
    
    # tpm data
    expData <- expData.full %>% 
      arrange(desc(TPM)) %>% 
      distinct(gene_symbol, .keep_all = TRUE) %>%
      dplyr::rename(!!patient_of_interest := TPM) %>%
      dplyr::select(!!patient_of_interest, gene_id, gene_symbol) %>%
      unique()
    
    # rownames(expData) <- expData$gene_symbol
    assign("expData", expData, envir = globalenv())
    
    # count data
    expData.counts <- expData.full %>%
      filter(expData.full$gene_id  %in% expData$gene_id) %>%
      dplyr::mutate(!!patient_of_interest := expected_count) %>%
      dplyr::select(!!patient_of_interest, gene_id, gene_symbol) %>%
      unique()
    rownames(expData.counts) <- expData.counts$gene_symbol
    assign("expData.counts", expData.counts, envir = globalenv())
  }
}

