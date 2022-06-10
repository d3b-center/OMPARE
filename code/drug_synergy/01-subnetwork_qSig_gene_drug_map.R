# Author: Run Jin
# Obtain Gene Drug Mapping for Subnetworks of Patients of Interest

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(biomaRt)
  library(drugTargetInteractions)
})

# function to get all directories

# Define Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# Parse command line options
option_list <- list(
  make_option(c("--interaction"),type="character",
              help="Gene interaction file for all modules (.tsv) "),
  make_option(c("--enrichment_nes"),type="character",
              help="cluster vs. modules enrichment score file (.tsv) "),
  make_option(c("--cluster"),type="character",
              help="File with cluster assignments for samples (.rds)"),
  make_option(c("--patient_of_interest"),type="character",
              help="Kids_First_Biospecimen_ID of patient of interest"),
  make_option(c("--chemblDb_path"),type="character",
              help="Path to chemblDb sqlite database"),
  make_option(c("--normal_qSig"),type="character",
              help="qSig output from CEMiTools comparing to normal tissues (.txt) "),
  make_option(c("--pediatric_qSig"),type="character",
              help="qSig output from CEMiTools comparing to all pediatric samples (.txt) "),
  make_option(c("--adult_qSig"),type="character",
              help="qSig output from CEMiTools comparing to adult samples (.txt)"),
  make_option(c("--subnetwork"),type="character",
              help="Output file for subnetworks of module of interest (.tsv) "),
  make_option(c("--subnetwork_mapped"),type="character",
              help="Output file for drug-mapped subnetworks (.tsv) "),
  make_option(c("--normal_mapped"),type="character",
              help="Output file for drug-mapped normal qSig subsetted subnetworks (.tsv) "),
  make_option(c("--pediatric_mapped"),type="character",
              help="Output file for drug-mapped pediatric qSig subsetted subnetworks (.tsv) "),
  make_option(c("--adult_mapped"),type="character",
              help="Output file for drug-mapped adult subsetted subnetworks (.tsv) ")
)

opt <- parse_args(OptionParser(option_list = option_list))
patient_of_interest <- opt$patient_of_interest
interaction_df <- readr::read_tsv(opt$interaction)
enrichment_nes_df <- readr::read_tsv(opt$enrichment_nes)
colnames(enrichment_nes_df)[2:ncol(enrichment_nes_df)] <- paste0('cluster_', colnames(enrichment_nes_df)[2:ncol(enrichment_nes_df)])
cluster_df <- readr::read_tsv(opt$cluster)
normal_qSig <- readr::read_tsv(opt$normal_qSig)
pediatric_qSig <- readr::read_tsv(opt$pediatric_qSig)
adult_qSig <- readr::read_tsv(opt$adult_qSig)

#### Run drug annnotation for subnetwork -------------------------------------

# Find the cluster assignment of samples of interest 
cluster_assignment <- cluster_df %>% 
  dplyr::filter(Kids_First_Biospecimen_ID == patient_of_interest) %>%
  pull(cluster_assigned_nb)
cluster_assignment <- paste0("cluster_", cluster_assignment)

# Find the modules positive correlated with the cluster
enrichment_nes_df <- enrichment_nes_df %>%
  dplyr::select(pathway, all_of(cluster_assignment))
colnames(enrichment_nes_df) = c("modules", "cluster_of_interest")

module_selected <- enrichment_nes_df %>% 
  filter(cluster_of_interest > 0) %>%
  pull(modules) %>% unique()

if((length(module_selected) == 0) | (nrow(normal_qSig) == 0 & nrow(pediatric_qSig) == 0 & nrow(adult_qSig) == 0)){
  # do nothing
} else {
  # Use the module selected to subset the gene interaction file for all modules and generate sub-network
  sub_network <- interaction_df %>%
    dplyr::filter(Module %in% module_selected)
  
  readr::write_tsv(sub_network, file = opt$subnetwork)
  
  # Define dataframes to store the results
  qresult2_combined <- data.frame()
  subnetwork_normal_qSig_combined <- data.frame()
  subnetwork_pediatric_qSig_combined <- data.frame()
  subnetwork_adult_qSig_combined <- data.frame()
  
  for(i in 1:length(module_selected)){
    module_each <- module_selected[i]
    sub_network_each <- interaction_df %>%
      dplyr::filter(Module == module_each)
    
    
    # Get the union of all genes in the subnetwork
    gene1_list <- sub_network_each %>% pull(Gene1) %>% unique()
    gene2_list <- sub_network_each %>% pull(Gene2) %>% unique()
    
    all_genes_in_subnetwork <- union(gene1_list, gene2_list) %>% unique()
    
    # use the vector of all genes in the subnetwork as input to getSymEnsUp function.
    chemblDbPath <- opt$chemblDb_path
    resultsPath <- system.file("extdata", "results", package = "drugTargetInteractions")
    config <- genConfig(chemblDbPath = chemblDbPath, resultsPath = resultsPath)
    
    # generate a mapping of uniprot ids and gene symbols only once save as rds
    map_file <- file.path(data_dir, "uniprot_genesymbol_map.rds")
    if(!file.exists(map_file)){
      hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
      annotation.uni <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol"), mart=hsmart)
      saveRDS(annotation.uni, file = map_file)
    } else {
      annotation.uni <- readRDS(map_file)
    }
    
    # filter to genes of interest
    annotation.uni <- annotation.uni %>%
      filter(hgnc_symbol %in% all_genes_in_subnetwork )
    queryBy <- list(molType = "protein", idType = "UniProt_ID", ids = annotation.uni$uniprotswissprot)
    qresult2 <- drugTargetAnnot(queryBy, config=config)
    qresult2 <- qresult2 %>%
      filter(!is.na(First_Approval))
    
    # add gene name back to table
    qresult2 <- qresult2 %>%
      inner_join(annotation.uni, by = c("UniProt_ID" = "uniprotswissprot")) %>%
      mutate(module = module_each)
    
    qresult2_combined <- rbind(qresult2_combined, qresult2)
    #### Intersect subnetwork and qSig output-------------------------------------
    
    # First select drugs that are negatively correlated (WTCS<0) and keep the WTCS scores
    normal_qSig_df <- normal_qSig %>% 
      dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
      dplyr::select(Drug_Name, WTCS)
    normal_qSig_drugs <- normal_qSig_df %>% pull(Drug_Name) %>% unique()
    
    pediatric_qSig_df <- pediatric_qSig %>% 
      dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
      dplyr::select(Drug_Name, WTCS)
    pediatric_qSig_drugs <- pediatric_qSig_df %>% pull(Drug_Name) %>% unique()
    
    adult_qSig_df <- adult_qSig %>% 
      dplyr::filter(WTCS<0) %>% mutate(Drug_Name = toupper(pert)) %>% 
      dplyr::select(Drug_Name, WTCS)
    adult_qSig_drugs <- adult_qSig_df %>% pull(Drug_Name) %>% unique()
    
    # Combine the results with WTCS scores and save the results
    subnetwork_normal_qSig <- qresult2 %>% 
      dplyr::filter(Drug_Name %in% normal_qSig_drugs) %>%
      dplyr::left_join(normal_qSig_df) %>% distinct() %>% 
      dplyr::mutate(module = module_each)
    subnetwork_normal_qSig_combined <- rbind(subnetwork_normal_qSig_combined, subnetwork_normal_qSig)
    
    subnetwork_pediatric_qSig <- qresult2 %>% 
      dplyr::filter(Drug_Name %in% pediatric_qSig_drugs) %>%
      dplyr::left_join(pediatric_qSig_df) %>% distinct() %>% 
      dplyr::mutate(module = module_each)
    subnetwork_pediatric_qSig_combined <- rbind(subnetwork_pediatric_qSig_combined, subnetwork_pediatric_qSig)
    
    subnetwork_adult_qSig <- qresult2 %>% 
      dplyr::filter(Drug_Name %in% adult_qSig_drugs) %>%
      dplyr::left_join(adult_qSig_df) %>% distinct() %>% 
      dplyr::mutate(module = module_each)
    subnetwork_adult_qSig_combined <- rbind(subnetwork_adult_qSig_combined, subnetwork_adult_qSig)
  }
  
  # write out results
  readr::write_tsv(qresult2_combined, file = opt$subnetwork_mapped)
  readr::write_tsv(subnetwork_normal_qSig_combined, file = opt$normal_mapped)
  readr::write_tsv(subnetwork_pediatric_qSig_combined, file = opt$pediatric_mapped)
  readr::write_tsv(subnetwork_adult_qSig_combined, file = opt$adult_mapped)
}
